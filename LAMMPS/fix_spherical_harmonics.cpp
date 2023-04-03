/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_spherical_harmonics.h"

#include "group.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "modify.h"
#include "pair.h"
#include "update.h"
#include "domain.h"

#include <cstring>
#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;
using namespace MathSpecial;

#ifdef DBL_EPSILON
  #define MY_EPSILON (10.0*DBL_EPSILON)
#else
  #define MY_EPSILON (10.0*2.220446049250313e-16)
#endif

#define QEPSILON 1.0e-6

/* ---------------------------------------------------------------------- */

FixSphericalHarmonics::FixSphericalHarmonics(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  ylist(nullptr), ylm_r(nullptr), ylm_i(nullptr), ylm_r_loc(nullptr), ylm_i_loc(nullptr), 
  data_histo(nullptr), data_histo_loc(nullptr)
{

  if (narg < 5 ) error->all(FLERR,"Illegal compute spherical harmonics command");
  
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  jgroup = group->find(arg[3]);
  jgroupbit = group->bitmask[jgroup];
  
  nevery = utils::numeric(FLERR,arg[4],false,lmp);
  nsteps = utils::numeric(FLERR,arg[5],false,lmp);  

  // set default values

  nylist = 5;
  memory->create(ylist,nylist,"sphericalharmonics:ylist");
  ylist[0] = 4;
  ylist[1] = 6;
  ylist[2] = 8;
  ylist[3] = 10;
  ylist[4] = 12;
  ymax = 12;
  nbin_r = 1;
  range_r = 5;

  // process optional args
    
  print = 0;

  int iarg = 6;
  //printf("%d\n", narg);
  while (iarg < narg) {
    if (strcmp(arg[iarg],"lvals") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute spherical harmonics command");
      nylist = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (nylist <= 0)
        error->all(FLERR,"Illegal compute spherical harmonics command");
      memory->destroy(ylist);
      memory->create(ylist,nylist,"sphericalharmonics:ylist");
      iarg += 2;
      if (iarg+nylist > narg)
        error->all(FLERR,"Illegal compute spherical harmonics command");
      ymax = 0;
      for (int il = 0; il < nylist; il++) {
        ylist[il] = utils::numeric(FLERR,arg[iarg+il],false,lmp);
        if (ylist[il] < 0)
          error->all(FLERR,"Illegal compute spherical harmonics command");
        if (ylist[il] > ymax) ymax = ylist[il];
      }
      iarg += nylist;
    } else if (strcmp(arg[iarg],"shells") == 0) {
      if (iarg+3 > narg)
        error->all(FLERR,"Illegal compute spherical harmonics command");
      nbin_r = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      range_r = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      iarg += 3;
    } else if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute spherical harmonics command");
      print = 1;
      if (me == 0) {
        out = fopen(arg[iarg+1],"w");
        if (out == NULL) {
          char str[128];
          sprintf(str,"Cannot open compute spherical harmonics/file file %s",arg[iarg+1]);
          error->one(FLERR,str);
        }
      }
      iarg += 2; 
    } else error->all(FLERR,"Illegal compute spherical harmonics command");
  }
  
  // allocate memory
  size_vector = nylist*nbin_r*(2*ymax+1);
  
  memory->create(data_histo,nbin_r,"sphericalharmonics:data_histo");
  memory->create(data_histo_loc,nbin_r,"sphericalharmonics:data_histo_loc");
  memory->create(data_histo_avg,nbin_r,"sphericalharmonics:data_histo_avg");
  
  memory->create(ylm_r,size_vector,"sphericalharmonics:ylm_r");
  memory->create(ylm_r_loc,size_vector,"sphericalharmonics:ylm_r_loc");
  memory->create(ylm_r_avg,size_vector,"sphericalharmonics:ylm_r_avg");
  memory->create(ylm_i,size_vector,"sphericalharmonics:ylm_i");
  memory->create(ylm_i_loc,size_vector,"sphericalharmonics:ylm_i_loc");
  memory->create(ylm_i_avg,size_vector,"sphericalharmonics:ylm_i_avg");
  
  count = nevery;

}

/* ---------------------------------------------------------------------- */

FixSphericalHarmonics::~FixSphericalHarmonics()
{
  memory->destroy(data_histo);
  memory->destroy(data_histo_loc);
  memory->destroy(data_histo_avg);
  memory->destroy(ylm_r);
  memory->destroy(ylm_r_loc);
  memory->destroy(ylm_r_avg);
  memory->destroy(ylm_i);
  memory->destroy(ylm_i_loc);
  memory->destroy(ylm_i_avg);
  memory->destroy(ylist);
  if (me==0) fclose(out);
}

/* ---------------------------------------------------------------------- */

int FixSphericalHarmonics::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSphericalHarmonics::init()
{

  numavg = 0;

  for (int il = 0; il < nylist*nbin_r*(2*ymax+1); il++) {
    ylm_r_avg[il] = 0.0;
    ylm_i_avg[il] = 0.0;
  }
  
  for (int n=0; n<nbin_r; n++) {
    data_histo_avg[n] = 0;
  }
}

/* ---------------------------------------------------------------------- */

void FixSphericalHarmonics::end_of_step() 
{
  if (count == nevery) {
    numavg ++; 
    // determine com of the colloid
    int n;
    double xcm[3],vcm[3],xhat[3],yhat[3],zhat[3];
    double masstotal = group->mass(igroup);
    group->xcm(igroup,masstotal,xcm);
    group->vcm(igroup,masstotal,vcm);
    
    
    //change coordinate system to \hat{z}=\vec{v}_{colloid}
    double zmag = sqrt(vcm[0]*vcm[0] + vcm[1]*vcm[1] + vcm[2]*vcm[2]);
    zhat[0] = vcm[0]/zmag;
    zhat[1] = vcm[1]/zmag;
    zhat[2] = vcm[2]/zmag;
    
    xhat[0] = 0;
    xhat[1] = -zhat[2];
    xhat[2] = zhat[1];
    double xmag = sqrt(xhat[0]*xhat[0] + xhat[1]*xhat[1] + xhat[2]*xhat[2]);
    xhat[0] /= xmag;
    xhat[1] /= xmag;
    xhat[2] /= xmag;
    
    yhat[0] = zhat[1]*zhat[1]+zhat[2]*zhat[2];
    yhat[1] = -zhat[0]*zhat[1];
    yhat[2] = -zhat[0]*zhat[2];
    double ymag = sqrt(yhat[0]*yhat[0] + yhat[1]*yhat[1] + yhat[2]*yhat[2]);
    yhat[0] /= ymag;
    yhat[1] /= ymag;
    yhat[2] /= ymag;
    
    //initialize ylm list for each shell    
    for (int il = 0; il < nylist*nbin_r*(2*ymax+1); il++) {
      ylm_r_loc[il] = 0.0;
      ylm_i_loc[il] = 0.0;
      ylm_r[il] = 0.0;
      ylm_i[il] = 0.0;
    }
      //int l = ylist[il];
      //for (int m = 0; m < 2*l+1; m++) {
        //ylm_r_loc[il][m] = 0.0;
        //ylm_i_loc[il][m] = 0.0;
        //ylm_r[il][m] = 0.0;
        //ylm_i[il][m] = 0.0;
      //}
    //}
    
    for (int n=0; n<nbin_r; n++) {
      data_histo[n] = 0;
      data_histo_loc[n] = 0;
    }

    double delx, dely, delz, xnew, ynew, znew, vxnew, vynew, vznew;
    int i;
    double **x = atom->x;
    int nlocal = atom->nlocal;
    int * mask = atom->mask;
    
    //double check = 0.0;
    //double check2 = 0.0;
    
    for (i=0; i<nlocal; i++) {
      if(mask[i] & jgroupbit) {
        delx = x[i][0] - xcm[0];
        dely = x[i][1] -  xcm[1];
        delz =  x[i][2] - xcm[2];  
        
        domain->minimum_image(delx,dely,delz);
        
        xnew = delx*xhat[0] + dely*xhat[1] + delz*xhat[2];        
        ynew = delx*yhat[0] + dely*yhat[1] + delz*yhat[2];        
        znew = delx*zhat[0] + dely*zhat[1] + delz*zhat[2];
        
        double data_r = sqrt(xnew*xnew + ynew*ynew + znew*znew);
        int bin_r = data_r*((double) nbin_r)/(range_r);
        if (bin_r >= 0 && bin_r <nbin_r) {
          //printf("%f %d\n",data_r,bin_r);
          data_histo_loc[bin_r]++;
          
          double costheta = znew / data_r;
          double expphi_r = xnew;
          double expphi_i = ynew;
          //check += (xnew*xnew-ynew*ynew)/(data_r*data_r);
          //check2 += (-2*ynew*xnew)/(data_r*data_r);
          double rxymag = sqrt(expphi_r*expphi_r+expphi_i*expphi_i);
          if (rxymag <= MY_EPSILON) {
            expphi_r = 1.0;
            expphi_i = 0.0;
            printf("Something happened \n");
          } else {
            double rxymaginv = 1.0/rxymag;
            expphi_r *= rxymaginv;
            expphi_i *= rxymaginv;
          }
          
          for (int il = 0; il < nylist; il++) {
            int l = ylist[il];
            int nmlist = 2*ymax+1;
            // calculate spherical harmonics
            // Ylm, -l <= m <= l
            // sign convention: sign(Yll(0,0)) = (-1)^l
            
            ylm_r_loc[bin_r*nylist*nmlist+il*nmlist+l] += polar_prefactor(l, 0, costheta);
            double expphim_r = expphi_r;
            double expphim_i = expphi_i;
            for (int m = 1; m <= +l; m++) {

              double prefactor = polar_prefactor(l, m, costheta);
              double part_ylm_r = prefactor * expphim_r;
              double part_ylm_i = prefactor * expphim_i;
              ylm_r_loc[bin_r*nylist*nmlist+il*nmlist+m+l] += part_ylm_r;
              ylm_i_loc[bin_r*nylist*nmlist+il*nmlist+m+l] += part_ylm_i;
              if (m & 1) {
                ylm_r_loc[bin_r*nylist*nmlist+il*nmlist-m+l] -= part_ylm_r;
                ylm_i_loc[bin_r*nylist*nmlist+il*nmlist-m+l] += part_ylm_i;
              } else {
                ylm_r_loc[bin_r*nylist*nmlist+il*nmlist-m+l] += part_ylm_r;
                ylm_i_loc[bin_r*nylist*nmlist+il*nmlist-m+l] -= part_ylm_i;
              }
              double tmp_r = expphim_r*expphi_r - expphim_i*expphi_i;
              double tmp_i = expphim_r*expphi_i + expphim_i*expphi_r;
              expphim_r = tmp_r;
              expphim_i = tmp_i;
            }            
          }
        }
      }
    }
    
    //printf("%f\n",ylm_r_loc[0]);
    
    MPI_Allreduce(data_histo_loc, data_histo, nbin_r, MPI_DOUBLE, MPI_SUM, world);
    MPI_Allreduce(ylm_r_loc, ylm_r, nylist*nbin_r*(2*ymax+1), MPI_DOUBLE, MPI_SUM, world);
    MPI_Allreduce(ylm_i_loc, ylm_i, nylist*nbin_r*(2*ymax+1), MPI_DOUBLE, MPI_SUM, world);
    
    int nmlist = 2*ymax+1;
    for (int r = 0; r < nbin_r; r++) {
      data_histo_avg[r] += data_histo[r];
      for (int il = 0; il < nylist; il++) {
        int l = ylist[il];
        for (int m = 0; m < 2*l+1; m++) {
          ylm_r_avg[r*nylist*nmlist+il*nmlist+m] += ylm_r[r*nylist*nmlist+il*nmlist+m];
          ylm_i_avg[r*nylist*nmlist+il*nmlist+m] += ylm_i[r*nylist*nmlist+il*nmlist+m];
        }
      }
    }
    //printf("%f\n",ylm_r_loc[0]);

    if (print && (me == 0) && (update->ntimestep == int(nsteps-nevery+1)) ) {
      //int nmlist = 2*ymax+1;
      //fprintf(out,"t=%d\n",update->ntimestep);
      for (int r = 0; r < nbin_r; r++) {
          fprintf(out,"R=%f %f\n",r/((double) nbin_r)*(range_r),data_histo_avg[r]/numavg);
          for (int il = 0; il < nylist; il++) {
            int l = ylist[il];
            for (int m = 0; m < 2*l+1; m++) {
              fprintf(out,"%d %d %f %f\n",l,m-l,ylm_r_avg[r*nylist*nmlist+il*nmlist+m]/numavg,ylm_i_avg[r*nylist*nmlist+il*nmlist+m]/numavg);
            }
            fprintf(out,"\n");
          }
          fprintf(out,"\n");
      }
      //fprintf(out,"\n");
    }
    count = 0;
  } 
  count ++;
}

/* ----------------------------------------------------------------------
   polar prefactor for spherical harmonic Y_l^m, where
   Y_l^m (theta, phi) = prefactor(l, m, cos(theta)) * exp(i*m*phi)
------------------------------------------------------------------------- */

double FixSphericalHarmonics::polar_prefactor(int l, int m, double costheta)
{
  const int mabs = abs(m);

  double prefactor = 1.0;
  for (int i=l-mabs+1; i < l+mabs+1; ++i)
    prefactor *= static_cast<double>(i);

  prefactor = sqrt(static_cast<double>(2*l+1)/(MY_4PI*prefactor))
    * associated_legendre(l,mabs,costheta);

  if ((m < 0) && (m % 2)) prefactor = -prefactor;

  return prefactor;
}

/* ----------------------------------------------------------------------
   associated legendre polynomial
   sign convention: P(l,l) = (2l-1)!!(-sqrt(1-x^2))^l
------------------------------------------------------------------------- */

double FixSphericalHarmonics::associated_legendre(int l, int m, double x)
{
  if (l < m) return 0.0;

  double p(1.0), pm1(0.0), pm2(0.0);

  if (m != 0) {
    const double msqx = -sqrt(1.0-x*x);
    for (int i=1; i < m+1; ++i)
      p *= static_cast<double>(2*i-1) * msqx;
  }

  for (int i=m+1; i < l+1; ++i) {
    pm2 = pm1;
    pm1 = p;
    p = (static_cast<double>(2*i-1)*x*pm1
         - static_cast<double>(i+m-1)*pm2) / static_cast<double>(i-m);
  }

  return p;
}
