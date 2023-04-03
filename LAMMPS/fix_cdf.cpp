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

#include "fix_cdf.h"
#include "string.h"
#include "update.h"
#include "group.h"
#include "error.h"
#include "force.h"
#include "pair.h"
#include "memory.h"
#include "modify.h"
#include "atom.h"
#include "domain.h"
#include "input.h"
#include "comm.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixCDF::FixCDF(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  
  if (narg < 9) error->all(FLERR,"Illegal compute cdf command");
  
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
  
  nevery = utils::numeric(FLERR,arg[3],false,lmp);
  
  nbin_z = utils::numeric(FLERR,arg[4],false,lmp);
  range_z = utils::numeric(FLERR,arg[5],false,lmp);
  nbin_r = utils::numeric(FLERR,arg[6],false,lmp);
  range_r = utils::numeric(FLERR,arg[7],false,lmp);
  
  jgroup = group->find(arg[8]);
  jgroupbit = group->bitmask[jgroup];
  
  print = 0;
  
  //printf("out: %s\n",arg[9]);
  
  if (strcmp(arg[9],"file") == 0) {
    if (11 > narg) error->all(FLERR,"Illegal compute cdf/file command");
    print = 1;
    //printf("me: %d\n",me);
    if (me == 0) {
      //printf("Fix cdf, option file\n");
      out = fopen(arg[10],"w");
      if (out == NULL) {
        char str[128];
        sprintf(str,"Cannot open compute cdf/file file %s",arg[10]);
        error->one(FLERR,str);
      }
    } 
  }
  
  // allocate memory
  size_vector = nbin_z*nbin_r;
  memory->create(data_histo,size_vector,"cdf/data_histo");
  memory->create(data_histo_loc,size_vector,"cdf/data_histo_loc");
  memory->create(data_velz,size_vector,"cdf/data_velz");
  memory->create(data_velr,size_vector,"cdf/data_velr");
  memory->create(data_vel_locz,size_vector,"cdf/data_vel_locz");
  memory->create(data_vel_locr,size_vector,"cdf/data_vel_locr");
  
  memory->create(data_temp,size_vector,"cdf/data_temp");    
  memory->create(data_temp_loc,size_vector,"cdf/data_temp_loc");
  
  count = nevery;
  
}

/* ---------------------------------------------------------------------- */

FixCDF::~FixCDF()
{
  memory->destroy(data_histo);
  memory->destroy(data_histo_loc);
  memory->destroy(data_velz);
  memory->destroy(data_velr);
  memory->destroy(data_vel_locz);
  memory->destroy(data_vel_locr);
  memory->destroy(data_temp);
  memory->destroy(data_temp_loc);
  if (me==0) fclose(out);
}

/* ---------------------------------------------------------------------- */

int FixCDF::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCDF::init()
{

}

/* ---------------------------------------------------------------------- */

void FixCDF::end_of_step() 
{
  if (count == nevery) {
    // determine com of the colloid
    int n;
    double xcm[3],vcm[3],xhat[3],yhat[3],zhat[3];
    double masstotal = group->mass(igroup);
    group->xcm(igroup,masstotal,xcm);
    group->vcm(igroup,masstotal,vcm);
    
    //printf("%f %f %f %f %f %f\n", xcm[0], xcm[1], xcm[2], vcm[0], vcm[1], vcm[2]);
    
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
    
    //printf("%f %f %f\n", sqrt(xhat[0]*xhat[0] + xhat[1]*xhat[1] + xhat[2]*xhat[2]),sqrt(yhat[0]*yhat[0] + yhat[1]*yhat[1] + yhat[2]*yhat[2]),sqrt(zhat[0]*zhat[0] + zhat[1]*zhat[1] + zhat[2]*zhat[2]));
        
    for (n=0; n<nbin_z*nbin_r; n++) {
      data_histo[n] = 0;
      data_histo_loc[n] = 0;
      data_velz[n] = 0.0;
      data_vel_locz[n] = 0.0;
      data_velr[n] = 0.0;
      data_vel_locr[n] = 0.0;
      
      data_temp[n] = 0.0;
      data_temp_loc[n] = 0.0;
    }
    
    double delx, dely, delz, xnew, ynew, znew, vxnew, vynew, vznew;
    int i;
    double **x = atom->x;
    double **v = atom->v;
    int nlocal = atom->nlocal;
    int * mask = atom->mask;
    
    for (i=0; i<nlocal; i++) {
      if(mask[i] & jgroupbit) {
        delx = x[i][0] - xcm[0];
        dely = x[i][1] -  xcm[1];
        delz =  x[i][2] - xcm[2];  
        
        domain->minimum_image(delx,dely,delz);
        
        xnew = delx*xhat[0] + dely*xhat[1] + delz*xhat[2];        
        ynew = delx*yhat[0] + dely*yhat[1] + delz*yhat[2];        
        znew = delx*zhat[0] + dely*zhat[1] + delz*zhat[2];
        
        vxnew = v[i][0]*xhat[0] + v[i][1]*xhat[1] + v[i][2]*xhat[2];        
        vynew = v[i][0]*yhat[0] + v[i][1]*yhat[1] + v[i][2]*yhat[2];        
        vznew = v[i][0]*zhat[0] + v[i][1]*zhat[1] + v[i][2]*zhat[2];
            
        double data_z = znew + range_z;
        double data_r = sqrt(xnew*xnew + ynew*ynew);
        int bin_z=-1;
        if (data_z > 0.0)
          bin_z = data_z*((double) nbin_z)/(2.0*range_z);
        int bin_r = data_r*((double) nbin_r)/(range_r);
        if (bin_z >= 0 && bin_z <nbin_z && bin_r >= 0 && bin_r <nbin_r) {
          data_histo_loc[bin_z*nbin_r + bin_r]++;
          data_vel_locz[bin_z*nbin_r + bin_r]+=vznew;
          double sign = vxnew*xnew+vynew*ynew;
          double dx = xnew*xnew + ynew*ynew;
          data_vel_locr[bin_z*nbin_r + bin_r]+=sign/sqrt(dx);
	  data_temp_loc[bin_z*nbin_r + bin_r]+=(vxnew*vxnew+vynew*vynew+vznew*vznew)/3.0;
	  
	  //printf("%f %f\n",(vxnew*vxnew+vynew*vynew+vznew*vznew)/3.0,(v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2])/3.0);
        }
      }
    }
    
    MPI_Allreduce(data_histo_loc, data_histo, nbin_z*nbin_r, MPI_DOUBLE, MPI_SUM, world);
    MPI_Allreduce(data_vel_locz, data_velz, nbin_z*nbin_r, MPI_DOUBLE, MPI_SUM, world);
    MPI_Allreduce(data_vel_locr, data_velr, nbin_z*nbin_r, MPI_DOUBLE, MPI_SUM, world);    
    MPI_Allreduce(data_temp_loc, data_temp, nbin_z*nbin_r, MPI_DOUBLE, MPI_SUM, world);
    
    int r,zloc;
    for (zloc=0; zloc<nbin_z; zloc++) {
      for (r=0; r<nbin_r; r++) {
        if ( data_histo[zloc*nbin_r + r] != 0 ) {
        
          data_velz[zloc*nbin_r + r] /= data_histo[zloc*nbin_r + r];
          data_velr[zloc*nbin_r + r] /= data_histo[zloc*nbin_r + r];	  
	  data_temp[zloc*nbin_r + r] /= data_histo[zloc*nbin_r + r];
        }
      }
    }

    if (print && (me == 0) ) {
      fprintf(out,"t=%d\n",update->ntimestep);
      for (zloc=0; zloc<nbin_z; zloc++) {
        for (r=0; r<nbin_r; r++) {
          fprintf(out,"%f %f %f %f %f %f\n",zloc/((double) nbin_z)*(2.0*range_z)-range_z,r/((double) nbin_r)*(range_r),data_histo[zloc*nbin_r + r],data_velz[zloc*nbin_r + r],data_velr[zloc*nbin_r + r],data_temp[zloc*nbin_r + r]);
        }
        fprintf(out,"\n");
      }
      fprintf(out,"\n\n");
    }
    count = 0;
  } 
  count ++;
}
