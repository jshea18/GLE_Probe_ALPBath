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

#include "fix_movereg.h"
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

FixMOVEREG::FixMOVEREG(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  
  if (narg < 6) error->all(FLERR,"Illegal fix move region command");
  
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
  
  nevery = utils::numeric(FLERR,arg[3],false,lmp);
  
  neq = utils::numeric(FLERR,arg[4],false,lmp);
  
  range = utils::numeric(FLERR,arg[5],false,lmp);
  
  jgroup = group->find(arg[6]);
  jgroupbit = group->bitmask[jgroup];
  
  print = 0;
  
  //printf("nevery: %d\n",arg[3]);
  
  if (strcmp(arg[7],"file") == 0) {
    if (9 > narg) error->all(FLERR,"Illegal move region/file command");
    print = 1;
    //printf("me: %d\n",me);
    if (me == 0) {
      //printf("Fix cdf, option file\n");
      out = fopen(arg[8],"w");
      if (out == NULL) {
        char str[128];
        sprintf(str,"Cannot open move region/file file %s",arg[8]);
        error->one(FLERR,str);
      }
    } 
  }
  
  // allocate memory
  memory->create(data_xcol,1,"movereg/data_xcol");
  memory->create(data_ycol,1,"movereg/data_ycol");
  memory->create(data_zcol,1,"movereg/data_zcol");
  memory->create(data_velx,1,"movereg/data_velx");
  memory->create(data_vely,1,"movereg/data_vely");
  memory->create(data_velz,1,"movereg/data_velz");
  memory->create(data_numpart,1,"movereg/data_numpart");
  memory->create(data_velx_global,1,"movereg/data_velx_global");
  memory->create(data_vely_global,1,"movereg/data_vely_global");
  memory->create(data_velz_global,1,"movereg/data_velz_global");
  memory->create(data_numpart_global,1,"movereg/data_numpart_global");
  memory->create(data_average,1,"movereg/data_average");
  
  count = nevery;
  
}

/* ---------------------------------------------------------------------- */

FixMOVEREG::~FixMOVEREG()
{
  memory->destroy(data_xcol);
  memory->destroy(data_ycol);
  memory->destroy(data_zcol);
  memory->destroy(data_velx);
  memory->destroy(data_vely);
  memory->destroy(data_velz);
  memory->destroy(data_numpart);
  memory->destroy(data_velx_global);
  memory->destroy(data_vely_global);
  memory->destroy(data_velz_global);
  memory->destroy(data_numpart_global);
  memory->destroy(data_average);

  if (me==0) fclose(out);
}

/* ---------------------------------------------------------------------- */

int FixMOVEREG::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMOVEREG::init()
{

}

/* ---------------------------------------------------------------------- */

void FixMOVEREG::end_of_step() 
{
  if (count == nevery) {
    // determine com of the colloid
    int n;
    double xcm[3],vcm[3];
    double masstotal = group->mass(igroup);
    group->xcm(igroup,masstotal,xcm);

    data_xcol[0] = xcm[0];
    data_ycol[0] = xcm[1];
    data_zcol[0] = xcm[2];
    
    data_velx[0] = 0.0;
    data_vely[0] = 0.0;
    data_velz[0] = 0.0;
    data_average[0] = 0.0;
    data_numpart[0] = 0;

    double max_dist = 0;
    
    double delx, dely, delz, max_x, max_y, max_z;
    int i;
    double **x = atom->x;
    double **v = atom->v;
    int nlocal = atom->nlocal;
    int * mask = atom->mask;
    for (i=0; i<nlocal; i++) {
      if(mask[i] & jgroupbit) {
        delx = x[i][0] - xcm[0];
        dely = x[i][1] -  xcm[1];
        delz = x[i][2] - xcm[2];

        //printf("before min im: %f %f %f\n",delx, dely, delz);
        domain->minimum_image(delx,dely,delz);
        //printf("after min im: %f %f %f\n",delx, dely, delz);
    
        double dist = sqrt(delx*delx + dely*dely + delz*delz);
        if (dist < range) {
          //printf("Position solv: %f %f %f\n", x[i][0],x[i][1],x[i][2]);
          //printf("Position coll: %f %f %f\n", xcm[0],xcm[1],xcm[2]);
          //printf("dist: %f\n",dist);
          data_velx[0] += v[i][0];
          data_vely[0] += v[i][1];
          data_velz[0] += v[i][2];
          data_numpart[0]++;

/*          if (dist > max_dist) {
            max_dist = dist;
            max_x = x[i][0];
            max_y = x[i][1];
            max_z = x[i][2];
          }
*/
        }
      }
    }

    //printf("Max dist: %f, max x: %f, max y: %f, max z: %f, col_x: %f, col_y: %f, col_z: %f \n", max_dist, max_x, max_y, max_z, xcm[0], xcm[1], xcm[2]);

    //printf("Number of atoms: %d\n", num_atoms);  
    
    
    //printf("test1\n");

    MPI_Allreduce(data_velx, data_velx_global, 1, MPI_DOUBLE, MPI_SUM, world);
    MPI_Allreduce(data_vely, data_vely_global, 1, MPI_DOUBLE, MPI_SUM, world);
    MPI_Allreduce(data_velz, data_velz_global, 1, MPI_DOUBLE, MPI_SUM, world);

    MPI_Allreduce(data_numpart, data_numpart_global, 1, MPI_INT, MPI_SUM, world);
    
    //MPI_Allreduce(data_histo_loc, data_histo, nbin_x*nbin_r, MPI_DOUBLE, MPI_SUM, world);
    //MPI_Allreduce(data_vel_locx, data_velx, nbin_x*nbin_r, MPI_DOUBLE, MPI_SUM, world);
    //MPI_Allreduce(data_vel_locy, data_vely, nbin_x*nbin_r, MPI_DOUBLE, MPI_SUM, world);
    //MPI_Allreduce(data_vel_locz, data_velz, nbin_x*nbin_r, MPI_DOUBLE, MPI_SUM, world);
    
    //MPI_Allreduce(data_press_locx, data_pressx, nbin_x*nbin_r, MPI_DOUBLE, MPI_SUM, world);
    //MPI_Allreduce(data_press_locyz, data_pressyz, nbin_x*nbin_r, MPI_DOUBLE, MPI_SUM, world);
    //MPI_Allreduce(data_temp_loc, data_temp, nbin_x*nbin_r, MPI_DOUBLE, MPI_SUM, world);
    
    //printf("test3\n");

    if (data_numpart_global[0] != 0) {
      data_velx_global[0] /= data_numpart_global[0];
      data_vely_global[0] /= data_numpart_global[0];
      data_velz_global[0] /= data_numpart_global[0];
      data_average[0] = data_velx_global[0]*data_velx_global[0] + data_vely_global[0]*data_vely_global[0] + data_velz_global[0]*data_velz_global[0];
    }

    for (i=0; i<nlocal; i++) {
      if(mask[i] & groupbit) {
        v[i][0]=data_velx_global[0];
        v[i][1]=data_vely_global[0];
        v[i][2]=data_velz_global[0];
      }
    }

    if (print && (me == 0) ) {
      fprintf(out,"%d %d %f %f %f %f %f %f %f\n",update->ntimestep,data_numpart_global[0],data_velx_global[0],data_vely_global[0],data_velz_global[0],data_average[0],data_xcol[0],data_ycol[0],data_zcol[0]);
    }
    count = 0;
  } 
  count ++;
}
