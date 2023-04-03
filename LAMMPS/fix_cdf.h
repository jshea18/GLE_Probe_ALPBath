/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(cdf,FixCDF)

#else

#ifndef LMP_FIX_CDF
#define LMP_FIX_CDF

#include "fix.h"

namespace LAMMPS_NS {

class FixCDF : public Fix {
 public:
  FixCDF(class LAMMPS *, int, char **);
  ~FixCDF();
  int setmask();
  void init();
  void end_of_step();

 private:
  int me,nprocs;
  int nevery,count;
  int nbin_z;
  double range_z;
  int nbin_r;
  double range_r;
  
  double * data_histo;
  double * data_histo_loc;
  double * data_velz;
  double * data_vel_locz;
    //double * data_vely;
  //double * data_vel_locy;
    double * data_velr;
  double * data_vel_locr;
    
            double * data_temp;
    double * data_temp_loc;
  
  int jgroup,jgroupbit;
  
  int print;
  FILE * out;
};

}

#endif
#endif
