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

FixStyle(velocitycorrelation,FixVelocityCorrelation)

#else

#ifndef LMP_FIX_VELOCITYCORRELATION
#define LMP_FIX_VELOCITYCORRELATION

#include "fix.h"

namespace LAMMPS_NS {

class FixVelocityCorrelation : public Fix {
 public:
  FixVelocityCorrelation(class LAMMPS *, int, char **);
  ~FixVelocityCorrelation();
  int setmask();
  void init();
  void end_of_step();
  int *ylist;
  int nylist;

 protected:
  int me,nprocs;
  int nevery,count,nsteps,numavg;
  int nbin_r;
  double range_r;
  int ymax;
  int style;
  
  double polar_prefactor(int, int, double);
  double associated_legendre(int, int, double);
  
  double * data_histo;
  double * data_histo_loc;
  double * data_histo_avg;
  double * data_velcorr;
  double * data_velcorr_loc;
  double * data_velcorr_avg;
  double * data_velcorrpar;
  double * data_velcorrpar_loc;
  double * data_velcorrpar_avg;
  double * data_velcorrperp;
  double * data_velcorrperp_loc;
  double * data_velcorrperp_avg;
  double * ylm_r;
  double * ylm_r_loc;
  double * ylm_r_avg;
  double * ylm_i;
  double * ylm_i_loc;
  double * ylm_i_avg;
      
  int jgroup,jgroupbit;
  
  int print;
  FILE * out;
  FILE * out2;
};

}

#endif
#endif
