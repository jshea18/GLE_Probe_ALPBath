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

FixStyle(movereg,FixMOVEREG)

#else

#ifndef LMP_FIX_MOVEREG
#define LMP_FIX_MOVEREG

#include "fix.h"

namespace LAMMPS_NS {

class FixMOVEREG : public Fix {
 public:
  FixMOVEREG(class LAMMPS *, int, char **);
  ~FixMOVEREG();
  int setmask();
  void init();
  void end_of_step();

 private:
  int me,nprocs;
  int nevery,count,neq;
  double range;

  double * data_xcol;
  double * data_ycol;
  double * data_zcol;  
  double * data_velx;
  double * data_velx_global;
  double * data_vely;
  double * data_vely_global;
  double * data_velz;
  double * data_velz_global;
  int * data_numpart;
  int * data_numpart_global;
  double * data_average;
  
  int jgroup,jgroupbit;
  
  int print;
  FILE * out;
};

}

#endif
#endif
