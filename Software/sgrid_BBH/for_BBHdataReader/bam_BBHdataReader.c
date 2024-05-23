/* bam_BBHdataReader.c */
/* Wolfgang Tichy 4/2012 */

#include "bam.h"
#include "BBHdataReader.h"

void bam_BBHdataReader() 
{
  if(!Getv("physics", "BBHdataReader")) return;
  /* functions */
  AddFun(INITIALDATA, BBHdataReader, "Import BBHdata initial data");
  AddFun(PRE_GRID, BBHdataPars, "set some pars important for making the grid");

  /* variables */

  /* parameters */
  printf("Adding BBHdataReader: read data from sgrid BBHdata\n");

  AddPar("bhmass1", "1.0", "mass of black hole 1");
  AddPar("bhx1", "0.0", "origin_x of black hole 1");
  AddPar("bhy1", "0.0", "origin_y of black hole 1");
  AddPar("bhz1", "0.0", "origin_z of black hole 1");
  AddPar("bhpx1", "0.0", "momentum_x of black hole 1");
  AddPar("bhpy1", "0.0", "momentum_y of black hole 1");
  AddPar("bhpz1", "0.0", "momentum_z of black hole 1");
  AddPar("bhsx1", "0.0", "spin_x of black hole 1");
  AddPar("bhsy1", "0.0", "spin_y of black hole 1");
  AddPar("bhsz1", "0.0", "spin_z of black hole 1");

  AddPar("bhmass2", "1.0", "mass of black hole 2");
  AddPar("bhx2", "0.0", "origin_x of black hole 2");
  AddPar("bhy2", "0.0", "origin_y of black hole 2");
  AddPar("bhz2", "0.0", "origin_z of black hole 2");
  AddPar("bhpx2", "0.0", "momentum_x of black hole 2");
  AddPar("bhpy2", "0.0", "momentum_y of black hole 2");
  AddPar("bhpz2", "0.0", "momentum_z of black hole 2");
  AddPar("bhsx2", "0.0", "spin_x of black hole 2");
  AddPar("bhsy2", "0.0", "spin_y of black hole 2");
  AddPar("bhsz2", "0.0", "spin_z of black hole 2");
  AddPar("punctures_lapse", "none", "lapse [none,rtoN_atPunc]");
  AddPar("punctures_lapse_rPower_atPunc", "4", "exponent N for rtoN_atPunc lapse");

  /* the sgrid data could be orgsanized like this:
     datadir
       sgrid
       name
       name2
  */
  AddPar("BBHdataReader_sgrid_exe", "./sgrid", "location of sgrid executable");
  AddPar("BBHdataReader_sgrid_datadir", "", "location of sgrid outdir with data");
  AddPar("BBHdataReader_keep_sgrid_output", "no",
         "whether we keep the output sgrid creates while running [yes,no]");
  AddPar("BBHdataReader_use_interpolator", "yes", 
         "whether we use BBHdataReader_interpolator");
  AddPar("BBHdataReader_IDfiles_dir", "",
         "Directory with ID .dat-files already interpolated to "
         "our current bam grid (needed if BBHdataReader_use_interpolator = no)");
  AddPar("BBHdataReader_sgrid_x_CM",  "0", "sgrid's CM of system");
  AddPar("BBHdataReader_m1", "0", "rest mass of star1");
  AddPar("BBHdataReader_m2", "0", "rest mass of star2");
  AddPar("BBHdataReader_xc1", "0", "pos. of max density in star1");
  AddPar("BBHdataReader_xc2", "0", "pos. of max density in star2");

  /* Robin boundary */
/*
  AddPar("BBHdataReader_boundary", "robin", "boundary condition");
  AddPar("BBHdataReader_robin_falloff", "1", "If zero, use Dirichlet boundary condition");

  if (Getv("BBHdataReader_boundary", "robin"))
  {
    AddVar("robinindex", "", "index for Robin boundary");
    AddVar("robinflag",  "", "temporary flag for Robin boundary");
  } 
*/
}
