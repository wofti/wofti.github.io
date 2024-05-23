/* BBHdataReader.c */
/* Wolfgang Tichy 4/2012 */


#include "bam.h"
#include "BBHdataReader.h"

#define NDATAMAX 26
#define STRLEN 16384

/* Variables shared by BBHdataReader subroutines */
double bh_pos[2][3];


/* position filepointer after the string str */
int BBH_position_fileptr_after_str(FILE *in, char *str)
{
  char line[STRLEN];
  
  while(fgets(line, STRLEN-1, in)!=NULL)
  {
    if(strstr(line, str)!=NULL) return 1; //break;
  }
  return EOF;
}


/* Compute BBHdataReader data */
int BBHdataReader(tL *level)
{
  tG *grid = level->grid;
  int pr=0;
  int i, j, n;
  int ndata = NDATAMAX;
  int Knonzero;
  double *x = level->v[Ind("x")];
  double *y = level->v[Ind("y")];
  double *z = level->v[Ind("z")];
  double *gxx, *gxy, *gxz, *gyy, *gyz, *gzz;
  double *Kxx, *Kxy, *Kxz, *Kyy, *Kyz, *Kzz;
  double *psi, *alpha, *betax, *betay, *betaz;
  double *dpsiopsi1, *dpsiopsi2, *dpsiopsi3;
  double *ddpsiopsi11, *ddpsiopsi12, *ddpsiopsi13;
  double *ddpsiopsi22, *ddpsiopsi23, *ddpsiopsi33;
  double r, rmax;
  int rank = bampi_rank();
  char *outdir = Gets("outdir");
  FILE *fp1, *fp2;
  char str[STRLEN];
  char gridfile[STRLEN], call_interpolator[STRLEN], IDfile[STRLEN];
  char sgridoutdir[STRLEN], sgridoutdir_previous[STRLEN];
  char sgridcheckpoint_indir[STRLEN], IDfile_new[STRLEN];
  double sgrid_x_CM = Getd("BBHdataReader_sgrid_x_CM");
  int use_interpolator, ret;

  /* say where we are */
  prdivider(0);
  printf("BBHdataReader: ");
  printf("(grid=%p  level=%p  level->l=%d)\n", grid, level, level->l);

  /* some info */
  printf("  sgrid_x_CM = %g\n", sgrid_x_CM);

  /* initialize file names */
  sprintf(gridfile, "%s/grid_level_%d_proc_%d.dat", outdir, level->l, rank);
  sprintf(sgridoutdir, "%s/sgrid_proc_%d", outdir, rank);
  sprintf(sgridoutdir_previous, "%s/sgrid_proc_%d_previous", outdir, rank);
  snprintf(sgridcheckpoint_indir, STRLEN-1, "%s",
           Gets("BBHdataReader_sgrid_datadir"));
  
  /* enable storage, and get pointers to ADM vars */
  i = Ind("psi"); enablevar(level, i);
  psi = level->v[i++];

  i = Ind("dpsiopsix"); enablevar(level, i);
  dpsiopsi1 = level->v[i+0];
  dpsiopsi2 = level->v[i+1];
  dpsiopsi3 = level->v[i+2];

  i = Ind("ddpsiopsixx"); enablevar(level, i);
  ddpsiopsi11 = level->v[i+0];
  ddpsiopsi12 = level->v[i+1];
  ddpsiopsi13 = level->v[i+2];
  ddpsiopsi22 = level->v[i+3];
  ddpsiopsi23 = level->v[i+4];
  ddpsiopsi33 = level->v[i+5];

  i = Ind("alpha"); enablevar(level, i);
  alpha = level->v[i++];

  i = Ind("betax"); enablevar(level, i);
  betax = level->v[i++]; betay = level->v[i++]; betaz = level->v[i++];

  i = Ind("gxx"); enablevar(level, i);
  gxx = level->v[i++]; gxy = level->v[i++]; gxz = level->v[i++];
  gyy = level->v[i++]; gyz = level->v[i++]; gzz = level->v[i++];

  i = Ind("Kxx"); enablevar(level, i);
  Kxx = level->v[i++]; Kxy = level->v[i++]; Kxz = level->v[i++];
  Kyy = level->v[i++]; Kyz = level->v[i++]; Kzz = level->v[i++];

  /* read parameters */
  //ReadBBHdataParameters(&Knonzero);


  /* dump the grid points covered by this processor on this level 
     in gridfile */
  fp1 = fopen(gridfile, "wb");
  if(fp1==NULL) errorexits("could not open %s", gridfile);
  fprintf(fp1, "%s", "# this pointsfile contains the (x+sgrid_x_CM, y, z) of the bam grid points\n");
  fprintf(fp1, "%s", "$BEGIN_data:\n");
  forallpoints(level, i)
  {
    double xs;

    /* fprintf(fp1, "%e %e %e\n",x[i], y[i], z[i]); */
    xs = x[i] + sgrid_x_CM; /* shift */
    fwrite(&xs, sizeof(double), 1, fp1);
    fwrite(y+i, sizeof(double), 1, fp1);
    fwrite(z+i, sizeof(double), 1, fp1);
  }
  fclose(fp1);
  fflush(stdout);

  /* decide if we call the interpolator, or if we just read data */
  if(!Getv("BBHdataReader_use_interpolator", "yes"))
  {
    /* If you want to use previously generated ID use these two lines
       Set BBHdataReader_IDfiles_dir to the dir where you put the "ID*.dat"
       files generated in a previous run */
    sprintf(IDfile, "%s/ID_level_%d_proc_%d.dat",
            Gets("BBHdataReader_IDfiles_dir"), level->l, rank);
    printf("Looking for file %s\n", IDfile);

    /* call interpolator if file cannot be opened */
    fp2 = fopen(IDfile, "rb");
    if(fp2==NULL)
    {
      printf("Cannot open %s\n", IDfile);
      /* check if other IDfile exists */
      sprintf(IDfile, "%s/ID_level_%d_proc_%d.dat", outdir, level->l, rank);
      fp2 = fopen(IDfile, "rb");
      if(fp2==NULL)
      {
        printf("Cannot open %s\n", IDfile);
        printf("Will call interpolator.\n");
        use_interpolator = 1;
      }
      else
      {
        printf("%s exists, no need to call interpolator.\n", IDfile);
        use_interpolator = 0;
        fclose(fp2);
      }
    }
    else
    {
      printf("%s exists, no need to call interpolator.\n", IDfile);
      use_interpolator = 0;
      fclose(fp2);
    }
  }
  else 
    use_interpolator = 1;

  /* if this is a checkpoint restart do not interpolate */
  if(Getv("checkpoint","restart"))
  {
    sprintf(IDfile, "%s_previous/ID_level_%d_proc_%d.dat",
            outdir, level->l, rank);
    use_interpolator = 0;
  }

  /* check if we call the interpolator, or if we just read data */
  if(use_interpolator)
  {
    sprintf(IDfile, "%s/ID_level_%d_proc_%d.dat", outdir, level->l, rank);
    sprintf(IDfile_new, "%s%s", IDfile, "_new");

    /* info */
    printf("BBHdataReader:\n");
    printf("  sgridoutdir = %s\n", sgridoutdir);
    printf("  sgridcheckpoint_indir = %s\n", sgridcheckpoint_indir);

    /* remove any old sgridoutdir and make a new one */
    system2("rm -rf", sgridoutdir);
    //system2("rm -rf", sgridoutdir_previous);
    system2("mkdir", sgridoutdir);
    //system2("mkdir", sgridoutdir_previous);
    fflush(stdout);

    /* call the interpolator that will generate the ID */
    sprintf(call_interpolator, "%s %s/%s.par "
            "--modify-par:BBHdata_Interpolate_pointsfile=%s "
            "--modify-par:BBHdata_Interpolate_output=%s "
            "--modify-par:outdir=%s "
            "--modify-par:checkpoint_indir=%s",
            Gets("BBHdataReader_sgrid_exe"),
            Gets("BBHdataReader_sgrid_datadir"),
            Gets("BBHdataReader_sgrid_datadir"),
            gridfile, IDfile_new, sgridoutdir, sgridcheckpoint_indir);
    if(!Getv("BBHdataReader_keep_sgrid_output", "yes"))
      strcat(call_interpolator, " > /dev/null");
    printf("System call:\n%s\n", call_interpolator);
    fflush(stdout);
    ret = system(call_interpolator);
    if(ret) errorexit("interpolator returned non-zero exit code!");

    /* move IDfile_new to IDfile */
    rename(IDfile_new, IDfile);

    /* clean up */
    system2("rm -rf", sgridoutdir_previous);
    if(!Getv("BBHdataReader_keep_sgrid_output", "yes"))
      system2("rm -rf", sgridoutdir);
  }
  else
  {
    printf("Skipping interpolator, reading data from %s\n", IDfile);
  }

  /* info about filenames */
  printf("BBHdataReader:\n");
  printf("  gridfile = %s\n", gridfile);
  printf("  IDfile = %s\n", IDfile);
  fflush(stdout);

  /* read ADM variables from files generated by the interpolator */
  fp2 = fopen(IDfile, "rb");
  if(fp2==NULL) errorexits("could not open %s", IDfile);
  ndata = 16;
  j=BBH_position_fileptr_after_str(fp2, "$BEGIN_data:\n");
  if(j==EOF) errorexits("could not find $BEGIN_data: in %s", IDfile);
  forallpoints(level, i)
  {
    double detg;
    double dataline[NDATAMAX];

    /* read ndata doubles, i.e. one line of:
       alpha betax betay betaz gxx gxy gxz gyy gyz gzz 
       Kxx Kxy Kxz Kyy Kyz Kzz */
    n=fread(dataline, sizeof(double), ndata, fp2);
    j=0;
    alpha[i]=dataline[j];  j++;
    betax[i]=dataline[j];  j++;
    betay[i]=dataline[j];  j++;
    betaz[i]=dataline[j];  j++;
    gxx[i]=dataline[j];  j++;
    gxy[i]=dataline[j];  j++;
    gxz[i]=dataline[j];  j++;
    gyy[i]=dataline[j];  j++;
    gyz[i]=dataline[j];  j++;
    gzz[i]=dataline[j];  j++;
    Kxx[i]=dataline[j];  j++;
    Kxy[i]=dataline[j];  j++;
    Kxz[i]=dataline[j];  j++;
    Kyy[i]=dataline[j];  j++;
    Kyz[i]=dataline[j];  j++;
    Kzz[i]=dataline[j];  j++;

    if(j!=n) errorexits("error reading dataline from %s", IDfile);

    /* print g_ij, K_ij, beta^i, alpha */
    if(pr)
    {
      printf("(x,y,z)=(%g,%g,%g)\n", x[i],y[i],z[i]);

      printf("alpha=%g beta=%g %g %g\n",
	alpha[i], betax[i], betay[i], betaz[i]);  

      printf("g=%g %g %g %g %g %g \n", 
	gxx[i], gxy[i], gxz[i], gyy[i], gyz[i], gzz[i]);
    
      printf("K=%g %g %g %g %g %g\n",
	Kxx[i], Kxy[i], Kxz[i], Kyy[i], Kyz[i], Kzz[i]);

      printf("\n");
    }
    /* check if data makes sense */
    detg=(2.*gxy[i]*gxz[i]*gyz[i] + gxx[i]*gyy[i]*gzz[i] -
             gzz[i]*gxy[i]*gxy[i] - gyy[i]*gxz[i]*gxz[i] -
             gxx[i]*gyz[i]*gyz[i]);

    if(detg<=0)
    {
      printf("det(g_ij)=%g at ccc=i=%d:  x=%g y=%g z=%g\n", 
             detg, i, x[i], y[i], z[i]);
       //errorexit("found a point with det(g_ij)<=0.");
    }

    /* set psi and its derivs */
    psi[i] = 1.0;
    dpsiopsi1[i] = dpsiopsi2[i] = dpsiopsi3[i] = 0.0;
    ddpsiopsi11[i] = ddpsiopsi12[i] = ddpsiopsi13[i] = 0.0;
    ddpsiopsi22[i] = ddpsiopsi23[i] = ddpsiopsi33[i] = 0.0;
  }

  /* reset lapse */
  if(Getv("punctures_lapse", "rtoN_atPunc"))
    Set_alpha_rtoN_atPunc(level);

  prdivider(0);
  return 0;
}

/* get BBHdata pars */
int BBHdataPars(tL *level)
{
  FILE *fp1;
  char str[STRLEN];
  char datadir[STRLEN];
  double x_CM, xc1, xc2;
  int j;

  snprintf(datadir, STRLEN-1, "%s", Gets("BBHdataReader_sgrid_datadir"));
  strcat(datadir, "/BBHdata_properties.txt");

  /* open file */
  fp1 = fopen(datadir, "r");
  if(fp1==NULL) errorexits("could not open %s", datadir);
  j=BBH_position_fileptr_after_str(fp1, "BBH data properties (time = 0):\n");
  if(j==EOF) errorexits("could not find (time = 0) in %s", datadir);

  /* get pars from file */
  fgetparameter(fp1, "x_CM", str);
  Sets("BBHdataReader_sgrid_x_CM", str);
  x_CM=atof(str);
  fgetparameter(fp1, "m1", str);
  Sets("BBHdataReader_m1", str);
  fgetparameter(fp1, "m2", str);
  Sets("BBHdataReader_m2", str);
  /* shift xc1/2 such that CM is at 0 */
  fgetparameter(fp1, "xc1", str);
  Setd("BBHdataReader_xc1", atof(str)-x_CM);
  fgetparameter(fp1, "xc2", str);
  Setd("BBHdataReader_xc2", atof(str)-x_CM);

  /* close file */
  fclose(fp1);

  /* info */
  printf("BBHdataPars:\n");
  printf("m1 = %s\n", Gets("BBHdataReader_m1"));
  printf("m2 = %s\n", Gets("BBHdataReader_m2"));
  printf("xc1 - sgrid_x_CM = %s\n", Gets("BBHdataReader_xc1"));
  printf("xc2 - sgrid_x_CM = %s\n", Gets("BBHdataReader_xc2"));
  printf("sgrid_x_CM = %s\n", Gets("BBHdataReader_sgrid_x_CM"));

  /* set some pars relevant for setting up bam's grid */
  printf("Setting some pars relevant for setting up bam's grid:\n");
  Sets("bhmass1", Gets("BBHdataReader_m1"));
  Sets("bhmass2", Gets("BBHdataReader_m2"));
//Seti("bhmass2", 0);
  Sets("bhx1", Gets("BBHdataReader_xc1"));
  Sets("bhx2", Gets("BBHdataReader_xc2"));
  printf("bhmass1 = %s\n", Gets("bhmass1"));
  printf("bhmass2 = %s\n", Gets("bhmass2"));    
  printf("bhx1 = %s\n", Gets("bhx1"));
  printf("bhx2 = %s\n", Gets("bhx2"));    

  return 0;
}
