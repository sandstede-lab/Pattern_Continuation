/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*                                                                           */
/*                                EZ-SPIRAL                                  */
/*                    A Code for Simulating Spiral Waves                     */
/*                                                                           */
/*           Copyright (C) 1992, 1994, 1997 Dwight Barkley                   */
/*                                          barkley@maths.warwick.ac.uk      */
/*                                                                           */
/*                                Version 2                                  */
/*                                                                           */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* 
 * $Revision: 2.4 $
 * $Date: 1997/11/13 23:22:00 $
 */
/*---------------------------------------------------------------------------*/

#include "ezspiral.h"
#include "ezstep.h"
#include "math.h"
#include "mex.h"

/* Global variables used throughout the EZ-Spiral Code */
precision  **u, **v,                    /* u and v fields */
           **u_lap[2], **v_lap[2];      /* u and v Laplacian fields */
  
precision  *tip_x, *tip_y,              /* spiral tip arrays */
           a_param, b_param,            /* model parameters */
           one_o_eps, length1, length2, Dv, /* model parameters */
           ts,                          /* numerical parameters */
           dx1, dx2, dt,t,                  /* numerical parameters */
           L1div2, L2div2, dx1dx2,   
           one_o_a, b_o_a,              /* useful parameter combinations */
           dt_o_eps, dt_o_2eps, dx1pow2, dx2pow2,         /* useful parameter combinations */
           dx1_o_dx2,
           one_o_dx2pow2, one_o_dx1pow2; 
precision  delta, dt, t, one_o_a, b_o_a, dt_o_eps, dt_o_2eps,
           dx1, dx2, L1div2, L2div2,  Dv, dx1pow2, dx2pow2, one_o_dx2pow2, one_o_dx1pow2, dx1_o_dx2,
           r, phi, sum;                                

int        M1, M2,                      /* # of grid points per direction */
           u_plot, v_plot,              /* flags indicating fields to plot */
           verbose,                     /* verbosity level */
           k1, k2, ntips, itip, iter,         /* useful integers */
           i2, twoi;

int        screenshot;                  /* make screenshot */
#if POLAR
int        N4;  /* M2/4 */
precision dx1pow2dx2pow2, bcfac, one_o_twodx1pow2, two_o_dx1pow2, one_o_dx1pow2dx2pow2, two_o_dx1pow2dx2pow2;
precision  twodx2_o_Pidx1pow2;
static void       polarbc();
#endif

#if NINEPOINT 
precision one_o_sixdxpow2;
#endif

/* Functions local to ezstep */
#if NBC
static void       nbc();
#else
static void       pbc();
#endif

/* constants */
const precision zero   = 0.;
const precision half   = 0.5;
const precision one    = 1.;
const precision two    = 2.;
const precision three  = 3.;
const precision four   = 4.;
const precision twenty = 20.;


/* Global variables and functions for this file only */
static int      nits, iplot, iwrt, hist_x, hist_y;
static FILE     *history, *path;
/* static void     initialize(); */
static void     initialize(double *U_in, double *V_in);
static void     finish_up ();
static void     write_dat (int iter);
static void     write_fc  (char *filename);
/* static void     read_ic   (char *filename); */
static void     read_ic   (double *U_in, double *V_in);
static void     next_line (FILE *filep);

static char     cmd[60];
static int      frame=0;

/*---------------------------------------------------------------------------*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

     /* 
      *  Main routine. Takes nits time steps of the model. 
      */
    double *U_in = mxGetPr( prhs[0] );
    double *V_in = mxGetPr( prhs[1] );
    printf("%f\n",U_in[0]);
    
    
    initialize(U_in, V_in);
    plhs[0] = mxCreateDoubleMatrix(  M2*M1+1, 1, mxREAL); 
    plhs[1] = mxCreateDoubleMatrix( M2*M1+1, 1, mxREAL); 
    double *U = mxGetPr(plhs[0]);
    double *V = mxGetPr(plhs[1]);
    for(iter=0;iter<nits;iter++) {
      /*if(keyboard_chk())         break;*/
      if(iwrt && (iter%iwrt)==0) write_dat(iter);

    /*  if(itip && (iter%itip)==0) 
	tip_finder(); */

    /*if((iter%iplot)==0)        
      {
	plot();*/
	
	/* make screenshot */
	/* if(screenshot) save_image(iter); */
     /* } */
    
    
    step();
  }
   
   finish_up(); 
   for(mwSize i=1;i<=M1;++i) {
      for(mwSize j=1;j<=M2;++j) {
        U[(i-1)*M2+j] = u[i][j]; 
        V[(i-1)*M2+j] = v[i][j]; 
      }
  }
   U[M1*M2+1] = u[0][1];
   V[M1*M2+1] = v[0][1];
  /* return 0; */
}
/*---------------------------------------------------------------------------*/

void initialize(double *U_in, double *V_in)

     /* 
      *  Opens files, reads data, initializes graphics and time stepping. 
      */
{
  double p_in;
  FILE *fp;
  int i;

  /* Read parameters from task.dat */
  if((fp=fopen("task.dat","r"))==NULL) {
    if(verbose) fprintf(stderr,"Cannot open task file: task.dat \n");
    exit(1);
  }
  else {
    fscanf(fp,"%lg",&p_in); next_line(fp); a_param=p_in;
    fscanf(fp,"%lg",&p_in); next_line(fp); b_param=p_in;
    fscanf(fp,"%lg",&p_in); next_line(fp); one_o_eps=p_in;
    fscanf(fp,"%lg",&p_in); next_line(fp); length1=p_in;
    fscanf(fp,"%lg",&p_in); next_line(fp); length2=p_in;
    fscanf(fp,"%lg",&p_in); next_line(fp); Dv=p_in;
    next_line(fp);
    fscanf(fp,"%d", &M1); next_line(fp); 
    fscanf(fp,"%d", &M2); next_line(fp); 
    fscanf(fp,"%lg",&p_in); next_line(fp); ts=p_in;
    next_line(fp);
    fscanf(fp,"%d",&nits);  next_line(fp);
    fscanf(fp,"%d",&iplot); next_line(fp);
    fscanf(fp,"%d",&itip);  next_line(fp);
    fscanf(fp,"%d",&iwrt);  next_line(fp);
    fscanf(fp,"%d,%d",&hist_x,&hist_y); next_line(fp);
    next_line(fp);
    fscanf(fp,"%d",&u_plot);next_line(fp);
    fscanf(fp,"%d",&v_plot);next_line(fp);
    fscanf(fp,"%d",&verbose);next_line(fp);
    next_line(fp);
    next_line(fp);
    next_line(fp);
    next_line(fp);
    fscanf(fp,"%d",&screenshot);
    fclose(fp);
  }

#define STABILITY_LIMIT ( 0.25*min(dx1*dx1, dx2*dx2) )

 
#if POLAR /* periodic top, bottom, Neumann right */
  length2 = 2*M_PI;
  dx1   = length1/M1;     /* 0..M1 -> 0..R */
  dx2   = length2/M2;     /* 1..M2+1 -> 0..2 Pi */
  twodx2_o_Pidx1pow2 = dx2 * two_o_dx1pow2 / M_PI;

#else /* cart. coords */
# if NBC /* if Neumann  boundary conditions are used */
  dx1   = length1/(M1-1);
  dx2   = length2/(M2-1);
# else /* Periodic boundary conditions are used */
  dx1   = length1/M1;
  dx2   = length2/M2;
# endif
#endif
 
#if  DEBUG
  dt=0.1; /* debug */
#else
    dt = ts*STABILITY_LIMIT;
#endif
  one_o_a   = 1.0/a_param; 
  b_o_a     = b_param/a_param; 
  dt_o_eps  = dt*one_o_eps; 
  dt_o_2eps = dt_o_eps/2.0;

  L1div2 = length1/2;
  L2div2 = length2/2;

  dx1pow2 = dx1*dx1;
  dx2pow2 = dx2*dx2;
  dx1dx2 = dx1*dx2;
  dx1_o_dx2 = dx1/dx2;
  one_o_dx2pow2 = 1/dx2pow2;
  one_o_dx1pow2 = 1/dx1pow2;

#if NINEPOINT
  one_o_sixdxpow2 = one_o_dx1pow2/6;
#endif

#if POLAR
  N4 = M2/4;
  dx1pow2dx2pow2 = dx1pow2*dx2pow2;
  bcfac = (2*M1+1)/(2*M1*dx1pow2);
  one_o_twodx1pow2 = 1/(2*dx1pow2);
  two_o_dx1pow2 = 2*one_o_dx1pow2;
  one_o_dx1pow2dx2pow2 = 1 / dx1pow2dx2pow2;
  two_o_dx1pow2dx2pow2 = 2 * one_o_dx1pow2dx2pow2;
#endif


  if(verbose) {
    printf("\n\nModel Parameters: \n");
    printf("a     = %g\n", a_param);
    printf("b     = %g\n", b_param);
    printf("eps   = 1/%g = %g\n", one_o_eps, 1.0/one_o_eps);
    printf("L1     = %g\n", length1);
    printf("L2     = %g\n", length2);
    printf("Dv    = %g\n", Dv);
    printf("\nNumerical Parameters: \n");
    printf("M1 = %-10d M2 = %-10d ts     = %-10g\n", 
	   M1, M2, ts);
    printf("dt   = %-10g dt/eps = %-10g dx1 = %-10g dx2 = %-10g\n", 
	   dt, dt_o_eps, dx1, dx2);
    printf("dt/(dx1*dx1) = %-10g dt/(dx2*dx2) = %-10g ", dt/(dx1*dx1),dt/(dx2*dx2)  );
    printf("\n\nNumber of time steps = %d\n", nits);
    if(u_plot||v_plot) printf("    time steps per plot = %d\n", iplot);
    if(itip)   printf("    time steps per locating tip = %d\n", itip);
    if(iwrt)   printf("    time steps per write to data files= %d\n", iwrt);
    if(hist_x) printf("    history point = (%d, %d)\n", hist_x, hist_y ); 
  }

 
  /*f( V_DIFF_ON && Dv==0.) {
    fprintf(stderr,"***** V_DIFF_ON is 1 and Dv == 0. ******\n");
    exit(1);
  }
  */

  if(!V_DIFF_ON && Dv!=0.) {
    fprintf(stderr,"***** V_DIFF_ON is 0 and Dv != 0. ******\n");
    exit(1);
  }
  if(hist_x<0 || hist_x>M1 || hist_y<0 || hist_y>M2) {
    fprintf(stderr,"***** history point out of range ******\n");
    exit(1);
  }
  if(ts > 1.0 ) {
    fprintf(stderr,"***** ts > 1 (the diffusion stability limit) ******\n");
    exit(1);
  }

#if POLAR
  if( div(M2,4).rem !=0) {
    fprintf(stderr,"***** M2 has to be a multiple of four for polar coords ******\n");
    exit(1);
  }
#endif

#if NINEPOINT
  if( dx1 != dx2 ) {
    fprintf(stderr,"***** Ninepoint Laplacian only defined if dx1 = dx2 ******\n");
    exit(1);
  }  
#endif

  /*  Allocate memory for u[M1+2][M2+2], v[M1+2][M2+2], 
   *  u_lap[2][M1+2][M2+2], and if V_DIFF_ON, for v_lap[2][M1+2][M2+2]. 
   *  u and v could be of size [M1][M2]; however, because u_lap and v_lap 
   *  MUST be of size [2][M1+2][M2+2], I think it best that u and v also
   *  be of size [M1+2][M2+2].  The memory for each field is allocated 
   *  with the (M1+2)*(M2+2) locations contiguous.
   */

  u = (precision **) malloc((unsigned)(M1+2)*sizeof(precision*));
  v = (precision **) malloc((unsigned)(M1+2)*sizeof(precision*));
  u[0] = (precision *) calloc((M1+2)*(M2+2), sizeof(precision));
  v[0] = (precision *) calloc((M1+2)*(M2+2), sizeof(precision));
  for(i=1;i<=M1+1;i++){
    u[i] = u[0]+i*(M2+2);
    v[i] = v[0]+i*(M2+2);
  }

  u_lap[0] = (precision **) malloc((unsigned)(M1+2)*sizeof(precision*));
  u_lap[1] = (precision **) malloc((unsigned)(M1+2)*sizeof(precision*));
  u_lap[0][0] = (precision *) calloc((M1+2)*(M2+2), sizeof(precision));
  u_lap[1][0] = (precision *) calloc((M1+2)*(M2+2), sizeof(precision));
  for(i=1;i<=M1+1;i++){
    u_lap[0][i] = u_lap[0][0]+i*(M2+2);
    u_lap[1][i] = u_lap[1][0]+i*(M2+2);
  }

#if V_DIFF_ON
  v_lap[0] = (precision **) malloc((unsigned)(M1+2)*sizeof(precision*));
  v_lap[1] = (precision **) malloc((unsigned)(M1+2)*sizeof(precision*));
  v_lap[0][0] = (precision *) calloc((M1+2)*(M2+2), sizeof(precision));
  v_lap[1][0] = (precision *) calloc((M1+2)*(M2+2), sizeof(precision));
  for(i=1;i<=M1+1;i++){
    v_lap[0][i] = v_lap[0][0]+i*(M2+2);
    v_lap[1][i] = v_lap[1][0]+i*(M2+2);
  }
#endif

  ntips = 0;
  if(itip) {  /* memory for tip path. */
    tip_x=(precision *)malloc((unsigned)(nits+1)*sizeof(precision));
    tip_y=(precision *)malloc((unsigned)(nits+1)*sizeof(precision));
  }
  if(iwrt) {  /* open history and tip-path files. */
    if(hist_x!=0) history = fopen("history.dat", "w");
    if(itip  !=0) path    = fopen("path.dat",    "w");
  }

  /* Read initial conditions and initialize time stepping and graphics */
  /*read_ic("ic.dat");  */
  read_ic(U_in, V_in);
  step_ini();
  /* plot_ini(); */
}
/*---------------------------------------------------------------------------*/

void finish_up()

     /* 
      *  Writes final conditions, possibly writes last history and/or 
      *  tip-path data, and closes data files.
      */
{
  write_fc("fc.dat");
  if(iwrt) {
    if((nits%iwrt)==0) write_dat(nits); 
    if(hist_x!=0) fclose(history);
    if(itip  !=0) fclose(path);
  }
}
/*---------------------------------------------------------------------------*/

void write_dat(int iter)

     /* 
      *  Writes history and tip-path data. 
      */
{
  FILE *fp; /* neu */
  char filename[80];
  int i,j;
  
  if(verbose>1) printf("iter = %d", iter);

  if(hist_x!=0) {
    fprintf(history, "%.5f %.5f %.5f \n", dt*(precision)iter, 
	    u[hist_x][hist_y], v[hist_x][hist_y]); 
    if(verbose>1) printf("; history (u,v) = %.5f, %.5f", 
			 u[hist_x][hist_y], v[hist_x][hist_y]); 
  }

  if(ntips>0) {
    fprintf(path, "%.5f %.5f %.5f \n", dt*(precision)iter, 
	    dx1*(tip_x[ntips-1]-1), dx2*(tip_y[ntips-1]-1) );
    if(verbose>1) printf("; tip (x,y) = %.5f %.5f", 
	     dx1*(tip_x[ntips-1]-1), dx2*(tip_y[ntips-1]-1) );
  }
  if(verbose>1) printf("\n"); 

  /* neu: schreibe frames in Datei */

#if 0
  sprintf(filename, "frame%d.dat", iter/iwrt);
  fp = fopen(filename, "w");
  
  fprintf(fp,"Model Parameters: a, b, 1/eps, L1, L2, Dv = %g, %g, %g, %g, %g, %g\n",
	  a_param, b_param, one_o_eps, length1, length2, Dv);
  fprintf(fp,"Numerical Parameters: M1, M2, ts = %d, %d, %g \n", 
	  M1, M2, ts);
  fprintf(fp,"%e \n", t); 

#if POLAR
  fprintf(fp,"r\tphi\tu\tv\n");
#else
  fprintf(fp,"x\ty\tu\tv\n");
#endif

#if POLAR
  for(i=0;i<=M1;i++) 
    for(j=0;j<=M2;j++) 
      fprintf (fp, "%e %e %e %e\n", dx1*i, dx2*(j-1), u[i][j], v[i][j]);
#else
  for(i=1;i<=M1;i++) 
    for(j=1;j<=M2;j++) 
      fprintf (fp, "%e %e %e %e\n", (i-1)*dx1 - L1div2, (j-1)*dx2 - L2div2, u[i][j], v[i][j]);
#endif
  
  fclose(fp);
  /* zip output */
  sprintf(cmd, "gzip -f frame%d.dat",iter/iwrt);
  system(cmd);
  
#endif
}
/*---------------------------------------------------------------------------*/

void write_fc (char *filename)

     /* 
      *  Write final conditions to a file. 
      */
{
  FILE *fp;
  time_t tt1;
  int i,j;

  time(&tt1); 
  fp = fopen(filename, "w");
  fprintf(fp,"Model Parameters: a, b, 1/eps, L1, L2, Dv = %g, %g, %g, %g, %g, %g\n",
	  a_param, b_param, one_o_eps, length1, length2, Dv);
  fprintf(fp,"Numerical Parameters: M1, M2, ts = %d, %d, %g \n", 
	  M1, M2, ts);
  fprintf(fp,"File written: %s",ctime(&tt1));
  fprintf(fp,"Comments:\n");
  fprintf(fp,"Values of u and v follow\n");
  for(i=1;i<=M1;i++) {
    for(j=1;j<=M2;j++) {
      fprintf (fp, "%e %e\n", u[i][j], v[i][j]);
      
    }
  }
#if POLAR
  fprintf (fp, "%e %e\n", u[0][1], v[0][1]); /* central point */
#endif
  fclose(fp);
  
}
/*---------------------------------------------------------------------------*/


void read_ic(double *U_in, double *V_in)    
     
{
  /*precision **u_tmp, **v_tmp, ratio1, ratio2; */
  double u_in, v_in;
  int M1_ic, M2_ic, i_tmp, j_tmp, dummy, i, j;
  
  
  
    /*u =((precision **) malloc((unsigned)(M1)*sizeof(precision*)))-1;
    v =((precision **) malloc((unsigned)(M1)*sizeof(precision*)))-1;
    for(i=1;i<=M1;i++){
      u[i]=((precision *) malloc((unsigned)(M2)*sizeof(precision)))-1;
      v[i]=((precision *) malloc((unsigned)(M2)*sizeof(precision)))-1;
    } */

    
    for(i=1;i<=M1;i++) {
      for(j=1;j<=M2;j++) {
	
    u[i][j] = U_in[(i-1)*M2+j]; v[i][j] = V_in[(i-1)*M2+j]; 
      }
    }
    
#if POLAR
    
     u[0][1] = U_in[M1*M2+1]; v[0][1] = V_in[M1*M2+1]; 
#endif
    
    
    
    /*t = 0;*/
}

/*void read_ic (double *U_in, double *V_in)
{
  precision **u_tmp, **v_tmp, ratio1, ratio2;
  double u_in, v_in;
  int M1_ic, M2_ic, i_tmp, j_tmp, dummy, i, j;
for(i=0;i<=M1+1;i++) {
      for(j=1;j<=M2;j++) {
#if POLAR
	if(j<2*N4+1) v[i][j] =  a_param/2; else v[i][j] = 0.;  
	if((j<N4+1) || (j>3*N4)) u[i][j] = 1.0; else u[i][j] = 0.;
	u[0][1] = 1; v[0][1] = a_param/2; 
#else
	if(i<(M1/2)) u[i][j] = 0.; else u[i][j] = 1.0; 
	if(j<M2/2)   v[i][j] = 0.; else v[i][j] = a_param/2.; 
#endif
	
      }
    }
   
    t=0;
    }
*/
/*---------------------------------------------------------------------------*/

void next_line(FILE *filep)
     
     /* 
      *  Skips to next input line. 
      */
{
  int dummy;
  while( (dummy=getc(filep)) != '\n');
}
/*---------------------------------------------------------------------------*/
void step()

     /* 
      *  Routine for taking one time step. 
      *  See ezstep.h for definitions of various macros used here.
      */
{
  precision u_thresh;
  register int i,j;
  
  /* Interchange k1 and k2 */
  k1 = 1-k1;
  k2 = 1-k2;
  
#if POLAR
  /* compute central point u[0][1] */
  u_thresh = U_THRESHOLD(v[0][1]);
  v[0][1]  = V_KINETICS(u[0][1],v[0][1]);
  u[0][1]  = U_KINETICS(u[0][1],u_thresh) + dt * u_lap[k1][0][1];
#  if V_DIFF_ON 
  v[0][1] += dt * v_lap[k1][0][1];
#  endif
#endif

  /* Main Loop (ALMOST ALL WORK DONE HERE !!!) */
  for(i=1; i<=M1; i++)
    for(j=1; j<=M2;j++) 
      {
	u_thresh = U_THRESHOLD(v[i][j]);
	v[i][j]  = V_KINETICS(u[i][j],v[i][j])  + V_DIFFUSION;
	u[i][j]  = U_KINETICS(u[i][j],u_thresh) + U_DIFFUSION;

	ADD_TO_U_LAPLACIAN 
	ADD_TO_V_LAPLACIAN
	ZERO_USED_LAPLACIANS
      }
  
#if POLAR  
  /* set boundary conditions for polar coords */
  /* for central point at the left, periodic up down */
  polarbc();
#endif 

#if NBC
    nbc(); 
#else    
    pbc();
#endif

    /* add boundary points to laplacian */
    /* top bottom */
    for(i=1;i<=M1;i++) 
      { 
	u_lap[k2][i][1]    +=  U_BOUNDARY_BOTTOM;
	u_lap[k2][i][M2]   +=  U_BOUNDARY_TOP; 
#if V_DIFF_ON 
	v_lap[k2][i][1]    +=  V_BOUNDARY_BOTTOM;
	v_lap[k2][i][M2]   +=  V_BOUNDARY_TOP; 
#endif
      } 
    /* left right */
    for(j=1;j<=M2; j++) 
      { 
	u_lap[k2][1][j]  +=  U_BOUNDARY_LEFT;  
	u_lap[k2][M1][j] +=  U_BOUNDARY_RIGHT; 
#if V_DIFF_ON 
	v_lap[k2][1][j]  +=  V_BOUNDARY_LEFT;  
	v_lap[k2][M1][j] +=  V_BOUNDARY_RIGHT; 
#endif
      }

}    
/*---------------------------------------------------------------------------*/

void step_ini()

     /* 
      *  Routine for initializing Laplacians. 
      */
{
  int i, j;
  
  /* Set initial k1 and k2 */
  k1 = 0;
  k2 = 1;
  
  /* Compute Laplacians of initial fields */
  for(i=1;i<=M1;i++) {
    for(j=1;j<=M2;j++) {
      ADD_TO_U_LAPLACIAN; 
      ADD_TO_V_LAPLACIAN;
      ZERO_USED_LAPLACIANS;
    }
  }

#if POLAR  
     /* set boundary conditions for polar coords */
     /* for central point at the left, periodic up down */
     polarbc();
#endif


#if NBC
     nbc(); 
#else
     pbc();
#endif

     /* add boundary points to laplacian */
     /* top bottom */
     for(i=1;i<=M1;i++) 
       { 
	 u_lap[k2][i][1]    +=  U_BOUNDARY_BOTTOM;
	 u_lap[k2][i][M2]   +=  U_BOUNDARY_TOP; 
#if V_DIFF_ON 
	 v_lap[k2][i][1]    +=  V_BOUNDARY_BOTTOM;
	 v_lap[k2][i][M2]   +=  V_BOUNDARY_TOP; 
#endif
       } 
     /* left right */
     for(j=1;j<=M2; j++) 
       { 
	 u_lap[k2][1][j]  +=  U_BOUNDARY_LEFT;  
	 u_lap[k2][M1][j] +=  U_BOUNDARY_RIGHT; 
#if V_DIFF_ON 
	 v_lap[k2][1][j]  +=  V_BOUNDARY_LEFT;  
	 v_lap[k2][M1][j] +=  V_BOUNDARY_RIGHT; 
#endif
       }
     
}

void nbc()

     /* 
      *  Imposes Neumann boundary conditions. 
      */
{
  int j;
   
#if POLAR

  for(j=1;j<=M2;j++) /* neumann right */
    {
      u[M1+1][j] = u[M1-1][j];
    }
  
#else /* cart coords */
  int i;
  
  /* update solution */
 for(i=1;i<=M1;i++) /* neumann bottom, top */
  {
      u[i][0] = u[i][2];
      v[i][0] = v[i][2];
      
      u[i][M2+1] = u[i][M2-1];
      v[i][M2+1] = v[i][M2-1];
    }

  for(j=1;j<=M2;j++) /* neumann left, right */
    {
      u[0][j] = u[2][j];
      v[0][j] = v[2][j];
    
      u[M1+1][j] = u[M1-1][j];
      v[M1+1][j] = v[M1-1][j];
    }
			       
#if NINEPOINT
  /* set edge points */
  u[0][0] = u[1][1]/2; u[M1+1][0] = u[M1][1]/2;
  u[0][M2+1] = u[1][M2]/2; u[M1+1][M2+1] = u[M1][M2]/2;
#endif

#   endif /* cart coord */

}
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

void pbc()

     /* 
      *  Imposes Periodic boundary conditions. 
      */
{
  int i, j;
  
  for(i=1;i<=M1;i++) /* periodic up down */
    {
      u[i][0] = u[i][M2]; v[i][0] = v[i][M2];
      u[i][M2+1] = u[i][1]; v[i][M2+1] = v[i][1]; 
    }
  
  for(j=1;j<=M2;j++) /* periodic left, right */
    {
      u[0][j] = u[M1][j]; v[0][j] = v[M1][j]; 
      u[M1+1][j] = u[1][j];  v[M1+1][j] = v[1][j];
    }
}

#if POLAR
/* set boundary conditions  for polar coords */
void polarbc()
{
  int i,j;
  double avg_u = 0;
#if V_DIFF_ON
  double avg_v = 0;
#endif

  for(i=1;i<=M1;i++) /* periodic up down */
    {
      u[i][0] = u[i][M2]; v[i][0] = v[i][M2];
      u[i][M2+1] = u[i][1]; v[i][M2+1] = v[i][1];
    }

  for(j=1;j<=M2;j++) /* center left */
    {
      u[0][j] = u[0][1]; v[0][j] = v[0][1];
    }
  
  /* set laplacian at central point */
  /* u_lap[k2][0][1] = one_o_dx1pow2 * (u[1][1] + u[1][2*N4+1] - four*u[0][1] + u[1][N4+1] + u[1][3*N4+1]); */

  for(j=1; j<=M2; j++)
    {
      avg_u += u[1][j];
#  if V_DIFF_ON 
      avg_v += v[1][j];
#  endif
    }
  
  u_lap[k2][0][1] = one_o_dx1pow2 * 4 * (avg_u/M2 - u[0][1]);
#  if V_DIFF_ON 
  v_lap[k2][0][1] = one_o_dx1pow2 * 4 * (avg_v/M2 - v[0][1]);
#  endif

}

#endif
/*---------------------------------------------------------------------------*/

