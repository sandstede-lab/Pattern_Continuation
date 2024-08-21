/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*                                                                           */
/*                                EZSTEP.C                                   */
/*                  Time-stepping Subroutines for EZ-SPIRAL                  */
/*                                                                           */
/*                                Version 2                                  */
/*                                                                           */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* 
 * $Revision: 2.3 $
 * $Date: 1997/11/13 23:18:21 $
 */
/*---------------------------------------------------------------------------*/

#include "ezspiral.h"
#include "ezstep.h"
#include "/Applications/MATLAB_R2021a.app/extern/include/mex.h"
/* External variables */
extern precision  **u, **v, **u_lap[2], **v_lap[2],
                  delta, dt, t, one_o_a, b_o_a, dt_o_eps, dt_o_2eps,
                  dx1, dx2, L1div2, L2div2,  Dv, dx1pow2, dx2pow2, one_o_dx2pow2, one_o_dx1pow2, dx1_o_dx2; 
                                     
extern int        M1, M2, k1, k2;

#if POLAR 
extern int        N4;
int                i2, twoi;
precision         r, phi, sum;
extern precision  dx1pow2dx2pow2, bcfac, one_o_twodx1pow2, two_o_dx1pow2, one_o_dx1pow2dx2pow2, two_o_dx1pow2dx2pow2;
extern precision  twodx2_o_Pidx1pow2;
static void       polarbc();
#endif

#if NINEPOINT
extern precision one_o_sixdxpow2;
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

