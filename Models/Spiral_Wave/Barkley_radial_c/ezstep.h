/*
 *  Define the macros U_THRESHOLD(v) and G(u,v) according to the model you 
 *  wish to simulate. In principle these can be almost anything you want.  
 *  Three examples are included.
 * 
 * $Revision: 2.3 $
 * $Date: 1997/11/13 23:17:22 $
 *---------------------------------------------------------------------------*/

#if 1  /* Standard Model */
#  define U_THRESHOLD(v)  ( one_o_a*(v) + b_o_a )
#  define G(u,v)          ( (u)-(v) )
#  define F(u, uth)        ( one_o_eps*(u)*(one-(u))*((u)-(uth)) ) 
#endif

#if 0  /* Simple Model for Turbulent Spirals */
#  define U_THRESHOLD(v)  ( one_o_a*(v) + b_o_a )
#  define G(u,v)          ( (u)*(u)*(u) - (v) )
#endif

#if 0  /*  Bar and Eiswirth Model (Phys. Rev. E V. 48, p. 1635, 1993).
	*  I do not include their case u>1 because it is not relevant.
	*/
#  define U_THRESHOLD(v)  ( one_o_a*(v) + b_o_a )
#  define G(u,v)          (-(v) + ( ((u)<1./3.) ? zero : \
                                    (1-6.75*(u)*(1-(u))*(1-(u))) ) )  
#endif

/*---------------------------------------------------------------------------*/
/*          In general you should not need to change anything below          */
/*---------------------------------------------------------------------------*/

/*  I rely on the compiler (gcc) to take care of most optimization.  
 *  The only place where I give the compiler help is with the nine-point 
 *  Laplacian formulas because apparently gcc does not do the obvious 
 *  thing.  I do not unroll the main loops in step() nor do I write explicit   
 *  pointer operations because gcc seems to handle all the indirections very 
 *  efficiently.
 */

/*---------------------------------------------------------------------------*/
/* 
 *  Defining the kinetics macros: U_KINETICS(u,uth) and V_KINETICS(u,v) 
 */

#define U_KINETICS(u,uth) ( (u)+dt_o_eps*(u)*(one-(u))*((u)-(uth)) )
#define V_KINETICS(u,v) ( (v)+dt*G(u,v) )

/*---------------------------------------------------------------------------*/
/* 
 *  Defining the diffusion macros: U_DIFFUSION,  V_DIFFUSION,  
 *                                 ADD_TO_U_LAPLACIAN, ADD_TO_V_LAPLACIAN, 
 *                                 and ZERO_USED_LAPLACIANS. 
 */

#define U_DIFFUSION   ( dt * u_lap[k1][i][j] ) 
#if V_DIFF_ON  
#  define V_DIFFUSION ( dt * Dv * v_lap[k1][i][j] )
#else   
#  define V_DIFFUSION zero
#endif

/* Polar Laplacian formulas */
#if POLAR
#  define ADD_TO_U_LAPLACIAN  \
     i2 = i*i; twoi=2*i;\
     u_lap[k2][i][j] -= (two_o_dx1pow2 + two_o_dx1pow2dx2pow2/i2)*u[i][j]; \
     u_lap[k2][i+1][j] += one_o_dx1pow2 * (twoi + 1) * u[i][j] / (twoi+2); \
     if (i>1) u_lap[k2][i-1][j] += one_o_dx1pow2 * (twoi - 1) * u[i][j] / (twoi-2); \
     u_lap[k2][i][j+1] += one_o_dx1pow2dx2pow2 * u[i][j] / i2; \
     u_lap[k2][i][j-1] += one_o_dx1pow2dx2pow2 * u[i][j] / i2;  
#  if V_DIFF_ON 
#    define ADD_TO_V_LAPLACIAN       \
       v_lap[k2][i][j]   -= (two_o_dx1pow2 + two_o_dx1pow2dx2pow2/i2)*v[i][j]; \
       v_lap[k2][i+1][j] += one_o_dx1pow2 * (twoi + 1) * v[i][j] / (twoi+2); \
       if (i>1) v_lap[k2][i-1][j] += one_o_dx1pow2 * (twoi - 1) * v[i][j] / (twoi-2); \
       v_lap[k2][i][j+1] += one_o_dx1pow2dx2pow2 * v[i][j] / i2; \
       v_lap[k2][i][j-1] += one_o_dx1pow2dx2pow2 * v[i][j] / i2; 
#  else
#    define ADD_TO_V_LAPLACIAN 
#  endif
/* for adding boundary terms to laplacian */
#define U_BOUNDARY_BOTTOM ( u[i][0] / (dx1pow2dx2pow2*i*i) )
#define U_BOUNDARY_TOP    ( u[i][M2+1] / (dx1pow2dx2pow2*i*i) )
#define U_BOUNDARY_RIGHT  ( bcfac * u[M1+1][j] )
#define U_BOUNDARY_LEFT   ( one_o_twodx1pow2 * u[0][j] )

#   if V_DIFF_ON
#define V_BOUNDARY_BOTTOM ( v[i][0] / (dx1pow2dx2pow2*i*i) )
#define V_BOUNDARY_TOP    ( v[i][M2+1] / (dx1pow2dx2pow2*i*i) )
#define V_BOUNDARY_RIGHT  ( bcfac * v[M1+1][j] )
#define V_BOUNDARY_LEFT   ( one_o_twodx1pow2 * v[0][j] )
#   endif

#else /* cartesian coords */

#if NINEPOINT
/* Nine-point Laplacian formulas */
/* assume dx1 == dx2 */
#  define ADD_TO_U_LAPLACIAN           \
    { \
      u_lap[k2][i][j] -= twenty*u[i][j] * one_o_sixdxpow2;  \
      u_lap[k2][i+1][j] +=  four*u[i][j] * one_o_sixdxpow2;  \
      u_lap[k2][i-1][j] +=  four*u[i][j] * one_o_sixdxpow2;  \
      u_lap[k2][i][j+1] +=  four*u[i][j] * one_o_sixdxpow2;  \
      u_lap[k2][i][j-1] +=  four*u[i][j] * one_o_sixdxpow2;  \
      u_lap[k2][i+1][j+1] += u[i][j] * one_o_sixdxpow2;  \
      u_lap[k2][i-1][j+1] += u[i][j] * one_o_sixdxpow2;  \
      u_lap[k2][i+1][j-1] += u[i][j] * one_o_sixdxpow2;  \
      u_lap[k2][i-1][j-1] += u[i][j] * one_o_sixdxpow2;  \
     }
#  if V_DIFF_ON 
#    define ADD_TO_V_LAPLACIAN            \
       { \
        v_lap[k2][i][j] -= twenty*v[i][j] * one_o_sixdxpow2;  \
        v_lap[k2][i+1][j] += four*v[i][j] * one_o_sixdxpow2;   \
        v_lap[k2][i-1][j] += four*v[i][j] * one_o_sixdxpow2;   \
        v_lap[k2][i][j+1] += four*v[i][j] * one_o_sixdxpow2;   \
        v_lap[k2][i][j-1] += four*v[i][j] * one_o_sixdxpow2;   \
        v_lap[k2][i+1][j+1] += v[i][j] * one_o_sixdxpow2;   \
        v_lap[k2][i-1][j+1] += v[i][j] * one_o_sixdxpow2;   \
        v_lap[k2][i+1][j-1] += v[i][j] * one_o_sixdxpow2;   \
        v_lap[k2][i-1][j-1] += v[i][j] * one_o_sixdxpow2;   \
       }
#  else
#    define ADD_TO_V_LAPLACIAN 
#  endif
/* for adding boundary terms to laplacian */
#define U_BOUNDARY_BOTTOM (one_o_sixdxpow2 * (u[i-1][0] + four * u[i][0] + u[i+1][0]))
#define U_BOUNDARY_TOP    (one_o_sixdxpow2 * (u[i-1][M2+1] + four * u[i][M2+1] + u[i+1][M2+1]))
#define U_BOUNDARY_RIGHT  (one_o_sixdxpow2 * (u[M1+1][j-1] + four * u[M1+1][j] + u[M1+1][j+1]))
#define U_BOUNDARY_LEFT   (one_o_sixdxpow2 * (u[0][j-1] + four * u[0][j] + u[0][j+1]))
#   if V_DIFF_ON
#define V_BOUNDARY_BOTTOM (one_o_sixdxpow2 * (v[i-1][0] + four * v[i][0] + v[i+1][0]))
#define V_BOUNDARY_TOP    (one_o_sixdxpow2 * (v[i-1][M2+1] + four * v[i][M2+1] + v[i+1][M2+1]))
#define V_BOUNDARY_RIGHT  (one_o_sixdxpow2 * (v[M1+1][j-1] + four * v[M1+1][j] + v[M1+1][j+1]))
#define V_BOUNDARY_LEFT   (one_o_sixdxpow2 * (v[0][j-1] + four * v[0][j] + v[0][j+1]))
#   endif
#else 
/* Five-point Laplacian formulas */
#  define ADD_TO_U_LAPLACIAN       \
     u_lap[k2][i][j] -= (two/dx1pow2+two/dx2pow2)*u[i][j]; \
     u_lap[k2][i+1][j] += u[i][j]/dx1pow2; \
     u_lap[k2][i-1][j] += u[i][j]/dx1pow2; \
     u_lap[k2][i][j+1] += u[i][j]/dx2pow2; \
     u_lap[k2][i][j-1] += u[i][j]/dx2pow2;  
#  if V_DIFF_ON 
#    define ADD_TO_V_LAPLACIAN       \
       v_lap[k2][i][j] -= (two/dx1pow2+two/dx2pow2)*v[i][j]; \
       v_lap[k2][i+1][j] += v[i][j]/dx1pow2; \
       v_lap[k2][i-1][j] += v[i][j]/dx1pow2; \
       v_lap[k2][i][j+1] += v[i][j]/dx2pow2; \
       v_lap[k2][i][j-1] += v[i][j]/dx2pow2; 
#  else
#    define ADD_TO_V_LAPLACIAN 
#  endif
/* for adding boundary terms to laplacian */
#define U_BOUNDARY_BOTTOM (one_o_dx2pow2 * u[i][0])
#define U_BOUNDARY_TOP    (one_o_dx2pow2 * u[i][M2+1])
#define U_BOUNDARY_RIGHT  (one_o_dx1pow2 * u[M1+1][j])
#define U_BOUNDARY_LEFT   (one_o_dx1pow2 * u[0][j])
#   if V_DIFF_ON
#define V_BOUNDARY_BOTTOM (one_o_dx2pow2 * v[i][0])
#define V_BOUNDARY_TOP    (one_o_dx2pow2 * v[i][M2+1])
#define V_BOUNDARY_RIGHT  (one_o_dx1pow2 * v[M1+1][j])
#define V_BOUNDARY_LEFT   (one_o_dx1pow2 * v[0][j])
#   endif
# endif /* nine/five point laplacian */
#endif /* cart coords */ 

#if V_DIFF_ON 
#  define ZERO_USED_LAPLACIANS \
    u_lap[k1][i][j] = zero;       \
    v_lap[k1][i][j] = zero; 
#else  
#  define ZERO_USED_LAPLACIANS \
    u_lap[k1][i][j] = zero; 
#endif

/*---------------------------------------------------------------------------*/
