/*
 *  The parameters below affect the compilation of EZ-Spiral.
 *  Define them according to your needs.  
 */

/*---------------------------------------------------------------------------*/
#define precision double  /* Precision of arithmetic (float or double)       */
                          /*                                                  */
#define V_DIFF_ON 1       /* if 1 then v-field diffuses                       */
#define GRAPHICS  0       /* if 1 then compile X11 graphics                   */
#define NINEPOINT 1       /* if 1 then 9 point Laplacian formulas             */
#define NBC       1       /* if 1 then Neumann BCs, if 0 periodic BCs         */
#define POLAR     1       /* use Polar coords                                 */
#define POLARPLOT 1       /* 1: plot (x,y)        0: plot (r,phi)             */
#define DEBUG     0
/*---------------------------------------------------------------------------*/

#if POLAR /* for polar coords pbc and ninepoint laplacian are nonsense */
#undef  NBC
#define NBC 1
#undef  NINEPOINT
#define NINEPOINT 0
#else     /* cart coords */
#undef  POLARPLOT
#define POLARPLOT 0
#endif

#define min(a,b) ((a)>(b) ? (b) : (a))
#define max(a,b) ((a)>(b) ? (a) : (b))

#define PRNT_WARN(verbose,level,message) if(verbose>=level) printf(message);

/* Include Standard C Headers */
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <mex.h>

/* Prototypes */
void     step          ();
void     step_ini      ();
void     plot          ();
void     plot_ini      ();
int      keyboard_chk  ();
void     tip_finder    ();
