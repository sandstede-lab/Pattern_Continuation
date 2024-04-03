/* X11 Headers */
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/keysym.h>
#include <X11/keysymdef.h>
#include <X11/xpm.h>
/*---------------------------------------------------------------------------*/

/*  Here are some macro definitions which can be set as you desire.  If you 
 *  prefer, these quantities can be changed to variables and set at run time. 
 */

#define PLOTSIZE     0.35    /* plot size as fraction of the screen height   */
#define WINX         20      /* window location, in pixels, from screen left */
#define WINY         20      /* window location, in pixels, from screen top  */
			     /* if either WINX or WINY < 0 then window       */
                             /* location is determined by the window manager */
#define BORDER       10      /* border width in pixels                       */
#define BACKGROUND   1       /* background color (0 BLACK or 1 WHITE)        */
/*---------------------------------------------------------------------------*/

#define maxColors 32
#define Color     1  /* Do not change this */
#define Mono      0  /* Do not change this */
#define FOCUS_CK  0

#if FOCUS_CK
/*  This case is necessary for window managers which do not set 
 *  keyboard focus to theWindow when the pointer enters theWindow.
 *  I have found FOCUS_CK should be 1 on SUN systems, but 0 on others. 
 *  0 is preferable because if FOCUS_CK is 1 then moving the pointer 
 *  into the graphics window will cause a paused simulation to continue 
 *  (slightly annoying).  See also keyboard_chk() in ezgraph.c.
 */
#  define EV_MASK (KeyPressMask | EnterWindowMask)
#else
#  define EV_MASK (KeyPressMask)
#endif

/* number of grid points moved per key stroke */
#define MOVESTEP1 ((int)(1+M1/50)) 
#define MOVESTEP2 ((int)(1+M2/50)) 
/*  Macros for converting from grid coordinates (i,j) to 
 *  pixel location (PX,PY) within one of the square plots.
 */
#if POLARPLOT
#define DY ( (int) (360*64/M2) )
#define PX(i) ( (int) (sqsize1*(i)/(2*M1)) )
#define PY(j) ( (j-1)*DY ) 
#define CARTX(i) ((int)(sqsize1/2+rsize1*((i)-1)))
#define CARTY(i) ((int)(sqsize1/2-rsize1*((i)-1)))
#else 
#define PX(i) ((int)(rsize1*((i)-1)))      
#define PY(j) ((int)(rsize2*(M2-(j))))   
#endif
/* Global variables within ezgraph.c */
static int            u_origin[2], v_origin[2];
static int            sqsize1, sqsize2, irsize1, irsize2, plot_type;
static float          rsize1, rsize2;
static Display       *theDisplay; 
static Window         theWindow;
static Pixmap         theUPixmap, theVPixmap; 
static GC             theGC;
static Colormap       theColormap;
static int            theScreen; 
static int            theDepth;
static int            thePWidth;
static int            thePHeight;
static unsigned long  theBlackPixel;
static unsigned long  theWhitePixel;
static unsigned long  theRedPixel;
static unsigned long  theBluePixel;
static unsigned long  theColors[maxColors];

#if POLARPLOT
static int            d;
static int            ix;
#endif

/* Functions local to ezgraph.c */
static int    setColor      (precision val, char cval);
static void   plot_tip_path ();
static void   initX         ();
static void   openWindow    (int x, int y, int width, int height);
static void   initColors    ();
static void   quitX         ();
static void   mover         (int n_up, int n_rt);


