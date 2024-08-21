/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*                                                                           */
/*                                 EZGRAPH.C                                 */
/*                   X11-Graphics Subroutines for EZ-SPIRAL                  */
/*                                                                           */
/*                                Version 2                                  */
/*                                                                           */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* 
 * $Revision: 2.2 $
 * $Date: 1994/06/27 23:31:12 $
 */
/*---------------------------------------------------------------------------*/

#include "ezspiral.h"

#include "/Applications/MATLAB_R2021a.app/extern/include/mex.h"

#if GRAPHICS   

#include "ezgraph.h"

/* External variables */
extern precision  **u, **v, **u_lap[2],*tip_x, *tip_y, a_param;
extern int        M1, M2, k2, ntips, u_plot, v_plot, verbose;

void write_fc (char *filename); /* from ezspiral */

/*---------------------------------------------------------------------------*/

void plot()

     /* 
      *  Main plotting routine. 
      *  Filled rectangles of appropriate colors are drawn to pixmaps.
      *  The spiral tip path is also drawn to the pixmaps if desired.
      *  After all graphics are drawn, the pixmaps are copied to the 
      *  X-window.  This is similar to double buffering and gives better
      *  animation than drawing directly to the X-window.
      */
{
  int i,j, iset;
 
  if(!(u_plot||v_plot)) 
       return;  /* neither field is plotted so return */
   
  if(u_plot) {
    /*  Plotting the u-field.
     *  First set the entire u square to the quiescent state
     *  (blue for color plot, white for bw plot).
     *  Filled rectangle will only be drawn if u is non-quiescent.  
     *  This is faster than drawing filled rectangles at every grid point.
     */
    switch (plot_type) {
    case Color :
      XSetForeground(theDisplay, theGC, theBluePixel);
      break;
    case Mono :
      XSetForeground(theDisplay, theGC, theWhitePixel);
      break;
    }

#if POLARPLOT
    XSetForeground(theDisplay, theGC, theWhitePixel);
#endif

    XFillRectangle(theDisplay, theUPixmap, theGC, 0, 0, sqsize1, sqsize2);
    
    for(i=M1;i>=1;i--) {
      for(j=1;j<=M2;j++) {
	iset = setColor(u[i][j],'u');
	
#if POLARPLOT /* draw arc at every point */
	ix= (int) (sqsize1*(1.0-(float)i/(float)M1))/2;
	d = 2*PX(i); 
	XFillArc(theDisplay, theUPixmap, theGC, ix, ix, d, d, PY(j), DY);
#else /*  Draw filled rectangle only if setColor returns non-zero value, 
       *  (i.e. if u is non-quiescent).  
       */
	if(iset) {
	  XFillRectangle(theDisplay, theUPixmap, theGC, PX(i), PY(j), irsize1, irsize2);
	}
#endif
      }
    }
  }

  if(v_plot) {
    /*  Plotting the v-field.
     *  Set the color and draw filled rectangles.
     */ 
#if POLARPLOT
    XSetForeground(theDisplay, theGC, theWhitePixel);
#endif

    XFillRectangle(theDisplay, theVPixmap, theGC, 0, 0, sqsize1, sqsize2);

    for(i=M1;i>=1;i--) {
      for(j=1;j<=M2;j++) {
	setColor(v[i][j],'v');
#if POLARPLOT
	ix=(int) (sqsize1*(1.0-(float)i/(float)M1))/2;
	d = 2*PX(i); 
	XFillArc(theDisplay, theVPixmap, theGC, ix, ix, d, d, PY(j), DY);
#else
	XFillRectangle(theDisplay, theVPixmap, theGC, 
		       PX(i), PY(j), irsize1, irsize2);
#endif      
      }
    }
  }

  /* If tip path is to be plotted, then do it here. */
  if(ntips) plot_tip_path();

  /*  Copy pixmaps to the window so that they become visible.
   *  The u_origin and v_origin determine where in the window the
   *  pixmaps are copied.  See plot_ini().
   */
  if(u_plot) {
    XCopyArea(theDisplay, theUPixmap, theWindow, theGC, 
	      0, 0, sqsize1, sqsize2, u_origin[0], u_origin[1]);
  }
    
  if(v_plot) {
    XCopyArea(theDisplay, theVPixmap, theWindow, theGC, 
	      0, 0, sqsize1, sqsize2, v_origin[0], v_origin[1]);
  }
    
  XFlush(theDisplay);
}
/*---------------------------------------------------------------------------*/

void save_image(int iter)

  /* Saves window image to file */
{
  static char filename[80], xpmfix[40];
  int filenumber;

  if (iter==0) {  sprintf(filename,"screen.xpm"); }
  else { filenumber = iter + 1000;
         sprintf(filename, "screen%4d.xpm", filenumber);
  }

  XpmWriteFileFromPixmap(theDisplay,filename,theVPixmap,NULL,NULL);
  sprintf(xpmfix,"@xpmfix ");
  strcat(xpmfix,filename);
  system(xpmfix);
}

/*---------------------------------------------------------------------------*/

int setColor (precision val, char cval)

     /* 
      *  Sets the foreground color.
      *  val is the value of the field of "type" cval (cval = 'u' or 'v').
      *  Returns 1 if cval == 'u' and the state is non-quiescent.
      */
{
  float scaled_v;
  int icol;

  switch (cval) {

  case 'u' :  /* Set the u color */
    if (val > 0.9) {
      /* Excited state: u > 0.9 */
      switch (plot_type) {
      case Color :
	/* For Color plot, red denotes excited. */
	XSetForeground(theDisplay, theGC, theRedPixel);
	break;
      case Mono :
	/* For Mono plot, excited states are black. */
	XSetForeground(theDisplay, theGC, theBlackPixel);
	break;
      }
      return(1);
    }
    else if(val > 0.1) {
      /* Interface: 0.1 < u < 0.9 
       * For both Color and Mono the interface states are black. */
      XSetForeground(theDisplay, theGC, theBlackPixel);
      return(1);
    }
    else {
      /* quiescent: u < 0.1 
       * Return 0 so that no rectangle is drawn. */ 
#if POLARPLOT
      XSetForeground(theDisplay, theGC, theBluePixel);
#endif
      return(0);
    }

  case 'v' :  /* Set the v_color */

    /* Scale the v-field to be between 0 and 1.
     * The maximum (VMAX) and minimum (VMIN) values used for
     * scaling are somewhat ad hoc.  The max and min functions
     * are used to ensure  0 <= scaled_v < 1.
     */
#define VMAX (0.8*a_param) 
#define VMIN 0.1 
    scaled_v = (val-VMIN)/(VMAX-VMIN);
    scaled_v = max(0.0,  scaled_v);
    scaled_v = min(0.999,scaled_v);

    switch (plot_type) {
    case Color :
      /* For color plot */
      icol = (int)maxColors*scaled_v;
      XSetForeground(theDisplay, theGC, theColors[icol] );
      break;
    case Mono :
      /* For Mono plot */
      if(scaled_v > 0.5) {
	XSetForeground(theDisplay, theGC, theBlackPixel);
      }
      else {
	XSetForeground(theDisplay, theGC, theWhitePixel);
      }
      break;
    }
  }
  return(0);
}
/*---------------------------------------------------------------------------*/

void plot_tip_path()

     /* 
      *  Plots the stored tip path.
      */
{
  XPoint *thePoints;
  int i;

  if(ntips<1) return;

  /* Set the color of the tip path */
  switch (plot_type) {
  case Color :
    XSetForeground(theDisplay, theGC, theWhitePixel);
    break;
  case Mono :
    XSetForeground(theDisplay, theGC, theBlackPixel);
    break;
  }

  /* Allocate space for thePoints */
  thePoints=(XPoint *)malloc((unsigned)(ntips)*sizeof(XPoint));

  /* Copy tip location to thePoints */
  for(i=0;i<ntips;i++) {
#if POLARPLOT 
    thePoints[i].x = (short) (sqsize1/2+rsize1*(tip_x[ntips-i-1]-1)*cos(2*M_PI*(tip_y[ntips-i-1]-1)/M2));
    thePoints[i].y = (short) (sqsize1/2-rsize1*(tip_x[ntips-i-1]-1)*sin(2*M_PI*(tip_y[ntips-i-1]-1)/M2));
#else
    thePoints[i].x = (short) PX(tip_x[ntips-i-1]);
    thePoints[i].y = (short) PY(tip_y[ntips-i-1]);
#endif
  }


  /* Draw the tip path */
  if(u_plot) {
    XDrawLines(theDisplay, theUPixmap, theGC, thePoints, ntips, CoordModeOrigin);     /* draw tip path */
  }

  if(v_plot) {
    XDrawLines(theDisplay, theVPixmap, theGC, thePoints, ntips, CoordModeOrigin);
  }

  free(thePoints);
}
/*---------------------------------------------------------------------------*/

void plot_ini()

     /* 
      *  Initializes the graphics. 
      */
{
  int width, height;

  if(!(u_plot||v_plot)) return;  /* neither field is plotted so return */

  /* Initialize X */
  initX();

  /*  Set up the graphics coordinates.
   *  sqsize is the length, in pixels, of one side of the square area in 
   *  which the u,v fields are plotted.  PLOTSIZE is the fraction of the
   *  full screen height for one of the square plots and thePHeight the 
   *  screen height in pixels.
   *  
   *  width and height are the width and height of the X-window in pixels.
   *  The height is always sqsize + 2*BORDER (to account for BORDER pixels
   *  above and BORDER pixels below the square).
   *  In the case of both fields plotted, the width is 2*sqsize (for the
   *  2 fields) + 3*BORDER (right BORDER, left BORDER and a BORDER 
   *  between the two plots).
   *  In the case of one field plotted, the width is sqsize + 2*BORDER 
   *  exactly as for the height.
   *
   *  The arrays u_origin and v_origin control where the u and v squares 
   *  appear within the X-window.  See also plot().  
   *  In the case of both fields, the upper-left corner of the u-field 
   *  is copied to the window BORDER pixels down and BORDER pixels to 
   *  the right of the upper-left corner of the window.  The upper-left 
   *  corner of the v-field is copied to the window BORDER pixels down and 
   *  sqsize+2*BORDER pixels to the right of the upper-left corner of the 
   *  window.  
   *  In the case of one field, the upper-left corner of that field
   *  is copied to the window BORDER pixels down and BORDER pixels to 
   *  the right of the upper-left corner of the window.  Both u_origin
   *  and v_origin are set but only one is used.
   */

  sqsize1 = PLOTSIZE*thePHeight; 
  sqsize2 = PLOTSIZE*thePHeight;

  if(u_plot && v_plot) {
    /* 
     * Both u and v fields plotted.
     */
    height =   sqsize2 + 2*BORDER;
    width  = 2*sqsize1 + 3*BORDER;
    u_origin[0] = BORDER;
    u_origin[1] = BORDER;
    v_origin[0] = sqsize1 + 2*BORDER;
    v_origin[1] = BORDER;
  }
  else {
    /* 
     * Only one field plotted. 
     */
    height = sqsize2 + 2*BORDER;
    width = sqsize1 + 2*BORDER;
    u_origin[0] = BORDER;
    u_origin[1] = BORDER;
    v_origin[0] = BORDER;
    v_origin[1] = BORDER;
  }

  /* rsize = the "exact" (i.e. floating point) size of the rectangles drawn 
   * at each grid point.  irsize is the integer equivalent, with the +1 
   * added to ensure that rectangles are overlapping.
   */
  
#if POLARPLOT
  rsize1  = sqsize1/(2*M1); /* 2R in Window */
  rsize2 = rsize1;
  irsize1 = (int)(rsize1);
  irsize2 = irsize1;
#else
  rsize1  = (float)sqsize1/(float)M1;
  rsize2 = (float)sqsize2/(float)M2;

  irsize1 = (int)(rsize1)+1;
  irsize2 = (int)(rsize2)+1;
#endif

  /* Open a window and call other necessary X routines */
  openWindow(WINX, WINY, width, height);

  /* If Color plot, then initialize colors */
  if(plot_type==Color) initColors();
}
/*---------------------------------------------------------------------------*/

void initX()

     /* 
      *  Initializes X.
      */
{
  char *theDisplayName = NULL;
  
  /* Open the display */
  if( (theDisplay = XOpenDisplay(NULL)) == NULL) {
    fprintf(stderr, 
	    "ERROR: Could not open a connection to X on display %s\n",
	    XDisplayName(theDisplayName));
    exit(1);
  }

  theScreen     = DefaultScreen(theDisplay);
  theDepth      = DefaultDepth (theDisplay, theScreen);
  thePWidth     = DisplayWidth (theDisplay, theScreen);
  thePHeight    = DisplayHeight(theDisplay, theScreen);
  theBlackPixel = BlackPixel   (theDisplay, theScreen);
  theWhitePixel = WhitePixel   (theDisplay, theScreen);
  theColormap   = DefaultColormap(theDisplay, theScreen);

  if(theDepth==1) {
    plot_type = Mono;
  }
  else {
    plot_type = Color;
  }

  /* Print useful information. I suggest printing this at this at least once */
  if(verbose>1) {
    printf("%s version %d of the X Window System, X%d R%d\n",
	   ServerVendor    (theDisplay),
	   VendorRelease   (theDisplay),
	   ProtocolVersion (theDisplay),
	   ProtocolRevision(theDisplay));

    if(theDepth==1) {
      printf("Color plane depth...........%d (monochrome)\n", theDepth);
    }
    else {
      printf("Color plane depth...........%d \n",             theDepth);
    }

    printf("Display Width...............%d \n", thePWidth);
    printf("Display Height..............%d \n", thePHeight);
    printf("The display %s\n", XDisplayName(theDisplayName));
  }
}
/*---------------------------------------------------------------------------*/

void openWindow(int winx, int winy, int width, int height)

     /* 
      *  Creates and maps the X window.
      *  Also creates the Graphics Context and Pixmaps.
      */
{
  XSizeHints    theSizeHints;
  XTextProperty theWindowName, theIconName;
  char          *Name="EZ-Spiral";
  unsigned long background;

  /* Set the background color (i.e. color of border around the square plots */
  background = BlackPixel(theDisplay,theScreen);
  if(BACKGROUND) background = WhitePixel(theDisplay,theScreen);

  /* Create the window */
  theWindow = XCreateSimpleWindow(theDisplay,
				  RootWindow(theDisplay, theScreen),
				  winx, winy, width, height, 0,
				  background, background);  

  /* Set standard window properties. 
   * theWindowName and theIconName are set to Name, which unless you
   * change it, it will be EZ-Spiral.
   * theSizeHints.flags are set such that if winx and winy both positive, 
   * they try to locate the window there, else let the window manager decide.
   */
  XStringListToTextProperty(&Name,1,&theWindowName);
  XStringListToTextProperty(&Name,1,&theIconName);

  if(winx>0 && winy>0) {
    theSizeHints.flags  = USPosition | PSize;
  } else {
    theSizeHints.flags  = PPosition | PSize;
  }

  XSetWMProperties(theDisplay, theWindow, &theWindowName, &theIconName,
		   NULL, 0, &theSizeHints, NULL, NULL);

  /* This selects which events the program will see.  This is used, e.g. to 
   * pass keyboard or mouse input to the program.  The event types are set 
   * in EV_MASK.  (See ezgraph.h)
   */
  XSelectInput(theDisplay, theWindow, EV_MASK);

  
  /* Create the Graphics Context. */
  theGC = XCreateGC(theDisplay, theWindow, (unsigned long) 0, NULL);
  
  
  if(theGC==0) {
    XDestroyWindow(theDisplay,theWindow);
    printf("Could not create the GC \n");
    exit(1);
  }

  /* Now that the Graphics Context has been created, setting and changing
   * most graphics characteristics (such as line width or line type) can 
   * easily be done with appropriate X routine calls.  I include an example 
   * of how to set the line width to 2 pixels.  LineSolid has the obvious
   * meaning.  CapRound and JoinRound are probably not important to you.
   *
   *  XSetLineAttributes(theDisplay, theGC, (unsigned int) 2, 
   *                     LineSolid, CapRound, JoinRound);
   */

   XSetLineAttributes(theDisplay, theGC, (unsigned int) 2, LineSolid, CapRound, JoinRound);
  /*  XSetArcMode(theDisplay, theGC, ArcChord); */

  /* Map the window, i.e. make it visible */
  XMapWindow(theDisplay, theWindow); 
  XFlush(theDisplay);

  /* Create pixmaps of size sqsize by sqsize for the u/v plots */
  if(u_plot) {
    theUPixmap = XCreatePixmap(theDisplay,theWindow,
			       sqsize1, sqsize2, theDepth);
    if(theUPixmap==False) printf("couldn't open UPixmap\n");
  }

  if(v_plot) {
    theVPixmap = XCreatePixmap(theDisplay,theWindow,
			       sqsize1, sqsize2, theDepth);
    if(theVPixmap==False) printf("couldn't open VPixmap\n");
  }
}
/*---------------------------------------------------------------------------*/

void initColors()

{
  XColor theRGBColor, theHardwareColor;
  int    theStatus;
  int    cv,i;

  /* Colors for the u-field: Red, Blue, (and black already defined) */
  
  if( XAllocNamedColor(theDisplay, theColormap, "Red",
		       &theRGBColor, &theHardwareColor) ) {
    theRedPixel = theHardwareColor.pixel; 
  }
  else {
    printf("Cannot allocate the Red Pixel \n");
  }
  
  if( XAllocNamedColor(theDisplay, theColormap, "Blue",
		       &theRGBColor, &theHardwareColor) ) {
    theBluePixel = theHardwareColor.pixel; 
  }
  else {
    printf("Cannot allocate the Blue Pixel \n");
  }
  
  /*  Colors for the v-field: Continuous range from Red to Blue.
   *  A few possibilities are included; you should modify according
   *  to taste.
   */
  
  for(i=0;i<maxColors;i++) {
    
    cv = (int)( i*65535./(maxColors-1.) );
#if 1
    theRGBColor.red   = theHardwareColor.red   = 65535 - cv;
    theRGBColor.green = theHardwareColor.green = cv/1.6; 
    theRGBColor.blue  = theHardwareColor.blue  = cv;
#else  /* Alternative color map. Try you own. */
    theRGBColor.red   = theHardwareColor.red   = 65535 - cv;
    theRGBColor.green = theHardwareColor.green = cv; 
    theRGBColor.blue  = theHardwareColor.blue  = cv;
#endif
    
    theStatus = XAllocColor(theDisplay, theColormap, 
			    &theHardwareColor);
    if(theStatus) {
      theColors[i] = theHardwareColor.pixel;
      if(verbose>2) {
	printf("Pixel and RGBvalues (requested/received) = ");
	printf("%d, %d, %d, %d, %d, %d, %d \n", 
	       (int)theHardwareColor.pixel, 
	       (int)theRGBColor.red, 
	       (int)theRGBColor.green, 
	       (int)theRGBColor.blue,
	       (int)theHardwareColor.red, 
	       (int)theHardwareColor.green, 
	       (int)theHardwareColor.blue);
      }
    }
    else {
      printf("Failed ! Pixel and RGBvalues (requested) = %d, %d, %d, %d \n", 
	     (int)theHardwareColor.pixel, 
	     (int)theRGBColor.red, 
	     (int)theRGBColor.green, 
	     (int)theRGBColor.blue);
    }
  }
}
/*---------------------------------------------------------------------------*/

void quitX()

     /* 
      *  Free memory and close the theDisplay.
      */
{
  XFreeGC(theDisplay,theGC);
  if(u_plot) XFreePixmap(theDisplay,theUPixmap);
  if(v_plot) XFreePixmap(theDisplay,theVPixmap);
  XCloseDisplay(theDisplay);
}
/*---------------------------------------------------------------------------*/

int keyboard_chk()

     /* 
      *  Checks keyboard event queue and takes appropriate action. 
      */
{
  XEvent         theEvent;
  KeySym         theKeySym;
  int            theKeyBufferMaxLen = 64;
  char           theKeyBuffer[65];
  XComposeStatus theComposeStatus;
  int            n_up=0, n_rt=0;

  if(!(u_plot||v_plot)) return(0);   /* neither field is plotted so return 0 */

  while(XCheckWindowEvent(theDisplay, theWindow, EV_MASK, &theEvent) ) {
    /* If an event then ... */

    switch(theEvent.type) {

#if FOCUS_CK
      /*  This case is necessary for window managers which do not set 
       *  keyboard focus to theWindow when the pointer enters theWindow.
       *  See also ezgraph.h.
       */
    case EnterNotify:
      XSetInputFocus(theDisplay, theWindow, RevertToPointerRoot, CurrentTime); 
      break;
#endif
      
    case KeyPress:
      /* A KeyPress event has occurred.  Take appropriate action */

      XLookupString((XKeyEvent *)&theEvent,   theKeyBuffer, theKeyBufferMaxLen,
		    &theKeySym, &theComposeStatus);
      switch(theKeySym) {
      case XK_Up:     n_up++; break;                           /* move up    */
      case XK_Down:   n_up--; break;                           /* move down  */
      case XK_Right:  n_rt++; break;                           /* move right */
      case XK_Left:   n_rt--; break;                           /* move left  */
      case XK_Escape: quitX(); return(1);                           /* quit  */
      case XK_Q:      quitX(); return(1);                           /* quit  */
      case XK_q:      quitX(); return(1);                           /* quit  */
      case XK_P:                                                    /* pause */  
	XWindowEvent(theDisplay, theWindow, EV_MASK, &theEvent);
	break;      
      case XK_p:                                                    /* pause */  
	XWindowEvent(theDisplay, theWindow, EV_MASK, &theEvent);
	break;      
      }
      break;
    }
  }

  /* If a move has been requested, call mover() */
  if(n_up||n_rt) mover(n_up, n_rt);  
  return(0);
}
/*---------------------------------------------------------------------------*/

void mover (int n_up, int n_rt)

     /* 
      *  Moves the spiral(s) in response to a key press. 
      *  To save time, v_lap is not moved.  I have not implemented "wrap 
      *  around" moving in the case of periodic boundary conditions.
      */
{
  int n_move,i,j;
  
  if( n_up > 0 ) {                                                     /* up */
    n_move = n_up*MOVESTEP2;
    for(j=M2;j>n_move;j--) {
      for(i=1;i<=M1;i++) {
	u[i][j] = u[i][j-n_move]; 
	v[i][j] = v[i][j-n_move];
	u_lap[k2][i][j] = u_lap[k2][i][j-n_move]; 
      }
    }
  } 
  if( n_up < 0 ) {                                                   /* down */
    n_move = n_up*MOVESTEP2;
    for(j=1;j<=M2+n_move;j++) {
      for(i=1;i<=M1;i++) {
	u[i][j] = u[i][j-n_move]; 
	v[i][j] = v[i][j-n_move];
	u_lap[k2][i][j] = u_lap[k2][i][j-n_move]; 
      }
    }
  } 
  if( n_rt > 0 ) {                                                  /* right */
    n_move = n_rt*MOVESTEP1;
    for(i=M1;i>n_move;i--) {
      for(j=1;j<=M2;j++) {
	u[i][j] = u[i-n_move][j]; 
	v[i][j] = v[i-n_move][j];
	u_lap[k2][i][j] = u_lap[k2][i-n_move][j]; 
      }
    }
  }
  if( n_rt < 0 ) {                                                   /* left */
    n_move = n_rt*MOVESTEP1;
    for(i=1;i<=M1+n_move;i++) {
      for(j=1;j<=M2;j++) {
	u[i][j] = u[i-n_move][j]; 
	v[i][j] = v[i-n_move][j];
	u_lap[k2][i][j] = u_lap[k2][i-n_move][j]; 
      }
    }
  } 

  ntips = 0;
  plot();
}
/*---------------------------------------------------------------------------*/

#else   /* else #if GRAPHICS */

/* if NOT GRAPHICS, define these dummy routines. */

void plot_ini    () {}
void plot        () {}
int  keyboard_chk() {return(0);}

#endif

/*---------------------------------------------------------------------------*/

