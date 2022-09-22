/*
  CSC D18 - Assignment 1 - 2D light propagation
  
  This function specifiees the content of the scene. I'm
  providing a very simple scene to test your code. Change
  it as needed to test the different components of the
  assignment.

  Once you're happy your code works, feel free to modify
  this function to create a more complex scene and render
  a cool picture!

  This function also defines the size of the box
  
  Starter by: F.J. Estrada, Aug. 2017
*/

// Coordinates of the corners of the box where our scene is
// defined
#define W_TOP -2.0
#define W_LEFT -2.0
#define W_RIGHT 2.0
#define W_BOTTOM 2.0

void buildWalls(void)
{
 // This sets the 4 walls of our box given the values of the box corners
 // defined above.
 // E.g. with the original values of W_TOP=-2.0, W_LEFT=-2.0,
 //      W_RIGHT=2.0, W_BOTTOM=2.0, you get
 //
 //   (-2,-2)                              (2,-2)
 //     *------------------------------------>*
 //     ^                                     |
 //     |                                     |
 //     |                                     |
 //     |                 *                   |
 //     |               (0,0)                 |
 //     |                                     |
 //     |                                     |
 //     *<------------------------------------*
 //   (-2,2)                               (2,2)
 //
 // Notice that I've defined object space so that the Y coordinate increases 
 // downward, this is so that it agrees with the direction of increasing Y in
 // the image.
 //
 // Each 'wall' is really a ray, the origin p is one of the endpoints of the wall,
 // and the direction is scaled so that p + lambda*d gives you the second
 // endpoint for lambda=1. The 4 walls make up a closed loop and the directions
 // are as shown above. This will be imporant for you when you need to determine
 // the normal vector for each wall - you must make sure the vector points toward
 // the *inside* of the box! 
 //
 // You must also make sure that all objects you define for the scene fall within
 // the window defined by W_TOP, W_LEFT, W_RIGHT, and W_BOTTOM.

 struct point2D p,d;
 
 p.px=W_LEFT;
 p.py=W_TOP;
 d.px=W_RIGHT-W_LEFT;
 d.py=0;
 walls[0].w.p=p;
 walls[0].w.d=d;
 walls[0].material_type=1;		// Change to 0 to make the walls mirror-reflecting!
 
 p.px=W_RIGHT;
 p.py=W_TOP;
 d.px=0;
 d.py=W_BOTTOM-W_TOP;
 walls[1].w.p=p;
 walls[1].w.d=d;
 walls[1].material_type=1;

 p.px=W_RIGHT;
 p.py=W_BOTTOM;
 d.px=W_LEFT-W_RIGHT;
 d.py=0;
 walls[2].w.p=p;
 walls[2].w.d=d;
 walls[2].material_type=1;
 
 p.px=W_LEFT;
 p.py=W_BOTTOM;
 d.px=0;
 d.py=W_TOP-W_BOTTOM;
 walls[3].w.p=p;
 walls[3].w.d=d;
 walls[3].material_type=1;
}

void buildScene(void)
{
 struct point2D c,p,d;
 struct ray2D l;

 // Two spheres and one lightsource
 c.px=1.5;
 c.py=-1.5;
 addCirc(&c,.55,1,1.6);

 c.px=.5;
 c.py=-.5;
 addCirc(&c,.25,0,2.47);

 p.px=0.0;			// Set a light source at the origin, pointing in direction
 p.py=0.0;			// (1.0, -.8)  (remember that Y increases downward!)
 d.px=1.0;			// This ray should hit both spheres, so you can test your
 d.py=-.8;			// intersection code. The closest sphere is reflecting, so
 normalize(&d);			// you can test your code for mirror-reflecting rays
 l.p=p;
 l.d=d;
 lightsource.l=l;
 lightsource.light_type=1;	// Change this to set the light source type to point or laser 
 lightsource.R=1.0;		// Lightsource colour in RGB space. All values must be in [0,1]
 lightsource.G=1.0;		// and (1,1,1) is white
 lightsource.B=1.0;
}
