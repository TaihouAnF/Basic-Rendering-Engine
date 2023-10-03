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
 walls[0].material_type=0;		// Change to 0 to make the walls mirror-reflecting!
 
 p.px=W_RIGHT;
 p.py=W_TOP;
 d.px=0;
 d.py=W_BOTTOM-W_TOP;
 walls[1].w.p=p;
 walls[1].w.d=d;
 walls[1].material_type=0;

 p.px=W_RIGHT;
 p.py=W_BOTTOM;
 d.px=W_LEFT-W_RIGHT;
 d.py=0;
 walls[2].w.p=p;
 walls[2].w.d=d;
 walls[2].material_type=0;
 
 p.px=W_LEFT;
 p.py=W_BOTTOM;
 d.px=0;
 d.py=W_TOP-W_BOTTOM;
 walls[3].w.p=p;
 walls[3].w.d=d;
 walls[3].material_type=0;
}

void buildScene(void) {
  struct point2D c,p,d;
  struct ray2D l;
  // make a circle full of circles
  double radius = 1;
  double radius2 = 1.5;
  double mini_radius = 0.2;
  double mini_radius2 = 0.25;
  double increment = (2 * PI)/12;
  int mattype = 0;
  for (double i=0; i < 2*PI; i += increment) {
    c.px = radius * cos(i);
    c.py = radius * sin(i);
    addCirc(&c, mini_radius, mattype, 2.74);
  }
  mattype = 2;
  for (double i=0; i < 2*PI; i += increment) {
    c.px = radius2 * cos(i);
    c.py = radius2 * sin(i);
    addCirc(&c, mini_radius2, mattype, 2.74);
  }
  p.px=0;
  p.py=0;
  d.px=0.4720;
  d.py=-0.8816;
  normalize(&d);
  l.p=p;
  l.d=d;
  lightsource.l=l;
  lightsource.light_type=0;
  lightsource.R=1.0;
  lightsource.G=1.0;
  lightsource.B=1.0;
}