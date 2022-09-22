/*
  CSC D18 - Assignment 1 - 2D light propagation
  
  Have a look below and be sure you understand the
  different data structures and functions provided
  to help you complete your work.
  
  By F. J. Estrada, Aug. 2017
*/

#ifndef __light2D_header
#define __light2D_header

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#define PI 3.1415926535
#define TOL .000001	
#define MAX_OBJECTS 25

/*****************************************************************************************************
* Useful Data Structures
* ***************************************************************************************************/
struct point2D{					// Points and vectors in 2D space
  double px;
  double py;
};

struct ray2D{					// Rays in 2D
  struct point2D p;	// Ray origin
  struct point2D d;	// Ray direction
  double R,G,B;		// Colour of this light ray
  int inside_out;	// Toggle flag
			//   0 - ray is travelling in air (outside an object)
			//   1 - ray is travelling inside an object
  int monochromatic;	// Flag to indicate whether this ray is monochromatic (true if set to 1)
  double H;		// For monochromatic rays, HUE value (used to obtain colour, and as
			//   a convenient substitute for wavelength) values in [0 1] go
			//   from deep red to purple
};

struct circ2D{					// Simple 2D objects (circles in this case)
  struct point2D c;	// Center location
  double r;		// Circle radius
  int material_type;	// Type of object:
			//   0 - Mirror
			//   1 - Scattering
			//   2 - Refracting (transparent)
  double r_idx;		// Index of refraction for transparent materials
};

struct wall2D{					// Walls in 2D
  struct ray2D w;	// The actual wall is just a ray, the origin is one endpoint, and origin + d 
			// gives the second endpoind (i.e. with lambda=1)
  int material_type;	// Type of material, can only be 0 or 1 with the same meaning as for circles
};


struct light2D{					// Basic Light Sources
  struct ray2D l;	// Contains the location and if needed the direction of the
			// lightsource
  int light_type;	// Type of light source
			// 0 - Point. Only the location matters, ignore the direction
			// 1 - Laser. Emits light exclusively along a single direction
			//            from its location.
  double R,G,B;		// Colour of this light source
};

/***********************************************************************************************************
* GLOBAL DATA (yes, it's actually useful!)
************************************************************************************************************/
struct circ2D obj_list;			// This will hold the list of objects in the scene
double *imRGB;				// Pointer to the array where the image will be rendered
unsigned char *im;			// Final image for storage
int sx,sy;				// Image resolution
int max_depth;				// Maximum recursion depth
int num_rays;				// Number of light-source rays to propagate
struct circ2D objects[MAX_OBJECTS];	// Array to hold scene objects (up to MAX_OBJECTS)
struct wall2D walls[4];			// Four walls stored in order TOP, RIGHT, BOTTOM, LEFT
struct light2D lightsource;		// The light source for the scene
double w_x,w_y;


/************************************************************************************************************
* Function definitions
*************************************************************************************************************/
void hue2RGB(double H, double *R, double *G, double *B);		// Returns colour for a given hue
void renderObjects(void);						// Renders objects onto image
// This function renders a ray onto the image between the two given points, with the specified RGB colour
void renderRay(struct point2D *p1, struct point2D *p2, double R, double G, double B);
void setPixel(double x, double y, double R, double G, double B);	// Set image pixel at (x,y) to the 
									// specified colour.
double dot(struct point2D *p, struct point2D *q);	// Dot product between vectors p and q 
void normalize(struct point2D *d);			// Normalize a vector to unit length
void addCirc(struct point2D *c, double r, int type, double n);	// Inserts a circle into the objects array
void buildWalls(void);		// This function makes the walls of the box. Implemented in buildScene.c
void buildScene(void);		// This creates the list of objects and sets up the
				// light source. Implemented in buildScene.c, so you can
				// change the scene by updating buildScene.c accordingly
int main(int argc, char *argv[]);

#endif

