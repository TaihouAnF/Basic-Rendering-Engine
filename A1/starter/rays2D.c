/*
  CSC D18 - Assignment 1 - 2D light propagation

  This is the place where you will be doing most of your work for solving this
  assignment. Before you start working here, you shold have become familiar
  with the data structures and functions available to you from light2D.h, and
  with the way a scene is built in buildScene.c

  Go over each of the functions below, and implement the different components
  of the solution in the sections marked

  /************************
  / TO DO:
  ************************ /

  Do not add or modify code outside these sections.

  Details about what needs to be implemented are described in the comments, as
  well as in the assignment handout. You must read both carefully. 

  Starter by: F.J. Estrada, Aug. 2017
*/

/****************************************************************************
 * Uncomment the #define below to enable debug code, add whatever you need
 * to help you debug your program between #ifdef - #endif blocks
 * ************************************************************************/
#define __DEBUG_MODE

/*****************************************************************************
* COMPLETE THIS TEXT BOX:
*
* 1) Student Name:		Anson Feng
* 2) Student Name:		Chou Lu Wei
*
* 1) Student number:  1004721955
* 2) Student number:  1004295693
* 
* 1) UtorID           fengdian
* 2) UtorID           choulu1
* 
* We hereby certify that the work contained here is our own
*
* ___Anson Feng_______             _____Choulu1______
* (sign with your name)            (sign with your name)
********************************************************************************/

struct ray2D makeLightSourceRay(void)
{
 /*
   This function should return a light ray that has its origin at the light
   source, and whose direction depends on the type of light source.

   For point light sources (which emit light in all directions) the direction
    has to be chosen randomly and uniformly over a unit circle (i.e. any
    direction is equally likely)

   For a laser type light source, the direction of the ray is exactly the same
    as the direction of the lightsource.
    
   Set the colour of the ray to the same colour as the light source, and
    set the inside_outside flag to 0 (we assume light sources are 
    outside objects)

   In either case, the direction vector *must be unit length*.
*/
 
 /************************************************************************
 *  TO DO: Complete this function so that we can sample rays from the
 *         lightsource for the light propagation process.
 ************************************************************************/
 
 struct ray2D ray = {0};

 // This creates a dummy ray (which won't go anywhere since the direction
 // vector d=[0 0]!. But it's here so you see which data values you have to
 // provide values for given the light source position, and type.  
 // ** REPLACE THE CODE BELOW with code that provides a valid ray that
 //    is consistent with the lightsource.
 
 ray.p.px=lightsource.l.p.px;			// Ray's origin(old)
 ray.p.py=lightsource.l.p.py;     // Changed its origin to lightsource.location

if (lightsource.light_type) {    // Ray's direction: if laser, then we set direction as lightsource;
//ray.d.px=lightsource.l.d.px;	// else, change direction to?
//ray.d.py=lightsource.l.d.py;
ray.d.px=0.711451;
ray.d.py=-0.702736;
} else {
double rand_ang=(double)rand()/RAND_MAX*2*PI;   // Using rand()/RAND_MAX to get number from [0, 1].
ray.d.px=lightsource.l.p.px+cos(rand_ang);
ray.d.py=lightsource.l.p.py+sin(rand_ang);
}

 ray.inside_out=0;		// Initially 0 since the ray starts outside an object
 ray.monochromatic=0;		// Initially 0 since the ray is white (from lightsource)
 ray.R=lightsource.R;			// Ray colour in RGB must be the same as the lightsource
 ray.G=lightsource.G;
 ray.B=lightsource.B;
 
 return(ray);			// Currently this returns dummy ray(old), now it should return correct ray
}

void propagateRay(struct ray2D *ray, int depth)
{
 /*
   This function carries out the light propagation process. It is provided with access
   to a ray data structure, and must perform the following steps (in order!):
   
   - Check if maximum recursion depth has been reached (in which case, it just returns)
   - Find the *closest* intersection between the ray and objects in the scene. This
     means you have to check against the 4 walls, and any circles added in buildScene,
     determine which intersection is closest, and obtain the intersection point, the
     normal at the intersection, and the lambda value at which the intersection happens.
   - Renders the ray onto the image from its starting point all the way up to the 
     intersection point.
   - At the intersection, use the material properties to determine how the propagation
     process proceeds:
         * For mirror materials, compute the mirror reflection direction, create a ray
           along that direction whose origin is the intersection point, and propagate
           that ray
         * For scattering materials, choose a random direction within +- 90 degrees of
           the normal, create a ray with that direction with origin at the intersection
           point, and propagate that ray
         * For refracting materials you will need to propagate two rays - one in the
           mirror reflection direction (just like for reflecting materials), and
           another in the refraction direction. Propagate both rays one after the other!
           
   NOTE: You should only care about intersections for which lambda is POSITIVE (in front
         of the ray), and greater than zero (e.g. if the ray is reflected from some
         object, you do not care about the intersection with the object itself which will
         have a lambda of *very close* to zero)
    
   In every case, make sure the ray's direction vector has unit length. You will need to
   complete other functions as part of your work here.
*/
  
 /*********************************************************************************
  * TO DO: Complete this function to implement light propagation!
  ********************************************************************************/
 
 // Define your local variables here
printf("depth: %d\n", depth); 
if (depth>=max_depth) return;	 	// Leave this be, it makes sure you don't
					// recurse forever

//  float lambda1 = (ray->p.px - walls[3].w.p.px) / ray->d.px;
//  struct point2D intersectionPt = {
//    walls[3].w.p.px,
//    ray->p.py + lambda1 * ray->d.py
//  };
double lambda_min = 0.0;
for (int i = 0; i < 4; i++) {
  if (i == 0 || i == 2) { // down/up-horizontal walls, which y is fixed.
    if (ray->d.py != 0) {
      double lambda = (walls[i].w.p.py - ray->p.py)/ray->d.py;
      if ((lambda_min == 0.0 && lambda > 0) || (lambda_min != 0.0 && lambda < lambda_min && lambda > 0)) {
        lambda_min = lambda;  // only update when it is less than the current lambda_min
      }
    }                     // Else: The Ray parallel to this line.
  } else {                // left/right-vertical walls, which x is fixed.
    if (ray->d.px != 0) {
      double lambda = (walls[i].w.p.px - ray->p.px)/ray->d.px;
      if ((lambda_min == 0.0 && lambda > 0) || (lambda_min != 0.0 && lambda < lambda_min && lambda > 0)) {
        lambda_min = lambda;  // only update when it is less than the current lambda_min
      }
    }                     // Else: The Ray parallel to this line.
  }
}
 



 // Step 1 - Find *closest* intersection with the 4 walls (the written part of A1
 //          should help you figure out how to do that.

 // How many walls can the ray intersect? how many walls can the ray intersect in the
 // forward direction?

 // Step 2 - Check for intersection against objects in the object array - you must
 //          complete the intersectRay() function, call it, and obtain the closest
 //          intersection (in the forward ray direction) with objects in the scene.
 //          Note that you must provide variables for intersectRay() to return
 //          the point of intersection, normal at intersection, lambda, material type,
 //          and refraction index for the closest object hit by the ray.

struct point2D intersectPt;
struct point2D normal;
double lambda_intersect = lambda_min;
int material_type = 0;
double refraction_index = 0.0;
intersectRay(ray, &intersectPt, &normal, &lambda_intersect, &material_type, &refraction_index);
struct point2D rayIntersect;
rayIntersect.px = ray->p.px + lambda_intersect * ray->d.px;
rayIntersect.py = ray->p.py + lambda_intersect * ray->d.py;

// Step 4 - Render the ray onto the image. Use renderRay(). Provide renderRay() with
//          the origin of the ray, and the intersection point (it will then draw a
//          ray from the origin to the intersection). You also need to provide the
//          ray's colour.
renderRay(&ray->p, &rayIntersect, ray->R, ray->G, ray->B);
//struct point2D reflection;
//reflection.px = -0.969992;
//reflection.px = 0.243136;
//renderRay(&ray->p, &reflection, ray->R, ray->G, ray->B);
// draw the normal at the intersection point
// draw the normal at the intersection point debugging
if (lambda_intersect != lambda_min) {
  // struct point2D normalPt;
  // normalPt.px = intersectPt.px + normal.px;
  // normalPt.py = intersectPt.py + normal.py;
  // renderRay(&rayIntersect, &normalPt, 1.0, 0.0, 0.0);
  // This is for reflection
  if (material_type == 0) {
    struct ray2D refl_ray;
    refl_ray.p.px=intersectPt.px;
    refl_ray.p.py=intersectPt.py;
    refl_ray.R=ray->R;
    refl_ray.G=ray->G;
    refl_ray.B=ray->B;
    refl_ray.d.px=-2*dot(&normal, &ray->d)*normal.px+ray->d.px;
    refl_ray.d.py=-2*dot(&normal, &ray->d)*normal.py+ray->d.py;
    refl_ray.inside_out=ray->inside_out;
    refl_ray.H=ray->H;
    refl_ray.monochromatic=ray->monochromatic;
    //if (abs(refl_ray.d.px- ray->d.px) < 10 && abs(refl_ray.d.px- ray->d.px) < 10) {
    //  printf("reflected ray is the same direction as incoming ray at (%f, %f)\n", abs(refl_ray.d.px- ray->d.px), abs(refl_ray.d.py- ray->d.py));
    //}
    //printf("originating ray x:(%f) y:(%f) dx:(%f) dy=(%f)\n", ray->p.px, ray->p.py, ray->d.px, ray->d.py);
    //printf("reflected ray x:(%f) y:(%f) dx:(%f) dy=(%f)\n", refl_ray.p.px, refl_ray.p.py, refl_ray.d.px, refl_ray.d.py);
    propagateRay(&refl_ray, depth+1);
  } else if (material_type == 1) {
    double rand_an = (double)rand()/RAND_MAX*PI;
    struct point2D tang;
    tang.px = -normal.py;
    tang.py = normal.px;
    struct ray2D scat_ray;
    scat_ray.p.px=intersectPt.px;
    scat_ray.p.py=intersectPt.py;
    scat_ray.d.px=tang.px+cos(rand_an);
    scat_ray.d.py=tang.py+sin(rand_an);
    scat_ray.R=ray->R;
    scat_ray.G=ray->G;
    scat_ray.B=ray->B;
    scat_ray.H=ray->H;
    scat_ray.inside_out=ray->inside_out;
    scat_ray.monochromatic=ray->monochromatic;
    propagateRay(&scat_ray, depth+1);
  } else if(material_type == 2) {

  }
}

// Step 5 - Decide how to handle the ray's bounce at the intersection. You will have
//          to provide code for 3 cases:
//          If material type = 0, you have a mirror-reflecting object. 
//                                Create a ray in the mirror reflection direction,
//                                with the same colour as the incoming ray, and
//                                with origin at the intersection point.
//                                Then call propagateRay() recursively to trace it.
//          if material type = 1, you have a scattering surface. 
//                                Choose a random direction within +- 90 degrees 
//                                from the normal at the intersection. Create a
//                                ray in this direction, with the same colour as
//                                the incoming ray, and origin at the intersection,
//                                then call propagateRay() recursively to trace it.
//          if material type = 2, you have a refracting (transparent) material.
// 				   Here you need to process two rays:
//                                * First, determine how much of the incoming light is
//                                  reflected and how much is transmitted, using 
//				     Schlick's approximation:
// 					 R0 = ((n1-n2)/(n1+n2))^2   
// 					 R(theta)=R0+((1-R0)*(1-cos(theta))^5)
//				     If the ray is travelling from air to the inside
//                                  of an object, n1=1, n2=object's index of refraction.
//                                  If the ray is travelling from inside an object
//                                  back onto air, n1=object's index of refraction, n2=1
//				     And 'theta' is the angle between the normal and the
// 				     ray direction.
//				     R(theta) gives the amount Rs of reflected light, 
//				     1.0-R(theta) gives the amount Rt of transmitted light.
//                                * Now, make a ray in the mirror-reflection direction
//				     (same as for material type 0), with the same colour
//				     as the incoming ray, but with intensity modulated
//				     by Rs. (e.g. if the incoming's colour is R,G,B,
//                                  the reflected ray's colour will be R*Rs, G*Rs, B*Rs)
//				     trace this ray.
//				   * Make a ray in the refracted-ray direction. The 
//				     angle for the transmitted ray is given by Snell's law
//				     n1*sin(theta1) = n2*sin(theta2). The colour of the
//				     transmitted ray is the same as the incoming ray but
//			             modulated by Rt. Trace this ray.
 //	That's it! you're done!

   
}

void intersectRay(struct ray2D *ray, struct point2D *p, struct point2D *n, double *lambda, int *type, double *r_idx)
{
 /*
  This function checks for intersection between the ray and any objects in the objects 
  array. The objects are circles, so we are in fact solving for the intersection
  between a ray and a circle.
  
  For a unit circle centered at the origin, we would have the equation
  
  x^2 + y^2 = 1
  
  Using vector notation, with C=[x y]', we get
  
  ||C||^2 = 1
  
  A point on the ray is given by p + lambda*d
  
  Substituting in the equation for the circle we have 
  
  (p + lambda*d)(p + lambda*d) - 1 = 0
  
  If we expand the product above (here the product of two vectors is a DOT product), 
  we can form a quadratic equation
  
  A*lambda^2 + B*lambda + C = 0
  
  Which as you know, has a very simple solution. 
  
  Your task is to 
  * Figure out A, B, and C, considering that your circles don't necessarily have r=1
  * Figure out how to deal with the fact that circles in the scene are likely
    *not* centered at the origin
    
  Then implement the code that will find the value(s) of lambda at the intersection(s).
  
  Note that you can have no intersections, 1 intersection, or 2 intersections
  
  This function *must* find the closest intersection (if any) and update the value
  of lambda, the intersection point p, the normal n at the intersection, 
  the corresponding object's material type (needed outside here to figure out how
  to handle the light's bouncing off this object), and the index of refraction for
  the object (needed if this is a transparent object). 
  
  You must make sure that n is a unit-length vector.
 */
 
 /**********************************************************************************
  * TO DO: Complete this function to find the closest intersection between the
  *        ray and any objects in the scene, as well as the values at the
  *        intersection that will be needed to determine how to bounce/refract the
  *	   ray.
  * *******************************************************************************/
 
 // Step 1: loop over all objects
for (int i = 0; i < MAX_OBJECTS; i++) {
    struct circ2D* c = &objects[i];
    if (c->r ==  -1) {
      continue;
    }
    double A = dot(&ray->d, &ray->d);
    double B =  2 * dot(&ray->d, &ray->p) - 2 * dot(&ray->d, &c->c);
    double C = dot(&c->c, &c->c) - 2 * dot(&ray->p, &c->c) + dot(&ray->p, &ray->p) - pow(c->r,2);
    // solve lambda using quadratic equation
    double discrim = pow(B, 2) - 4 * A * C;
    if (discrim < 0) {
        continue;
    }
    double lambda1 = (-B + sqrt(discrim)) / (2 * A);
    double lambda2 = (-B - sqrt(discrim)) / (2 * A);
    // check if lambda is positive and smaller than current lambda
    if (lambda1 > 0 && lambda1 < *lambda) {
        // update lambda
        *lambda = lambda1;
        // update intersection point
        p->px = ray->p.px + lambda1 * ray->d.px;
        p->py = ray->p.py + lambda1 * ray->d.py;
        // update normal
        n->px = p->px - c->c.px;
        n->py = p->py - c->c.py;
        normalize(n);
        // update material type
        *type = c->material_type;
        // update index of refraction
        *r_idx = c->r_idx;
    }
    if (lambda2 > 0 && lambda2 < *lambda) {
        // update lambda
        *lambda = lambda2;
        // update intersection point
        p->px = ray->p.px + lambda2 * ray->d.px;
        p->py = ray->p.py + lambda2 * ray->d.py;
        // update normal
        n->px = p->px - c->c.px;
        n->py = p->py - c->c.py;
        normalize(n);
        // update material type
        *type = c->material_type;
        // update index of refraction
        *r_idx = c->r_idx;
    }
}
}
