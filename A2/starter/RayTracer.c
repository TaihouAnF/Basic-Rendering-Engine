/*
  CSC D18 - RayTracer code.

  Written Dec. 9 2010 - Jan 20, 2011 by F. J. Estrada
  Freely distributable for adacemic purposes only.

  Uses Tom F. El-Maraghi's code for computing inverse
  matrices. You will need to compile together with
  svdDynamic.c

  You need to understand the code provided in
  this file, the corresponding header file, and the
  utils.c and utils.h files. Do not worry about
  svdDynamic.c, we need it only to compute
  inverse matrices.

  You only need to modify or add code in sections
  clearly marked "TO DO" - remember to check what
  functionality is actually needed for the corresponding
  assignment!

  Last updated: Aug. 2017   - F.J.E.
*/

/*****************************************************************************
* COMPLETE THIS TEXT BOX:
*
* 1) Student Name:		Anson Feng
* 2) Student Name:		Peter Chou
*
* 1) Student number:  1004721955
* 2) Student number:  1004295693
* 
* 1) UtorID:          fengdian
* 2) UtorID:          choulu1
* 
* We hereby certify that the work contained here is our own
*
* _______Anson Feng___             _____Peter Chou_____
* (sign with your name)            (sign with your name)
********************************************************************************/

#include "utils.h"	// <-- This includes RayTracer.h

// A couple of global structures and data: An object list, a light list, and the
// maximum recursion depth
struct object3D *object_list;
struct pointLS *light_list;
struct textureNode *texture_list;
int MAX_DEPTH;

void buildScene(void)
{
#include "buildscene.c"		// <-- Import the scene definition!
}

void rtShade(struct object3D *obj, struct point3D *p, struct point3D *n, struct ray3D *ray, int depth, double a, double b, struct colourRGB *col)
{
    // This function implements the shading model as described in lecture. It takes
    // - A pointer to the first object intersected by the ray (to get the colour properties)
    // - The coordinates of the intersection point (in world coordinates)
    // - The normal at the point
    // - The ray (needed to determine the reflection direction to use for the global component, as well as for
    //   the Phong specular component)
    // - The current racursion depth
    // - The (a,b) texture coordinates (meaningless unless texture is enabled)
    //
    // Returns:
    // - The colour for this ray (using the col pointer)
    //

    struct colourRGB tmp_col;	// Accumulator for colour components
    double R,G,B;			// Colour for the object in R G and B

    // This will hold the colour as we process all the components of
    // the Phong illumination model
    tmp_col.R=0;
    tmp_col.G=0;
    tmp_col.B=0;

    if (obj->texImg==NULL)		// Not textured, use object colour
    {
        R=obj->col.R;
        G=obj->col.G;
        B=obj->col.B;
    }
    else
    {
        // Get object colour from the texture given the texture coordinates (a,b), and the texturing function
        // for the object. Note that we will use textures also for Photon Mapping.
        obj->textureMap(obj->texImg,a,b,&R,&G,&B);
    }

    //////////////////////////////////////////////////////////////
    // TO DO: Implement this function. Refer to the notes for
    // details about the shading model.
    //////////////////////////////////////////////////////////////

    // Local component, setting up a shadow ray
    // and perfect reflection of LS, camera_dir
    struct pointLS *curr_ls = light_list;
    struct ray3D shadow_ray;
    struct point3D *m = newPoint(n->px, n->py, n->pz);
    struct point3D *camera_dir = newPoint(ray->p0.px, ray->p0.py, ray->p0.pz);
    subVectors(p, camera_dir);
    normalize(camera_dir);
    
    // Setting up shadow ray direction
    struct point3D shadow_direction = curr_ls->p0;
    // b = b - a, where a is p and b is shadow_direction
    // we store curr_ls-p0 in shadow_direction first
    subVectors(p, &shadow_direction);
    shadow_direction.pw = 0;

    initRay(&shadow_ray, p, &shadow_direction);

    // Initialize shadow lambda,
    // and temp_varobj for findFirstHit to return
    double shadow_lambda;
    double shadow_a, shadow_b; // might need to initialize in later
    // assignment, same for all those a and b
    struct object3D *temp_var_obj;
    struct point3D shadow_temp_p, shadow_temp_n;
    findFirstHit(&shadow_ray, &shadow_lambda, obj, &temp_var_obj,
                 &shadow_temp_p, &shadow_temp_n, &shadow_a, &shadow_b);

    if (shadow_lambda > 0.0 && shadow_lambda < 1.0) {
        tmp_col.R += obj->alb.ra;
        tmp_col.G += obj->alb.ra;
        tmp_col.B += obj->alb.ra;
        m->px = -1;
        m->py = -1;
        m->pz = -1;
    } else {
        normalize(&shadow_direction);
        double dot_intensity_max = dot(n, &shadow_direction);
        if (obj->frontAndBack && dot_intensity_max < 0.0) {
            m->px = -m->px;
            m->py = -m->py;
            m->pz = -m->pz;
            dot_intensity_max = fabs(dot_intensity_max);
        }
        double coeff = max(0, 2 * dot_intensity_max);
        m->px = coeff * m->px;
        m->py = coeff * m->py;
        m->pz = coeff * m->pz;
        subVectors(&shadow_direction, m);
        normalize(m);
        double specular =  pow(max(0.0, dot(camera_dir, m)), obj->shinyness);
        // Add all component in phong model
        tmp_col.R += obj->alb.ra + (obj->alb.rd * curr_ls->col.R * max(0.0, dot_intensity_max)) + obj->alb.rs * curr_ls->col.R * specular;
        tmp_col.G += obj->alb.ra + (obj->alb.rd * curr_ls->col.G * max(0.0, dot_intensity_max)) + obj->alb.rs * curr_ls->col.G * specular;
        tmp_col.B += obj->alb.ra + (obj->alb.rd * curr_ls->col.B * max(0.0, dot_intensity_max)) + obj->alb.rs * curr_ls->col.B * specular;
    }
    // Global component
    // if (depth < MAX_DEPTH) {
    //   if (obj has specular)
    // }
    //   if (obj has refraction)
    //      if (not fully refraction)
    // Ig = Ise + Ire
    // else Ig = 0
    // tmp_col += Ig;
    /* Currently use this to generate signature. normal */
    // Be sure to update 'col' with the final colour computed here!
    tmp_col.R = R * tmp_col.R;
    tmp_col.G = G * tmp_col.G;
    tmp_col.B = B * tmp_col.B;
    col->R = (tmp_col.R > 1.0) ? 1.0 : tmp_col.R;
    col->G = (tmp_col.G > 1.0) ? 1.0 : tmp_col.G;
    col->B = (tmp_col.B > 1.0) ? 1.0 : tmp_col.B;
    free(m);
    free(camera_dir);
    return;
}

void findFirstHit(struct ray3D *ray, double *lambda, struct object3D *Os, struct object3D **obj, struct point3D *p, struct point3D *n, double *a, double *b)
{
    // Find the closest intersection between the ray and any objects in the scene.
    // Inputs:
    //   *ray    -  A pointer to the ray being traced
    //   *Os     -  'Object source' is a pointer toward the object from which the ray originates. It is used for reflected or refracted rays
    //              so that you can check for and ignore self-intersections as needed. It is NULL for rays originating at the center of
    //              projection
    // Outputs:
    //   *lambda -  A pointer toward a double variable 'lambda' used to return the lambda at the intersection point
    //   **obj   -  A pointer toward an (object3D *) variable so you can return a pointer to the object that has the closest intersection with
    //              this ray (this is required so you can do the shading)
    //   *p      -  A pointer to a 3D point structure so you can store the coordinates of the intersection point
    //   *n      -  A pointer to a 3D point structure so you can return the normal at the intersection point
    //   *a, *b  -  Pointers toward double variables so you can return the texture coordinates a,b at the intersection point

    /////////////////////////////////////////////////////////////
    // TO DO: Implement this function. See the notes for
    // reference of what to do in here. Done
    /////////////////////////////////////////////////////////////

    *lambda = -1.0;
    struct object3D *obj_head = object_list;
    while (obj_head != NULL) {
        // Check if current obj is the object source
        if (obj_head != Os) {

            // Find intersection with current obj
            double temp_lambda = -1.0;     // temp lambda to store the result
            struct point3D temp_p;  // temp intersect point, only used when temp_lambda is valid
            struct point3D temp_n;  // temp normal, same as above
            double temp_a, temp_b;  // temp a and b, same as above
            obj_head->intersect(obj_head, ray, &temp_lambda, &temp_p, &temp_n, &temp_a, &temp_b);

            // If current smallest lambda is small than 0
            // OR new lambda is smaller than current smallest.
            // And the new lambda should be > 0
            if ((*lambda < 0.0 || temp_lambda < *lambda) && temp_lambda > 0.0) {

                // Update lambda
                *lambda = temp_lambda;

                // Update intersection point
                p->px = temp_p.px;
                p->py = temp_p.py;
                p->pz = temp_p.pz;
                p->pw = 1;

                // Update normal
                n->px = temp_n.px;
                n->py = temp_n.py;
                n->pz = temp_n.pz;
                n->pw = 1;

                // Update texture a and b
                *a = temp_a;
                *b = temp_b;

                // Update the obj
                *obj = obj_head;
            }
        }
        obj_head = obj_head->next;
    }

}

void rayTrace(struct ray3D *ray, int depth, struct colourRGB *col, struct object3D *Os)
{
    // Trace one ray through the scene.
    //
    // Parameters:
    //   *ray   -  A pointer to the ray being traced
    //   depth  -  Current recursion depth for recursive raytracing
    //   *col   - Pointer to an RGB colour structure so you can return the object colour
    //            at the intersection point of this ray with the closest scene object.
    //   *Os    - 'Object source' is a pointer to the object from which the ray
    //            originates so you can discard self-intersections due to numerical
    //            errors. NULL for rays originating from the center of projection.

    double lambda;		// Lambda at intersection
    double a,b;		// Texture coordinates
    struct object3D *obj;	// Pointer to object at intersection
    struct point3D p;	// Intersection point
    struct point3D n;	// Normal at intersection
    struct colourRGB I;	// Colour returned by shading function

    if (depth>MAX_DEPTH)	// Max recursion depth reached. Return invalid colour.
    {
        col->R=-1;
        col->G=-1;
        col->B=-1;
        return;
    }

    ///////////////////////////////////////////////////////
    // TO DO: Complete this function. Refer to the notes
    // if you are unsure what to do here. Done
    ///////////////////////////////////////////////////////
    findFirstHit(ray, &lambda, Os, &obj, &p, &n, &a, &b);

    if (lambda > 0) {

        // Shading to obtain the color
        rtShade(obj, &p, &n, ray, depth, a, b, &I);

        // Update the color
        col->R = I.R;
        col->G = I.G;
        col->B = I.B;
    }
        // if lambda <=0, col = background/originate obj,
        // as we set before passing in rayTrace()
    else {
        col->R = col->R;
        col->G = col->G;
        col->B = col->B;
    }

}

int main(int argc, char *argv[])
{
    // Main function for the raytracer. Parses input parameters,
    // sets up the initial blank image, and calls the functions
    // that set up the scene and do the raytracing.
    struct image *im;	// Will hold the raytraced image
    struct view *cam;	// Camera and view for this scene
    int sx;		// Size of the raytraced image
    int antialiasing;	// Flag to determine whether antialiaing is enabled or disabled
    char output_name[1024];	// Name of the output file for the raytraced .ppm image
    struct point3D e;		// Camera view parameters 'e', 'g', and 'up'
    struct point3D g;
    struct point3D up;
    double du, dv;			// Increase along u and v directions for pixel coordinates
    struct point3D pc,d;		// Point structures to keep the coordinates of a pixel and
    // the direction or a ray
    struct ray3D ray;		// Structure to keep the ray from e to a pixel
    struct colourRGB col;		// Return colour for raytraced pixels
    struct colourRGB background;   // Background colour
    int i,j;			// Counters for pixel coordinates
    unsigned char *rgbIm;

    if (argc<5)
    {
        fprintf(stderr,"RayTracer: Can not parse input parameters\n");
        fprintf(stderr,"USAGE: RayTracer size rec_depth antialias output_name\n");
        fprintf(stderr,"   size = Image size (both along x and y)\n");
        fprintf(stderr,"   rec_depth = Recursion depth\n");
        fprintf(stderr,"   antialias = A single digit, 0 disables antialiasing. Anything else enables antialiasing\n");
        fprintf(stderr,"   output_name = Name of the output file, e.g. MyRender.ppm\n");
        exit(0);
    }
    sx=atoi(argv[1]);
    MAX_DEPTH=atoi(argv[2]);
    if (atoi(argv[3])==0) antialiasing=0; else antialiasing=1;
    strcpy(&output_name[0],argv[4]);

    fprintf(stderr,"Rendering image at %d x %d\n",sx,sx);
    fprintf(stderr,"Recursion depth = %d\n",MAX_DEPTH);
    if (!antialiasing) fprintf(stderr,"Antialising is off\n");
    else fprintf(stderr,"Antialising is on\n");
    fprintf(stderr,"Output file name: %s\n",output_name);

    object_list=NULL;
    light_list=NULL;
    texture_list=NULL;

    // Allocate memory for the new image
    im=newImage(sx, sx);
    if (!im)
    {
        fprintf(stderr,"Unable to allocate memory for raytraced image\n");
        exit(0);
    }
    else rgbIm=(unsigned char *)im->rgbdata;

    ///////////////////////////////////////////////////
    // TO DO: You will need to implement several of the
    //        functions below. For Assignment 2, you can use
    //        the simple scene already provided. But
    //        for Assignment 3 you need to create your own
    //        *interesting* scene.
    ///////////////////////////////////////////////////
    buildScene();		// Create a scene. This defines all the
    // objects in the world of the raytracer

    //////////////////////////////////////////
    // TO DO: For Assignment 2 you can use the setup
    //        already provided here. For Assignment 3
    //        you may want to move the camera
    //        and change the view parameters
    //        to suit your scene.
    //////////////////////////////////////////

    // Mind the homogeneous coordinate w of all vectors below. DO NOT
    // forget to set it to 1, or you'll get junk out of the
    // geometric transformations later on.

    // Camera center is at (0,0,-1)
    e.px=0;
    e.py=0;
    e.pz=-1;
    e.pw=1;

    // To define the gaze vector, we choose a point 'pc' in the scene that
    // the camera is looking at, and do the vector subtraction pc-e.
    // Here we set up the camera to be looking at the origin.
    g.px=0-e.px;
    g.py=0-e.py;
    g.pz=0-e.pz;
    g.pw=1;
    // In this case, the camera is looking along the world Z axis, so
    // vector w should end up being [0, 0, -1]

    // Define the 'up' vector to be the Y axis
    up.px=0;
    up.py=1;
    up.pz=0;
    up.pw=1;

    // Set up view with given the above vectors, a 4x4 window,
    // and a focal length of -1 (why? where is the image plane?)
    // Note that the top-left corner of the window is at (-2, 2)
    // in camera coordinates.
    cam=setupView(&e, &g, &up, -1, -2, 2, 4);

    if (cam==NULL)
    {
        fprintf(stderr,"Unable to set up the view and camera parameters. Our of memory!\n");
        cleanup(object_list,light_list, texture_list);
        deleteImage(im);
        exit(0);
    }

    // Set up background colour here
    background.R=0;
    background.G=0;
    background.B=0;

    // Do the raytracing
    //////////////////////////////////////////////////////
    // TO DO: You will need code here to do the raytracing
    //        for each pixel in the image. Refer to the
    //        lecture notes, in particular, to the
    //        raytracing pseudocode, for details on what
    //        to do here. Make sure you undersand the
    //        overall procedure of raytracing for a single
    //        pixel. Done
    //////////////////////////////////////////////////////
    du=cam->wsize/(sx-1);		// du and dv. In the notes in terms of wl and wr, wt and wb,
    dv=-cam->wsize/(sx-1);		// here we use wl, wt, and wsize. du=dv since the image is
    // and dv is negative since y increases downward in pixel
    // coordinates and upward in camera coordinates.

    fprintf(stderr,"View parameters:\n");
    fprintf(stderr,"Left=%f, Top=%f, Width=%f, f=%f\n",cam->wl,cam->wt,cam->wsize,cam->f);
    fprintf(stderr,"Camera to world conversion matrix (make sure it makes sense!):\n");
    printmatrix(cam->C2W);
    fprintf(stderr,"World to camera conversion matrix:\n");
    printmatrix(cam->W2C);
    fprintf(stderr,"\n");

    fprintf(stderr,"Rendering row: ");
    for (j=0;j<sx;j++)		// For each of the pixels in the image
    {
        fprintf(stderr,"%d/%d, ",j,sx);
        for (i=0;i<sx;i++)
        {
            ///////////////////////////////////////////////////////////////////
            // TO DO - complete the code that should be in this loop to do the
            //         raytracing! Done
            ///////////////////////////////////////////////////////////////////

            // Setting up the point on the image plane in camera coordinate
            pc.px = cam->wl + (i * du); // u, which is x in pt on image plane
            pc.py = cam->wt + (j * dv); // v, which is y in pt on image plane
            pc.pz = -1;                 // w, focal length, the z distance from eye point to plane
            pc.pw = 1;                  // homogeneous point so the w is 1

            // Convert the point to the World coordinate and pc.pw is still 1,
            // check C2W and matVecMult to convince
            matVecMult(cam->C2W, &pc);

            // Getting direction vector in world coordinate by pc(World) - e(World)
            d.px = pc.px - e.px;        // x of direction vector in World coordinate
            d.py = pc.py - e.py;        // y of direction vector in World coordinate
            d.pz = pc.pz - e.pz;        // z of direction vector in World coordinate
            d.pw = 0;                   // Homogeneous, vector as 0, otherwise it's incorrect
            // and pc.pw - e.pw supposed to be 0
            // normalize it
            normalize(&d);

            // Create the ray
            initRay(&ray, &pc, &d);

            // Initialize the color, we might use it to handle
            // rayTrace has a negative lambda and we need to set
            // the color to background.
            col.R = background.R;
            col.G = background.G;
            col.B = background.B;

            // Trace the ray to get the color, if lambda > 0, color will be changed
            rayTrace(&ray, 1, &col, NULL);  // NULL as it is the beginning of the tracing

            /*
            *   Set pixel color, TODO: needs double check:
            *   We can see rgbIm as a flatten version of the image plane,
            *   we have sx columns and sy rows, as i represents column and j rows,
            *   we put all the indices into one row =>
            *   [0 1 2 ... i ... sx sx+1 ... 2sx ... j*sx+i ... sysx]
            *   we can see sx+1 will be the first index of the second row, so
            *   ij-th index would be j*sx+i in the array. And further, we expand
            *   each index to have 3 index to store R, G and B.
            */

            *(rgbIm + (j * sx + i) * 3) = (unsigned char)(col.R * 255);
            *(rgbIm + (j * sx + i) * 3 + 1) = (unsigned char)(col.G * 255);
            *(rgbIm + (j * sx + i) * 3 + 2) = (unsigned char)(col.B * 255);

        } // end for i
    } // end for j

    fprintf(stderr,"\nDone!\n");

    // Output rendered image
    imageOutput(im,output_name);

    // Exit section. Clean up and return.
    cleanup(object_list,light_list,texture_list);		// Object, light, and texture lists
    deleteImage(im);					// Rendered image
    free(cam);						// camera view
    exit(0);
}
