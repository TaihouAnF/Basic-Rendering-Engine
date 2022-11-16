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
// how many points to randomly sample from an area light source
int K = 50;

void buildScene(void)
{
#include "buildscene.c"		// <-- Import the scene definition!
}

void reflectionDirection(struct point3D *d, struct point3D *n)
{
    normalize(n);;
    //reflect d about n
    double dn = dot(d, n);
    d->px = -(d->px - 2 * dn * n->px);
    d->py = -(d->py - 2 * dn * n->py);
    d->pz = -(d->pz - 2 * dn * n->pz);
    d->pw = 1;
    normalize(d);
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

    if (obj->normalMapped)
    {
        // handle if the object is normal mapped
        image *normalmap = obj->normalMap;
        int x = (int) (a * (normalmap->sx - 1));
        int y = (int) (b * (normalmap->sy - 1));
        double* normaldata = (double *) normalmap->rgbdata;
        n->px = (double) *(normaldata + ((y * normalmap->sx) + x) * 3);
        n->py = (double) *(normaldata + ((y * normalmap->sx) + x) * 3 + 1);
        n->pz = (double) *(normaldata + ((y * normalmap->sx) + x) * 3 + 2);
    }

    //////////////////////////////////////////////////////////////
    // TO DO: Implement this function. Refer to the notes for
    // details about the shading model.
    //////////////////////////////////////////////////////////////

    // Local component, setting up a shadow ray
    // and perfect reflection of LS, camera_dir
    struct pointLS *curr_ls = light_list;
    struct ray3D shadow_ray;
    struct point3D *camera_dir = newPoint(-ray->d.px, -ray->d.py, -ray->d.pz);
    // direction pointing to reflection of camera ray
    struct point3D *r_camera = newPoint(camera_dir->px, camera_dir->py, camera_dir->pz);
    reflectionDirection(r_camera, n);
    // Initialize shadow lambda,
    // and temp_varobj for findFirstHit to return
    double shadow_lambda;
    double shadow_a, shadow_b; // might need to initialize in later
    // assignment, same for all those a and b
    struct object3D *temp_var_obj;
    struct point3D shadow_temp_p, shadow_temp_n;
    // Create temp color for global component
    struct colourRGB reflection_col, refraction_col;
    reflection_col.R = 0.0;
    reflection_col.G = 0.0;
    reflection_col.B = 0.0;

    refraction_col.R = 0.0;
    refraction_col.G = 0.0;
    refraction_col.B = 0.0;

    // objects will have the ambient component no matter what
    tmp_col.R += R * obj->alb.ra;
    tmp_col.G += G * obj->alb.ra;
    tmp_col.B += B * obj->alb.ra;
    struct object3D *obj_head = object_list;
    double x;
    double y;
    double z;
    while (obj_head != NULL) {
        if (obj_head != obj && obj_head->isLightSource) {
            // Sample K points from the area light source
            for (int i = 0; i < K; i++) {
                obj_head->randomPoint(obj_head, &x, &y, &z);
                struct point3D* shadow_direction = newPoint(x, y, z);
                //printf("x: %f, y: %f, z: %f\n", x, y, z);
                subVectors(p, shadow_direction);
                // Setting up shadow ray
                double dot_intensity_max = dot(n, shadow_direction);
                struct point3D* normal = newPoint(n->px, n->py, n->pz);
                if (obj->frontAndBack && dot_intensity_max < 0.0) {
                    normal->px = -normal->px;
                    normal->py = -normal->py;
                    normal->pz = -normal->pz;
                }
                initRay(&shadow_ray, p, shadow_direction);
                struct point3D *r_light = newPoint(shadow_direction->px, shadow_direction->py, shadow_direction->pz);
                reflectionDirection(r_light, normal);
                // // direction pointing to light source
                struct point3D *s = newPoint(shadow_direction->px, shadow_direction->py, shadow_direction->pz);
                normalize(s);
                // Check if the shadow ray intersects with any object
                findFirstHit(&shadow_ray, &shadow_lambda, obj, &temp_var_obj,
                             &shadow_temp_p, &shadow_temp_n, &shadow_a, &shadow_b);
                // NOTE: for some reason setting shadow_lambda < 1.0 doesn't work
                if (shadow_lambda > 0.0 && shadow_lambda < 0.999999) {}
                else {
                    // Flip the normal if object is supposed to be illuminated front and back
                    // direction pointing to reflection of light source ray
                    double specular =  pow(max(0.0, dot(camera_dir, r_light)), obj->shinyness);
                    double ns = dot(normal, s);
                    // Add the phong model and finish local
                    double intensity =  ((double) 1) / ((double)K);
                    tmp_col.R += R * intensity * (obj->alb.rd * obj_head->col.R * max(0.0, ns))
                                  + intensity * obj->alb.rs * obj_head->col.R * specular;
                    tmp_col.G += G * intensity * (obj->alb.rd * obj_head->col.G * max(0.0, ns))
                                  + intensity * obj->alb.rs * obj_head->col.G * specular;
                    tmp_col.B += B * intensity * (obj->alb.rd * obj_head->col.B * max(0.0, ns))
                                  + intensity * obj->alb.rs * obj_head->col.B * specular;
                }
                free(r_light);
                free(s);
                free(normal);
                free(shadow_direction);
            }
        }
        obj_head = obj_head->next;
    }

    // Global
    struct ray3D reflection_ray, refraction_ray;
    if (depth < MAX_DEPTH) {
        if (obj->alb.rs > 0) {
            initRay(&reflection_ray, p, r_camera);
            reflection_ray.inside = ray->inside;
            reflection_ray.ref_ind_stack = stackCopy(ray->ref_ind_stack);
            reflection_ray.first_ray = ray->first_ray;
            rayTrace(&reflection_ray, depth + 1, &reflection_col, obj);
            // I = Il + Ig, Il could be only ambient or full phong model
            tmp_col.R += obj->alb.rg * reflection_col.R;
            tmp_col.G += obj->alb.rg * reflection_col.G;
            tmp_col.B += obj->alb.rg * reflection_col.B;
        }

        // refraction
        if (obj->alpha < 1) {
            // Calculate the refraction direction
            // dt = rb + (rc - sqrt(1 - r^2*(1 - c^2)))n
            // First, we check the normal and incident ray direction
            struct point3D *refraction_direction = newPoint(n->px, n->py, n->pz);
            double c = dot(&(ray->d), n);
            double r;
            int inside = ray->inside;
            int going_out = 0;
            int going_in = 0;
            if (c <= 0) {
                // The ray is *Outside* of the surface, it is about to enter
                // so normal is the same, and we need to make c to be positive
                c = -c;
                r = ray->ref_ind_stack->current_index / obj->r_index;
                going_in = 1;
                inside = 1;
            } else if(c > 0 && inside) {
                // The ray is *Inside* of the surface since c is positive, and we 
                // need to flip the normal(refraction_direction) in this case
                if (!ray->ref_ind_stack || !ray->ref_ind_stack->next) return;
                double leaving = ray->ref_ind_stack->current_index;
                r = leaving / ray->ref_ind_stack->next->current_index;
                refraction_direction->px = -refraction_direction->px;
                refraction_direction->py = -refraction_direction->py;
                refraction_direction->pz = -refraction_direction->pz;
                going_out = 1;
            } else {
                stackFree(ray->ref_ind_stack);
                return;
            }
            //  1 - sqrt(1 - r^2*(1 - c^2))
            double intercheck = 1 - ((r * r) * (1 - (c * c)));
            if (intercheck >= 0) {
                // (rc - sqrt(1 - r^2*(1 - c^2)))n in our case refraction_direction 
                // is the n,  we need to use addVector() 
                double coeff = r * c - sqrt(intercheck);
                refraction_direction->px *= coeff;
                refraction_direction->py *= coeff;
                refraction_direction->pz *= coeff;
                // rb
                struct point3D *temp_incident_ray_d = newPoint(r * ray->d.px, r * ray->d.py, r * ray->d.pz);
                // dt = rb + (rc - sqrt(1 - r^2*(1 - c^2)))n b=a+b
                addVectors(temp_incident_ray_d, refraction_direction);
                normalize(refraction_direction);
                initRay(&refraction_ray, p, refraction_direction);
                if (going_out) {
                    refraction_ray.inside = 0;
                    refraction_ray.ref_ind_stack = stackPop(ray->ref_ind_stack);
                } else if (going_in) {
                    refraction_ray.inside = inside;
                    struct refraction_ind_stk *new_stack_top = newStackInstance(obj->r_index);
                    refraction_ray.ref_ind_stack = stackInsert(new_stack_top, ray->ref_ind_stack);
                } else {
                    stackFree(ray->ref_ind_stack);
                    return;
                }
                refraction_ray.first_ray = 0;
                rayTrace(&refraction_ray, depth + 1, &refraction_col, obj);
                tmp_col.R *= obj->alpha;
                tmp_col.G *= obj->alpha;
                tmp_col.B *= obj->alpha;
                
                tmp_col.R += ((1 - obj->alpha) * refraction_col.R);
                tmp_col.G += ((1 - obj->alpha) * refraction_col.G);
                tmp_col.B += ((1 - obj->alpha) * refraction_col.B);
                free(temp_incident_ray_d);
            }
            free(refraction_direction);
        }
    }


    // Update the final color of the pixel
    col->R = (tmp_col.R > 1.0) ? 1.0 : tmp_col.R;
    col->G = (tmp_col.G > 1.0) ? 1.0 : tmp_col.G;
    col->B = (tmp_col.B > 1.0) ? 1.0 : tmp_col.B;
    free(r_camera);
    free(camera_dir);
    stackFree(ray->ref_ind_stack);
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
        if (obj_head != Os || ray->inside) {

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

        // If the obj the ray intersected is the areaLS
        // it should just return the color of LS. Needs checking
        if (obj->isLightSource) {
            I.R = obj->col.R;
            I.G = obj->col.G;
            I.B = obj->col.B;
        }
        else { // Shading to obtain the color
            if (ray->first_ray) {
                ray->ref_ind_stack = newStackInstance(1.0); // Assuming initialized ray starts from vacuum
            }
            rtShade(obj, &p, &n, ray, depth, a, b, &I);
        }

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
    struct colourRGB sum_color;  // Summary of color when antialiasing
    struct colourRGB temp_color; // Container to store temp RGB when antialising
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

#pragma omp parallel for schedule(dynamic,16) private(ray, col, i, j, sum_color, temp_color)
    for (j=0;j<sx;j++)		// For each of the pixels in the image
    {
        fprintf(stderr,"%d/%d, ",j,sx);
        for (i=0;i<sx;i++)
        {
            ///////////////////////////////////////////////////////////////////
            // TO DO - complete the code that should be in this loop to do the
            //         raytracing! Done
            ///////////////////////////////////////////////////////////////////
            if (antialiasing) {
                // Initialize the sum color for summing all four location
                // and a temp color container TODOs: need checking
                sum_color.R = 0;
                sum_color.G = 0;
                sum_color.B = 0;

                for (double sub_i = .25; sub_i <= .75; sub_i += .5) {
                    for (double sub_j = .25; sub_j <= .75; sub_j += .5) {
                        // cam->wl+(i*du) and cam->wt+(j*dv)
                        // is the x and y of the center of pixel; while
                        // ((sub_i-.5)*du) and ((sub_j-.5)*dv) are offset
                        pc.px = cam->wl + (i * du) + ((sub_i - 0.5) * du);
                        pc.py = cam->wt + (j * dv) + ((sub_j - 0.5) * dv);
                        pc.pz = -1;
                        pc.pw = 1;
                        matVecMult(cam->C2W, &pc);

                        d = pc;
                        subVectors(&e, &d);
                        d.pw = 0;
                        normalize(&d);

                        initRay(&ray, &pc, &d);
                        temp_color.R = background.R;
                        temp_color.G = background.G;
                        temp_color.B = background.B;
                        ray.first_ray = 1;
                        rayTrace(&ray, 1, &temp_color, NULL);

                        sum_color.R += temp_color.R;
                        sum_color.G += temp_color.G;
                        sum_color.B += temp_color.B;
                    }
                }
            }

            // Setting up the point on the image plane in camera coordinate
            pc.px = cam->wl + (i * du); // u, which is x in pt on image plane
            pc.py = cam->wt + (j * dv); // v, which is y in pt on image plane
            pc.pz = -1;                 // w, focal length, the z distance from eye point to plane
            pc.pw = 1;                  // homogeneous point so the w is 1

            // Convert the point to the World coordinate and pc.pw is still 1,
            // check C2W and matVecMult to convince
            matVecMult(cam->C2W, &pc);

            // Getting direction vector in world coordinate by pc(World) - e(World)
            d = pc;
            subVectors(&e, &d);
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

            ray.first_ray = 1;

            // Trace the ray to get the color, if lambda > 0, color will be changed
            rayTrace(&ray, 1, &col, NULL);  // NULL as it is the beginning of the tracing

            if (antialiasing) {
                col.R = (col.R + sum_color.R) / 5;
                col.G = (col.G + sum_color.G) / 5;
                col.B = (col.B + sum_color.B) / 5;
            }

            /*
            *   Set pixel color:
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