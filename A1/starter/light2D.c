/*
  CSC D18 - Assignment 1 - 2D light propagation.
  Starter code: F.J. Estrada, Aug. 2017

  This assignment is intended to help you practice with
  vector and point operations, as well as with simple
  transformations in 2D. 

  You will do this in the context of light propagation
  inside a 2D box that contains objects (circles) with
  different properties - some reflect light like a mirror,
  others scatter light, and yet others behave like
  transparent media.

  Your goal is to write the fundamental components of 
  the light propagation computation, namely:

  - Determining the starting point and direction of rays
  - Determining whether a ray intersects an object or
    hits the wall of the box
  - Whenever a ray hits an object or a wall, bouncing
    or refracting the ray appropriately.

  To complete this assignment you need a solid
  understanding of:

  - Points, vectors, and rays (parametric equation)
  - 2D curves (circles) and their mathematical
    representation
  - 2D transformations (particularly rotation and
    translation)
  
  You will need also to have a solid handle on basic
  C programming with recursion. If you should feel the
  need to brush up on your C programming, do spend some
  time studying on your favorite on-line C programming
  resource.

  ****************************************************
  * You must not modify any part of the code outside
  * sections marked TO DO
  ****************************************************
  
  Note: This project is based on the Tantalum 2D light
        transport project by Benedikt Betterli. 
        Really cool stuff. Have a look here:
        
        https://benedikt-bitterli.me/tantalum/tantalum.html
        (needs a modern browser and a WebGL capable device)
*/

#include "light2D.h"	// READ THIS HEADER FILE CAREFULLY!

// Import scene definition - change buildScene.c to move
// things around inside the box.
#include "buildScene.c"

// The lines below IMPORT YOUR CODE, since this is a 
// very simple program, we won't bother with Makefiles,
// and we'll just have the pre-compiler paste your code
// right here.
#include "rays2D.h"
#include "rays2D.c"

// Let's get some work done!

void hue2RGB(double H, double *R, double *G, double *B)
{
 /* Given a HUE value in [0 1], with 0 for deep red and
  * 1 for purple, obtains the corresponding RGB values
  * that give the equivalent colour
  */

 double C,X;
 
 C=1.0;
 X=C*(1.0-fabs(fmod(6.0*H,2.0)-1.0));

 if (H<1.0/6.0)
 {
   *R=1.0;
   *G=X;
   *B=0;
 }
 else if (H<2.0/6.0)
 {
   *R=X;
   *G=C;
   *B=0;
 }
 else if (H<3.0/6.0)
 {
   *R=0;
   *G=C;
   *B=X;
 }
 else if (H<4.0/6.0)
 {
  *R=0;
  *G=X;
  *B=C;
 }
 else if (H<5.0/6.0)
 {
   *R=X;
   *G=0;
   *B=C;
 }
 else
 {
   *R=C;
   *G=0;
   *B=X;
 }
  
}

void renderObjects(void)
{
 /*
  * Useful for debugging - will overlay the objects in the scene onto the
  * final image in green - turn on/off with the debug flag in rays2D.c
  */
 
 double x,y;
 int xx,yy;
 
 for (int i=0;i<MAX_OBJECTS;i++)
 {
  if (objects[i].r<=0) break;
  for (double ang=0; ang<2*PI; ang+=.001)
  {
   x=objects[i].c.px+(cos(ang)*objects[i].r);
   y=objects[i].c.py+(sin(ang)*objects[i].r);
   x-=W_LEFT;
   y-=W_TOP;
   x=x/(W_RIGHT-W_LEFT);
   y=y/(W_BOTTOM-W_TOP);
   x=x*(sx-1);
   y=y*(sy-1);
   xx=(int)round(x);
   yy=(int)round(y);
   if (xx>=0&&xx<sx&&yy>=0&&yy<sy)
   {
    *(im+((xx+(yy*sx))*3)+0)=0;
    *(im+((xx+(yy*sx))*3)+1)=255;
    *(im+((xx+(yy*sx))*3)+2)=0;   
   }
  }
 }
}

void renderRay(struct point2D *p1, struct point2D *p2, double R, double G, double B)
{
 /*
   This function renders a ray onto the image from point p1 to point p2, with the
   specified colour R,G,B
 */

 double x1,y1,x2,y2,xt,yt;
 int xx,yy;
 double dx,dy;
 double inc;
 
 if (p1->px<W_LEFT-TOL||p1->px>W_RIGHT+TOL||p1->py<W_TOP-TOL||p1->py>W_BOTTOM+TOL||\
     p2->px<W_LEFT-TOL||p2->px>W_RIGHT+TOL||p2->py<W_TOP-TOL||p2->py>W_BOTTOM+TOL)
 {
  fprintf(stderr,"renderRay() - at least one endpoint is outside the image bounds, somewhere there's an error...\n");
  fprintf(stderr,"p1=(%f,%f)\n",p1->px,p1->py);
  fprintf(stderr,"p2=(%f,%f)\n",p2->px,p2->py);
 }
 
 x1=p1->px-W_LEFT;
 y1=p1->py-W_TOP;
 x2=p2->px-W_LEFT;
 y2=p2->py-W_TOP;
 
 x1=x1/(W_RIGHT-W_LEFT);
 y1=y1/(W_BOTTOM-W_TOP);
 x2=x2/(W_RIGHT-W_LEFT);
 y2=y2/(W_BOTTOM-W_TOP);

 x1=x1*(sx-1);
 y1=y1*(sy-1);
 x2=x2*(sx-1);
 y2=y2*(sy-1);

 dx=x2-x1;
 dy=y2-y1;
 
 if (abs(dx)>=abs(dy))
 {
  if (x2<x1)
  {
   xt=x1;
   yt=y1;
   x1=x2;
   y1=y2;
   x2=xt;
   y2=yt;
  }
  
   yt=y1;
   inc=(y2-y1)/abs(x2-x1);
   for (double xt=x1; xt<=x2; xt+=1)
   {
    xx=(int)round(xt);
    yy=(int)round(yt);
    if (xx>=0&&xx<sx&&yy>=0&&yy<sy)
    {
      (*(imRGB+((xx+(yy*sx))*3)+0))+=R;
      (*(imRGB+((xx+(yy*sx))*3)+1))+=G;
      (*(imRGB+((xx+(yy*sx))*3)+2))+=B;
    }
    yt+=inc;
   }
      
 }
 else
 {
  if (y2<y1)
  {
   xt=x1;
   yt=y1;
   x1=x2;
   y1=y2;
   x2=xt;
   y2=yt;
  }

   xt=x1;
   inc=(x2-x1)/abs(y2-y1);
   for (double yt=y1; yt<=y2; yt+=1)
   {
    xx=(int)round(xt);
    yy=(int)round(yt);
    if (xx>=0&&xx<sx&&yy>=0&&yy<sy)
    {
      (*(imRGB+((xx+(yy*sx))*3)+0))+=R;
      (*(imRGB+((xx+(yy*sx))*3)+1))+=G;
      (*(imRGB+((xx+(yy*sx))*3)+2))+=B;
    }
    xt+=inc;
   }
   
 }
 
}

void setPixel(double x, double y, double R, double G, double B)
{
 /* 
  * This function updates the image so that the location corresponding to (x,y) has
  * the specified colour. It handles conversion of coordinates from scene coordinates
  * to image pixel coordinates, and performs a bit of antialiasing
  */
 
 int xx,yy;
 int ii,jj;
 double W;

 if (R<0||G<0||B<0||R>1||G>1||B>1)
   fprintf(stderr,"Invalid RGB colours passed to setPixel() - image will have artifacts!\n");
 
 // Convert to image coordinates
 x-=W_LEFT;
 y-=W_TOP;
 x=x/(W_RIGHT-W_LEFT);
 y=y/(W_BOTTOM-W_TOP);
 x=x*(sx-1);
 y=y*(sy-1);
  
 xx=(int)round(x);
 yy=(int)round(y);
 
 for (int i=xx-1;i<=xx+1;i++)
   for (int j=yy-1;j<=yy+1;j++)
   {
    W=exp(-(((x-i)*(x-i))+((y-j)*(y-j)))*.5);
    if (i>=0&&j>=0&&i<sx&&j<sy)
    {
     (*(imRGB+((i+(j*sx))*3)+0))+=W*R;
     (*(imRGB+((i+(j*sx))*3)+1))+=W*G;
     (*(imRGB+((i+(j*sx))*3)+2))+=W*B;
    }
   }
}

double dot(struct point2D *p, struct point2D *q)
{
 return((p->px*q->px)+(p->py*q->py)); 
}

void normalize(struct point2D *d)
{
 double l;
 l=d->px*d->px;
 l+=(d->py*d->py);
 if (l>0)
 {
  l=sqrt(l);
  d->px=d->px/l;
  d->py=d->py/l;
 }
}

void addCirc(struct point2D *c, double r, int type, double r_idx)
{
 // This adds an object to the object array. Parameters specify
 // the circle's center c, the radius r, type of material, and
 // index of refraction (only meaningful for transparent objects)
 static int num_obj=0;
 
 if (num_obj>=MAX_OBJECTS)
 {
  fprintf(stderr,"List of objects is full!\n");
  return;
 }
 objects[num_obj].c=*c;
 objects[num_obj].r=r;
 objects[num_obj].material_type=type;
 objects[num_obj].r_idx=r_idx;
 num_obj++;
}

int main(int argc, char *argv[])
{
 struct ray2D ray;
 struct point2D p,d;
 double mx,mi,rng;
 FILE *f;

 // Parse command line arguments and validate their range 
 if (argc<5)
 {
  fprintf(stderr,"USAGE: light2D  sx   sy   num_samples   max_depth\n");
  fprintf(stderr,"  sx, sy - image resolution in pixels (in [256 4096])\n");
  fprintf(stderr,"  num_samples - Number of light rays to propagate  (in [1 10,000,000])\n");
  fprintf(stderr,"  max_depth - Maximum recursion depth (in [1 25])\n");
  exit(0);
 }
 sx=atoi(argv[1]);
 sy=atoi(argv[2]);
 num_rays=atoi(argv[3]);
 max_depth=atoi(argv[4]);
 if (sx<256||sy<256||sx>4096||sy>4096||num_rays<1||num_rays>10000000||max_depth<1||max_depth>25)
 {
  fprintf(stderr,"USAGE: light2D  sx   sy   num_samples   max_depth\n");
  fprintf(stderr,"  sx, sy - image resolution in pixels (in [256 4096])\n");
  fprintf(stderr,"  num_samples - Number of light rays to propagate  (in [1 10,000,000])\n");
  fprintf(stderr,"  max_depth - Maximum recursion depth (in [1 25])\n");
  exit(0);
 }
 fprintf(stderr,"Working with:\n");
 fprintf(stderr,"Image size (%d, %d)\n",sx,sy);
 fprintf(stderr,"Number of samples: %d\n",num_rays);
 fprintf(stderr,"Max. recursion depth: %d\n",max_depth);

 // Initialize a blank image for our render
 imRGB=(double *)calloc(sx*sy*3,sizeof(double));
 if (imRGB==NULL)
 {
  fprintf(stderr,"Out of memory?! is this a Commodore 64 you're running on???\n");
  exit(0);
 }

 // Reset walls and objects arrays
 memset(&objects[0],0,MAX_OBJECTS*sizeof(struct circ2D));
 memset(&walls[0],0,4*sizeof(struct ray2D));
 for (int i=0;i<MAX_OBJECTS;i++) objects[i].r=-1; // Make sure we can tell any objects 
						  // not added by buildScene() have
						  // negative radius!
 
 // Initialize the walls, scene objects, and lightsource
 buildWalls();
 buildScene();
 
 // Set image width for coordinate conversion
 w_x=(W_RIGHT-W_LEFT);
 w_y=(W_BOTTOM-W_TOP);
 
 // READY - The loop below will generate rays in a way that suits the type of lightsource
 // defined in buildScene.c, and will initiate the propagation process.

#pragma omp parallel for schedule(dynamic,32) private(ray) 
 for (int i=0; i<num_rays; i++)
 {
  if (num_rays>10)
   if (i%(num_rays/10)==0) fprintf(stderr,"Progress=%f\n",(double)i/(double)(num_rays));
#ifdef __DEBUG_MODE
  fprintf(stderr,"Propagating ray %d of %d\n",i,num_rays);
#endif
  ray=makeLightSourceRay();
#ifdef __DEBUG_MODE
  fprintf(stderr,"Ray is at (%f,%f), direction (%f, %f)\n",ray.p.px,ray.p.py,ray.d.px,ray.d.py);
#endif
  propagateRay(&ray,0);
 }

 // Done with light propagation. Process the image array to create a final rendered image
 
 // First adjust gamma (simple log transform)
 for (int i=0; i<sx*sy*3; i++)
   *(imRGB+i)=log((*(imRGB+i))+1.5);
 
 im=(unsigned char *)calloc(sx*sy*3,sizeof(unsigned char));
 mx=-1;
 mi=10e15;
 for (int i=0; i<sx*sy*3; i++)
 {
   if (*(imRGB+i)<mi) mi=*(imRGB+i);
   if (*(imRGB+i)>mx) mx=*(imRGB+i);
 } 
 rng=mx-mi;
 fprintf(stderr,"Image range: mi=%f, mx=%f, range=%f\n",mi,mx,rng);
  
 for (int i=0; i<sx*sy*3; i++)
  *(im+i)=(unsigned char)(255.0*((*(imRGB+i)-mi)/rng));
 
#ifdef __DEBUG_MODE
 renderObjects();
#endif
 
 f=fopen("light2D_output.ppm","w");
 if (f!=NULL)
 {
  fprintf(f,"P6\n");
  fprintf(f,"# Output from Light2D.c\n");
  fprintf(f,"%d %d\n",sx,sy);
  fprintf(f,"255\n");
  fwrite(im,sx*sy*3*sizeof(unsigned char),1,f);
  fclose(f);
 }
 else fprintf(stderr,"Can not create output image file\n");
 
 // Release resources
 free(imRGB);
 free(im);
  
}
