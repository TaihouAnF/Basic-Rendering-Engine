void buildScene(void)
{
 // Sets up all objets in the scene. This involves creating each object,
 // defining the transformations needed to shape and position it as
 // desired, specifying the reflectance properties (albedos and colours)
 // and setting up textures where needed.
 //
 // NOTE: Light sources are now EXCLUSIVELY area light sources. They are
 //       defined by regular objects whose 'isLightSource' flag is set
 //       to 1. Therefore, you can create light sources with any shape
 //       and colour using the same object primitives and transforms
 //       you're using to set up the scene.
 //
 // To create hierarchical objects:
 //    You must keep track of transformations carried out by parent objects
 //    as you move through the hierarchy. Declare and manipulate your own
 //    transformation matrices (use the provided functions in utils.c to
 //    compound transformations on these matrices). When declaring a new
 //    object within the hierarchy
 //    - Initialize the object
 //    - Apply any object-level transforms to shape/rotate/resize/move
 //      the object using regular object transformation functions
 //    - Apply the transformations passed on from the parent object
 //      by pre-multiplying the matrix containing the parent's transforms
 //      with the object's own transformation matrix.
 //    - Compute and store the object's inverse transform as usual.
 //
 // NOTE: After setting up the transformations for each object, don't
 //       forget to set up the inverse transform matrix!

struct object3D *o;
struct point3D p;
// Cornell box
int objnum = 10;
for (int j=0; j < objnum; j++) {
    // x between -5 and 5
    // y between -5 and 5
    // z between 0 and 10
    float x = (rand() % 1000) / 100.0 - 5.0;
    float y = (rand() % 1000) / 100.0 - 5.0;
    float z = 10 + j * 5;
    // random rgb value
    float r = (rand() % 1000) / 1000.0;
    float g = (rand() % 1000) / 1000.0;
    float b = (rand() % 1000) / 1000.0;
    // generate 0 and 1 and 2
    int type = rand() % 3;

    if (type == 0) {
        o=newSphere(1.0,0.0,0.0,r,g,b,.05,2.47);		// Diffuse
        Translate(o,x,y,z);
        invert(&o->T[0][0],&o->Tinv[0][0]);
        insertObject(o,&object_list, LS_num);
    } else if (type == 1) {
        o=newSphere(0.0,1.0,0.0,r,g,b,.05,2.47);		// Reflective
        Translate(o,x,y,z);
        invert(&o->T[0][0],&o->Tinv[0][0]);
        insertObject(o,&object_list, LS_num);
    } else {
        o=newSphere(0.0,0.0,1.0,r,g,b,.05,2.47);		// Refractive
        Translate(o,x,y,z);
        invert(&o->T[0][0],&o->Tinv[0][0]);
        insertObject(o,&object_list, LS_num);
    }
}

o=newPlane(1.00,0.00,0.0,.25,.25,.75,0.0,1.54);
Scale(o,45,45,45);
RotateZ(o,PI/2);
Translate(o,0,5,40);
invert(&o->T[0][0],&o->Tinv[0][0]);
loadTexture(o, "./Texture/skybox.ppm", 1, &texture_list);
o->isLightSource=1;
insertObject(o,&object_list, LS_num);

o=newSphere(1.00,0.00,0.0,1.0,1.0,1.0,0.0,1.54);
Translate(o,0,10,5);
invert(&o->T[0][0],&o->Tinv[0][0]);
o->isLightSource=1;
insertObject(o,&object_list, LS_num);

int stars = 120;
for (int i=0; i<50;i++) {
    // Planar light source at top
    o=newSphere(1.00,0.00,0.0,1.0,1.0,1.0,0.0,1.54);
    //generate y between 10 and 20
    float y = (rand() % 1000) / 100.0 + 10.0;
    Translate(o,0,10,5);
    invert(&o->T[0][0],&o->Tinv[0][0]);
    o->isLightSource=1;
    insertObject(o,&object_list, LS_num);
}

for (int i=0; i<stars;i++) {
    // Planar light source at top
    o=newSphere(1.00,0.00,0.0,1.0,1.0,1.0,0.0,1.54);
    //generate x between 10 and 20
    float x = (rand() % 1000) / 100.0 + 10.0;
    //generate y between -20 and 20
    float y = (rand() % 50) - 50;
    Scale(o,0.1,0.1,0.1);
    Translate(o,x,y+15,10);
    invert(&o->T[0][0],&o->Tinv[0][0]);
    o->isLightSource=1;
    insertObject(o,&object_list, LS_num);
}

}