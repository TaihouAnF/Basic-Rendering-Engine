 // Sets up all objects in the scene. This involves creating each object,
 // defining the transformations needed to shape and position it as
 // desired, specifying the reflectance properties (albedos and colours)
 // and setting up textures where needed.
 // Light sources must be defined, positioned, and their colour defined.
 // All objects must be inserted in the object_list. All light sources
 // must be inserted in the light_list.
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
 struct pointLS *l;
 struct point3D p;

 // Simple scene for Assignment 3:
 // Insert a couple of objects. A plane and two spheres
 // with some transformations.

 // Note the parameters: ra, rd, rs, rg, R, G, B, alpha, r_index, and shinyness)

o=newSphere(.05,.95,.95,.75,1,1,1,0.1,2,6);
Translate(o,0,0,1);
invert(&o->T[0][0],&o->Tinv[0][0]);
insertObject(o,&object_list);

o=newSphere(.05,.95,.95,.75,1,0,1,0.5,1.5,6);
Translate(o,-2.5,2.5,2);
invert(&o->T[0][0],&o->Tinv[0][0]);
insertObject(o,&object_list);

o=newSphere(.05,.95,.95,.75,1,0,1,0.5,1.5,6);
Translate(o,2.5,-2.5,2);
invert(&o->T[0][0],&o->Tinv[0][0]);
insertObject(o,&object_list);

o=newSphere(.05,.95,.95,.75,0,1,0,0.5,1.5,6);
Translate(o,2.5,2.5,2);
invert(&o->T[0][0],&o->Tinv[0][0]);
insertObject(o,&object_list);

o=newSphere(.05,.95,.95,.75,0,1,0,0.5,1.5,6);
Translate(o,-2.5,-2.5,2);
invert(&o->T[0][0],&o->Tinv[0][0]);
insertObject(o,&object_list);

o=newPlane(.05,.75,.05,.05,1,0,0,1,1,2);
Scale(o,11,11,11);
Translate(o,0,0,6);
invert(&o->T[0][0],&o->Tinv[0][0]);
loadTexture(o, "./Texture/textureG.ppm", 1, &texture_list);
insertObject(o,&object_list);

 // add a cylinder
o=newCyl(.05,.95,.95,.75,.75,.95,.55,1,1,6);
RotateX(o,PI/2);
Translate(o,0,10,-5);
invert(&o->T[0][0],&o->Tinv[0][0]);
o->isLightSource = 1;
insertObject(o,&object_list);

 // End of simple scene for Assignment 2
 // Keep in mind that you can define new types of objects such as cylinders and parametric surfaces,
 // or, you can create code to handle arbitrary triangles and then define objects as surface meshes.
 //
 // Remember: A lot of the quality of your scene will depend on how much care you have put into defining
 //           the relflectance properties of your objects, and the number and type of light sources
 //           in the scene.
 
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////
 // TO DO: For Assignment 3 you *MUST* define your own cool scene.
 //	   We will be looking for the quality of your scene setup, the use of hierarchical or composite
 //	   objects that are more interesting than the simple primitives from A2, the use of textures
 //        and other maps, illumination and illumination effects such as soft shadows, reflections and
 //        transparency, and the overall visual quality of your result. Put some work into thinking
 //        about these elements when designing your scene.
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////
