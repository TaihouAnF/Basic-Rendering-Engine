<a name="readme-top"></a>
<div align="center">
  <h1 align="center">2D/3D Rendering Engine</h1>
  <p>
    A 2D & 3D rendering engine supporting basic raytracing and other features
  </p>
</div>

## Table of Contents</summary>
  <ol>
    <li><a href="#about-the-project">About The Project</a></li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#2d-renderer">2D Renderer</a></li>
    <li><a href="#3d-renderer">3D Renderer</a></li>
    <li><a href="#advanced-renderer">Advanced Renderer</a></li>
    <li><a href="#path-tracer">Path Tracer</a></li>
  </ol>

## About The Project
This is a project to showcase our understanding of 3D graphics and maths, ray-tracing techniques, and other simulations to generate nice-looking scenes. 

## Usage
You will need to have an environment that can compile **C** language. (Linux or minGW)
There are compile shell scripts for each feature, you can directly run the script or simply type the script yourself and run it.
To run a feature and see the demo:
* 2D Renderer
  ```sh
  g++ -g -O3 -fopenmp light2D.c -lm -o light2D
  ```
  
* 3D Renderer
  ```sh
  g++ -O4 -g svdDynamic.c RayTracer.c utils.c -lm -o RayTracer
  ```
  
* Advanced Renderer
  ```sh
  g++ -O4 -g -fopenmp svdDynamic.c RayTracer.c utils.c -lm -o RayTracer
  ```
  
* Path Tracer
  ```sh
  g++ -O3 -g -fopenmp svdDynamic.c PathTracer.c utils_path.c -lm -o PathTracer
  ```
  

## 2D Renderer

The first feature of this project is a simulation of light rays in 2D.  
It would have at least one single-point light source for each scene. And in the scene, we can set up multiple round objects with different properties.  
The Renderer supports **light ray emission, diffusion, perfect reflection, and refraction**.

### Demo
This Demo shows a *1024 * 1024* scene of:
  * a point light source at the center;
  * 12 reflective circles and
  * 12 refractive circles around;
  * light rays will diffuse on walls;  

<div align="center">

  ![2D Simulation Demo](https://github.com/TaihouAnF/Basic-Rendering-Engine/blob/main/Demo/light2D_output.png)

</div>

## 3D Renderer

This is where the fun begins.  
This feature projects the 3D space into a 2D scene by utilizing the [**Backward Ray Tracing**](https://en.wikipedia.org/wiki/Ray_tracing_(graphics)) technique.  
During the process of generation, the engine would shoot a ray for each pixel on the 2D screen/plane, check whether it collides with any objects in the 3D space, and calculate the RGB value for that pixel once the ray terminated.  
The ray would terminate under: 1. it hits nothing; 2. it reaches the limit of recursion time.  
And for the calculation of RGB, we used a local illumination model called [**"Phong Model"**](https://en.wikipedia.org/wiki/Phong_reflection_model).  

### Example
This Demo shows a *512 * 512* scene of:
* one point light source behind the camera,
* together with two affinely transform sphere with different factor of reflection, diffsusion & RGBs,
* and one full diffusive plane(no reflection at all).  

As you can see, the engine supports specular reflection, diffuse, and ambient lighting, it also supports **recusive reflection**.

<div align="center">
  
  ![3D Raytracer Demo](https://github.com/TaihouAnF/Basic-Rendering-Engine/blob/main/Demo/full.png)
  
</div>


## Advanced Renderer

This feature is based on the basic [**3D Renderer**](#3d-renderer), we modified and added implementation to supports new features, including:
  *  Multi-Threading To Boost The Performance
  *  Random Sampling Light Ray Directions
  *  Anti-Aliasing
  *  Area Light Source & Soft Shadows 
  *  Mapping:
      *  Alpha Mappings
      *  Texture Mappings
      *  Normal Mappings
  *  Refraction on Arbitrarily Nubmer of Nested Refractive Objects (An Object inside on another Object)  

We were planned to finish other feature but we didn't have enough time for those, which includes but not limited to:
  *  Depth Of Field Effect(WIP)
  *  Octree Structure & Ray collide Detection
  *  Photon Mapping
  *  And Others

### Example

<div align="center">

  ![Advanced_Renderer_Demo](https://github.com/TaihouAnF/Basic-Rendering-Engine/blob/main/Demo/Advanced.png)

</div>

## Path Tracer

The last feature of this engine would be the [**Path Tracer**](https://en.wikipedia.org/wiki/Path_tracing) which would sample light rays in each pixel, then use the average value of RGBs of each light ray as the final color for that pixel.  
It will contians *noises* as we're sampling light rays with arbitrary direction would cannot guarantee each light ray can hit a light source.

### Example

#### Cornell Box
<div align="center">

  ![Cornell Box](https://github.com/TaihouAnF/Basic-Rendering-Engine/blob/main/Demo/Cornell_IS_ES.png)
  
</div>

#### Final Demo
<div align="center">

  ![Final Demo](https://github.com/TaihouAnF/Basic-Rendering-Engine/blob/main/Demo/WE_HAVE_CONQUERED_CG.png)
  
</div>

