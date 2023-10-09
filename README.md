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
    <li><a href="#adv-renderer">Advanced Renderer</a></li>
    <li><a href="#path-tracing">Path Tracer</a></li>
  </ol>

## About The Project
This is a project to showcase our understanding of 3D graphics and maths, ray-tracing techniques, and other simulations to generate nice-looking scenes. 

## Usage
You will need to have an environment that can compile **C** language.(Linux or minGW)
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
This Demo shows a 1024 * 1024 scene with a point light source at the center with 12 reflective circles and 12 refractive circles around, light rays will diffuse on walls.  

![2D Simulation Demo](https://github.com/TaihouAnF/Basic-Rendering-Engine/blob/main/Demo/light2D_output.png)

## 3D Renderer

This is where the fun begins.  
This feature projects the 3D space into a 2D scene by utilizing the [**Backward Ray Tracing**](https://en.wikipedia.org/wiki/Ray_tracing_(graphics)) technique.  
During the process of generation, the engine would shoot a ray for each pixel on the 2D screen/plane, check whether it collides with any objects in the 3D space, and calculate the RGB value for that pixel once the ray terminated.  
The ray would terminate under: 1. it hits nothing; 2. it reaches the limit of recursion time.  
And for the calculation of RGB, we used a local illumination model called [**"Phong Model"**](https://en.wikipedia.org/wiki/Phong_reflection_model).  

### Demo

![3D Raytracer Demo](https://github.com/TaihouAnF/Basic-Rendering-Engine/blob/main/Demo/full.png)
