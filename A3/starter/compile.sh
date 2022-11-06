#!/bin/sh
# g++ -O4 -g svdDynamic.c RayTracer.c utils.c -lm -o RayTracer
g++ -O4 -g -fopenmp svdDynamic.c RayTracer.c utils.c -lm -o RayTracer