#!/usr/bin/env bash
reset
ispc -O3 --target=avx2-i32x16 raytracer.ispc -o objs/raytracer_ispc.o -h objs/raytracer_ispc.h
g++ raytracer2.cpp -Iobjs/ -O3 -Wall -c -o objs/raytracer2.o
g++ -Iobjs/ -O3 -Wall -o raytracer objs/raytracer2.o objs/raytracer_ispc.o 
