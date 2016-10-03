#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

// include ispc's timing.h for timing
#include "timing.h"
// include ispc generated stub header
#include "spheres_ispc.h"
using namespace ispc;

typedef unsigned int unit;


int main() {
	printf("running!\n");
	const int width = 640;
	const int height = 480;
	
	Sphere *spheres = new Sphere[6];
	spheres[0].center[0] = 0.0;
	spheres[0].center[1] = -10004;
	spheres[0].center[2] = -20;
	spheres[0].radius = 10000;
	spheres[0].radius2 = 10000 * 10000;
	spheres[0].surfaceColor[0] = 0.20;
	spheres[0].surfaceColor[1] = 0.20;
	spheres[0].surfaceColor[2] = 0.20;
	spheres[0].reflection = 0.0;
	spheres[0].transparency = 0.0;
	spheres[0].emissionColor[0] = 0.0;
	spheres[0].emissionColor[1] = 0.0;
	spheres[0].emissionColor[2] = 0.0;

	spheres[1].center[0] = 0.0;
	spheres[1].center[1] = 0;
	spheres[1].center[2] = -20;
	spheres[1].radius = 4;
	spheres[1].radius2 = 16;
	spheres[1].surfaceColor[0] = 1;
	spheres[1].surfaceColor[1] = 0.32;
	spheres[1].surfaceColor[2] = 0.36;
	spheres[1].reflection = 1;
	spheres[1].transparency = 0.5;
	spheres[1].emissionColor[0] = 0.0;
	spheres[1].emissionColor[1] = 0.0;
	spheres[1].emissionColor[2] = 0.0;
	
	spheres[2].center[0] = 5;
	spheres[2].center[1] = -1;
	spheres[2].center[2] = -15;
	spheres[2].radius = 2;
	spheres[2].radius2 = 4;
	spheres[2].surfaceColor[0] = .90;
	spheres[2].surfaceColor[1] = 0.76;
	spheres[2].surfaceColor[2] = 0.46;
	spheres[2].reflection = 1;
	spheres[2].transparency = 0;
	spheres[2].emissionColor[0] = 0.0;
	spheres[2].emissionColor[1] = 0.0;
	spheres[2].emissionColor[2] = 0.0;
	
	spheres[3].center[0] = 5;
	spheres[3].center[1] = 0;
	spheres[3].center[2] = -25;
	spheres[3].radius = 3;
	spheres[3].radius2 = 9;
	spheres[3].surfaceColor[0] = 0.65;
	spheres[3].surfaceColor[1] = 0.77;
	spheres[3].surfaceColor[2] = 0.97;
	spheres[3].reflection = 1;
	spheres[3].transparency = 0;
	spheres[3].emissionColor[0] = 0.0;
	spheres[3].emissionColor[1] = 0.0;
	spheres[3].emissionColor[2] = 0.0;

	spheres[4].center[0] = -5.5;
	spheres[4].center[1] = 0;
	spheres[4].center[2] = -15;
	spheres[4].radius = 3;
	spheres[4].radius2 = 9;
	spheres[4].surfaceColor[0] = 0.90;
	spheres[4].surfaceColor[1] = 0.90;
	spheres[4].surfaceColor[2] = 0.90;
	spheres[4].reflection = 1;
	spheres[4].transparency = 0;
	spheres[4].emissionColor[0] = 0.0;
	spheres[4].emissionColor[1] = 0.0;
	spheres[4].emissionColor[2] = 0.0;
	
	// light
	spheres[5].center[0] = 0.0;
	spheres[5].center[1] = 20;
	spheres[5].center[2] = -30;
	spheres[5].radius = 3;
	spheres[5].radius2 = 9;
	spheres[5].surfaceColor[0] = 0.0;
	spheres[5].surfaceColor[1] = 0.0;
	spheres[5].surfaceColor[2] = 0.0;
	spheres[5].reflection = 0;
	spheres[5].transparency = 0;
	spheres[5].emissionColor[0] = 3.0;
	spheres[5].emissionColor[1] = 3.0;
	spheres[5].emissionColor[2] = 3.0;
	

	float *pixel_out_r = new float[width * height]; 
	float *pixel_out_g = new float[width * height]; 
	float *pixel_out_b = new float[width * height]; 
	reset_and_start_timer();
	render(pixel_out_r, pixel_out_g, pixel_out_b, width, height, spheres, 6);
	double dt = get_elapsed_mcycles();
	printf("@time of ISPC run:\t\t\t[%.3f] million cycles\n", dt);
	// print pixels
	for(int i = 0; i < width * height; i++){
		//printf("pixel %d: %f\n", i, pixel_out[i]);
	}

    // Save result to a PPM image (keep these flags if you compile under Windows)
    std::ofstream ofs("./untitled.ppm", std::ios::out | std::ios::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (unsigned i = 0; i < width * height; ++i) {
        ofs << (unsigned char)(std::min(float(1), pixel_out_r[i]) * 255) <<
               (unsigned char)(std::min(float(1), pixel_out_g[i]) * 255) <<
               (unsigned char)(std::min(float(1), pixel_out_b[i]) * 255);
    }   
    ofs.close();
    delete [] pixel_out_r;
    delete [] pixel_out_g;
    delete [] pixel_out_b;

	
	return 0;
}

	
	
