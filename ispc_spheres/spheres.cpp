#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

// include rand_sphere for sphere generation
#include "../common/rand_sphere.h"
// include ispc's timing.h for timing
#include "timing.h"
// include ispc generated stub header
#include "spheres_ispc.h"

using namespace ispc;

typedef unsigned int unit;


int main(int argc, char **argv) {

    // parse command line args 
    if (argc != 4 && argc != 1) {  
        std::cerr << "Usage: " << argv[0] << " <width> <height> <sphere_count>" << std::endl; 
        return EXIT_FAILURE;  
    } 
	
	int width, height, sphere_count;
    Sphere *spheres;
	
	if (argc == 1) {	
        width = 640;
        height = 480;
        sphere_count = 6;
        spheres = new Sphere[sphere_count];
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
    }	
    else {
        srand(time(NULL));
        height = std::atoi(argv[1]);
        width = std::atoi(argv[2]);
        sphere_count = std::atoi(argv[3]) + 1;
        spheres = new Sphere[sphere_count];
        unsigned radius = 2;
        float x,y,z;
        int count = 0;
        // generate spheres within frustrum
        while (count < sphere_count - 1) {
            rand_sphere(x, y, z, height, width, radius);
            spheres[count].center[0] = x;
            spheres[count].center[1] = y;
            spheres[count].center[2] = z;
            spheres[count].radius = 2;
            spheres[count].radius2 = 4;
            spheres[count].surfaceColor[0] = random_float(0,1);
            spheres[count].surfaceColor[1] = random_float(0,1);
            spheres[count].surfaceColor[2] = random_float(0,1);
            spheres[count].reflection = 0.5;
            spheres[count].transparency = 0.5;
            spheres[count].emissionColor[0] = 0.0;
            spheres[count].emissionColor[1] = 0.0;
            spheres[count].emissionColor[2] = 0.0;
            count++;
        }
    }
    
	// Add light sphere
    int last = sphere_count-1;
	spheres[last].center[0] = 0.0;
	spheres[last].center[1] = 20;
	spheres[last].center[2] = -30;
	spheres[last].radius = 3;
	spheres[last].radius2 = 9;
	spheres[last].surfaceColor[0] = 0.0;
	spheres[last].surfaceColor[1] = 0.0;
	spheres[last].surfaceColor[2] = 0.0;
	spheres[last].reflection = 0;
	spheres[last].transparency = 0;
	spheres[last].emissionColor[0] = 3.0;
	spheres[last].emissionColor[1] = 3.0;
	spheres[last].emissionColor[2] = 3.0;
	
    // Initialize rgb arrays
	float *pixel_out_r = new float[width * height]; 
	float *pixel_out_g = new float[width * height]; 
	float *pixel_out_b = new float[width * height]; 
	reset_and_start_timer();
	render(pixel_out_r, pixel_out_g, pixel_out_b, width, height, spheres, sphere_count);
	double dt = get_elapsed_mcycles();
	printf("@time of ISPC run:\t\t\t[%.3f] million cycles\n", dt);
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

	
	
