#include <cmath>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

// include ispc's timing.h for timing
#include "../ispc_spheres/timing.h"
// include rand_sphere for generating spheres
#include "../common/rand_sphere.h"
#define MAX_RAY_DEPTH 5


struct Sphere {
    float center[3];
    float radius, radius2;
    float transparency, reflection;
    float emissionColor[3], surfaceColor[3];
};


__device__ void Normalize(float3 &v1) {
    float nor2 = v1.x * v1.x + v1.y * v1.y + v1.z * v1.z;
    if (nor2 > 0) {
        float invNor = 1 / sqrt(nor2);
        v1.x *= invNor, v1.y *= invNor, v1.z *= invNor;
    }
}

__device__ float3 Cross(const float3 v1, const float3 v2) {
    float v1x = v1.x, v1y = v1.y, v1z = v1.z;
    float v2x = v2.x, v2y = v2.y, v2z = v2.z;
    float3 ret;
    ret.x = (v1y * v2z) - (v1z * v2y);
    ret.y = (v1z * v2x) - (v1x * v2z);
    ret.z = (v1x * v2y) - (v1y * v2x);
    return ret;
}

__device__ float Dot(const float3 a, const float3 b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

__device__ float mix(const float &a, const float &b, const float &mix) {
    return b * mix + a * (1 - mix);
}

__device__ float3 operator+(const float3 &a, const float3 &b) {
	return make_float3(a.x+b.x, a.y+b.y, a.z+b.z);
}

__device__ float3 operator-(const float3 &a, const float3 &b) {
	return make_float3(a.x-b.x, a.y-b.y, a.z-b.z);
}

__device__ float3 operator*(const float3 &a, const float3 &b) {
	return make_float3(a.x*b.x, a.y*b.y, a.z*b.z);
}

__device__ float3 operator*(const float3 &a, const float &b) {
	return make_float3(a.x*b, a.y*b, a.z*b);
}

__device__ bool SphereIntersect(const Sphere &sphere,
                            const float3 &rayorig,
                            const float3 &raydir,
                            float &t0, float &t1)  {
    float3 center = make_float3(sphere.center[0], sphere.center[1], sphere.center[2]);
    float3 l = center - rayorig;
    float tca = Dot(l, raydir);
    if (tca < 0) return false;
    float d2 = Dot(l, l) - tca * tca;
    if (d2 > sphere.radius2) return false;
    float thc = sqrt(sphere.radius2 - d2);
    t0 = tca - thc;
    t1 = tca + thc;
    return true;
}


 __device__ float3 trace(const float3 &rayorig, const float3 &raydir, const Sphere *spheres, const int depth, const int size) {
    float tnear = INFINITY;
    Sphere sphere;
    bool intersect_found = false;
    for (int i = 0; i < size; ++i) {
        float t0 = INFINITY, t1 = INFINITY;
        if (SphereIntersect(spheres[i], rayorig, raydir, t0, t1)) {
            if (t0 < 0) t0 = t1;
            if (t0 < tnear) {
                tnear = t0;
                sphere = spheres[i];
                intersect_found = true;
            }
        }
    }
    // if there's no intersection return black or background color
    if (!intersect_found) {
        float3 ret = make_float3(2, 2, 2);
        return ret;
    }
    float3 surfaceColor = make_float3(0, 0, 0);
    float3 phit = rayorig + raydir * tnear; // point of intersection
    float3 center = make_float3(sphere.center[0], sphere.center[1], sphere.center[2]);
    float3 nhit = phit - center; // normal at the insection point

    Normalize(nhit); // normalize the direction
    // If the normal and the view direction are not opposite to each other
    // reverse the normal direction. That also means we are inside the sphere so set
    // the inside bool to true. Finally reverse the sign of IdotN which we want
    // positive.
    float bias = 1e-4;
    bool inside = false;

    if (Dot(raydir, nhit) > 0) nhit = nhit * -1.0, inside = true;
    if ((sphere.transparency > 0 || sphere.reflection > 0) && depth < MAX_RAY_DEPTH) {
        float facingratio = -(Dot(raydir, nhit));
        // change the mix value to the tweak effect
        const float fresneleffect = mix(powf(1 - facingratio, 3), 1.0, 0.1);
        // compute reflection direction (no need to normalize because all vectors
        // are already normalized)
        float3 refldir = raydir - nhit * 2 * Dot(raydir, nhit);
        Normalize(refldir);
        float3 reflection = trace(phit + nhit * bias, refldir, spheres, depth + 1, size);
        float3 refraction = make_float3(0,0,0);
        // if sphere is also transparent, compute refraction ray
        if (sphere.transparency) {
            float ior = 1.1;
            float eta;
            if (inside) eta = ior;
            else eta = 1 / ior; 
            float cosi = -(Dot(nhit, raydir));
            float k = 1 - eta * eta * (1 - cosi * cosi);
            float3 refrdir = raydir * eta + nhit * (eta * cosi - sqrtf(k));
            Normalize(refrdir);
            refraction = trace(phit - nhit * bias, refrdir, spheres, depth + 1, size);
        }
        // the result is a mix of reflection and refraction 
        float3 sphereSurfaceColor = make_float3(sphere.surfaceColor[0],
                                     			sphere.surfaceColor[1],
                                     			sphere.surfaceColor[2]);
        surfaceColor = (
            reflection * fresneleffect +
            refraction * (1 - fresneleffect) * sphere.transparency) * sphereSurfaceColor;
    }
    else {
        // it's a diffuse object, no need to raytace any further
        for (int i = 0; i < size; ++i) {
            if (spheres[i].emissionColor[0] > 0) {
                // this is a light
                float3 transmission = make_float3(1,1,1);
                float3 centeri = make_float3(spheres[i].center[0],
                                  			 spheres[i].center[1],
                                  		     spheres[i].center[2]);
                float3 lightDirection = centeri - phit;
                Normalize(lightDirection);
                for (int j = 0; j < size; ++j) {
                    if (i != j) {
                        float t0, t1;
                        if (SphereIntersect(spheres[j], phit + nhit * bias,
                                                    lightDirection, t0, t1)) {
                            transmission = make_float3(0,0,0);
                            break;
                        }
                    }
                }
                float3 sphereSurfaceColor = make_float3(
											 sphere.surfaceColor[0],
                                             sphere.surfaceColor[1],
                                             sphere.surfaceColor[2]);

                float3 sphereEmissionColor = make_float3(
                                             spheres[i].emissionColor[0],
                                             spheres[i].emissionColor[1],
                                             spheres[i].emissionColor[2]);
                surfaceColor = surfaceColor + sphereSurfaceColor * transmission * max(0., Dot(nhit, lightDirection)) * sphereEmissionColor;
            }
        }
    }
    float3 sphereEmissionColor = make_float3(
                             sphere.emissionColor[0],
                             sphere.emissionColor[1],
                             sphere.emissionColor[2]);
    return surfaceColor + sphereEmissionColor;
}


__global__ void render(float *gpu_pixels_r,
					   float *gpu_pixels_g,
				       float *gpu_pixels_b,
					   const int width,
					   const int height,
					   Sphere * spheres,
					   const int sphere_count,
					   const float invWidth,
					   const float invHeight,
					   const float fov,
					   const float aspectratio,
					   const float angle){
    int pos_x = (blockIdx.x * blockDim.x) + threadIdx.x;
    int pos_y = (blockIdx.y * blockDim.y) + threadIdx.y;
    if(pos_x >= width || pos_y >= height)
        return;        
    int pixel_index = pos_x + pos_y * width;
	float xx = (2 * ((pos_x + 0.5) * invWidth) - 1) * angle * aspectratio;
	float yy = (1 - 2 * ((pos_y + 0.5) * invHeight)) * angle;
	float3 raydir = make_float3(xx, yy, -1);
	Normalize(raydir);
	float3 ray_orig = make_float3(0, 0, 0);
	float3 pixel_color = trace(ray_orig, raydir, spheres, 0, sphere_count);
	gpu_pixels_r[pixel_index] = pixel_color.x;
	gpu_pixels_g[pixel_index] = pixel_color.y;
	gpu_pixels_b[pixel_index] = pixel_color.z;
}


int main(int argc, char **argv) {
    // parse command line args
    if (argc != 4 && argc != 1) {
        std::cerr << "Usage: " << argv[0] << " <width> <height> <sphere_count>" << std::endl;
        return EXIT_FAILURE;
    }

	Sphere *cpuSpheres;
	Sphere *gpuSpheres;
	int sphereSize = sizeof(Sphere);
    int width, height, numSpheres, numBytes;
		
	// Define spheres
    if (argc == 1) { 
	    // Allocate memory for Sphere arrays
        numSpheres = 6;
        numBytes = numSpheres * sphereSize;
        cpuSpheres = (Sphere*)malloc(numBytes);
        cudaMalloc((void**)&gpuSpheres, numBytes);

        cpuSpheres[0].center[0] = 0.0;
        cpuSpheres[0].center[1] = -10004;
        cpuSpheres[0].center[2] = -20;
        cpuSpheres[0].radius = 10000;
        cpuSpheres[0].radius2 = 10000 * 10000;
        cpuSpheres[0].surfaceColor[0] = 0.20;
        cpuSpheres[0].surfaceColor[1] = 0.20;
        cpuSpheres[0].surfaceColor[2] = 0.20;
        cpuSpheres[0].reflection = 0.0;
        cpuSpheres[0].transparency = 0.0;
        cpuSpheres[0].emissionColor[0] = 0.0;
        cpuSpheres[0].emissionColor[1] = 0.0;
        cpuSpheres[0].emissionColor[2] = 0.0;

        cpuSpheres[1].center[0] = 0.0;
        cpuSpheres[1].center[1] = 0;
        cpuSpheres[1].center[2] = -20;
        cpuSpheres[1].radius = 4;
        cpuSpheres[1].radius2 = 16;
        cpuSpheres[1].surfaceColor[0] = 1;
        cpuSpheres[1].surfaceColor[1] = 0.32;
        cpuSpheres[1].surfaceColor[2] = 0.36;
        cpuSpheres[1].reflection = 1;
        cpuSpheres[1].transparency = 0.5;
        cpuSpheres[1].emissionColor[0] = 0.0;
        cpuSpheres[1].emissionColor[1] = 0.0;
        cpuSpheres[1].emissionColor[2] = 0.0;
        
        cpuSpheres[2].center[0] = 5;
        cpuSpheres[2].center[1] = -1;
        cpuSpheres[2].center[2] = -15;
        cpuSpheres[2].radius = 2;
        cpuSpheres[2].radius2 = 4;
        cpuSpheres[2].surfaceColor[0] = .90;
        cpuSpheres[2].surfaceColor[1] = 0.76;
        cpuSpheres[2].surfaceColor[2] = 0.46;
        cpuSpheres[2].reflection = 1;
        cpuSpheres[2].transparency = 0;
        cpuSpheres[2].emissionColor[0] = 0.0;
        cpuSpheres[2].emissionColor[1] = 0.0;
        cpuSpheres[2].emissionColor[2] = 0.0;
        
        cpuSpheres[3].center[0] = 5;
        cpuSpheres[3].center[1] = 0;
        cpuSpheres[3].center[2] = -25;
        cpuSpheres[3].radius = 3;
        cpuSpheres[3].radius2 = 9;
        cpuSpheres[3].surfaceColor[0] = 0.65;
        cpuSpheres[3].surfaceColor[1] = 0.77;
        cpuSpheres[3].surfaceColor[2] = 0.97;
        cpuSpheres[3].reflection = 1;
        cpuSpheres[3].transparency = 0;
        cpuSpheres[3].emissionColor[0] = 0.0;
        cpuSpheres[3].emissionColor[1] = 0.0;
        cpuSpheres[3].emissionColor[2] = 0.0;

        cpuSpheres[4].center[0] = -5.5;
        cpuSpheres[4].center[1] = 0;
        cpuSpheres[4].center[2] = -15;
        cpuSpheres[4].radius = 3;
        cpuSpheres[4].radius2 = 9;
        cpuSpheres[4].surfaceColor[0] = 0.90;
        cpuSpheres[4].surfaceColor[1] = 0.90;
        cpuSpheres[4].surfaceColor[2] = 0.90;
        cpuSpheres[4].reflection = 1;
        cpuSpheres[4].transparency = 0;
        cpuSpheres[4].emissionColor[0] = 0.0;
        cpuSpheres[4].emissionColor[1] = 0.0;
        cpuSpheres[4].emissionColor[2] = 0.0;
    }
    else {
        // Read in argv
        height = std::atoi(argv[1]);
        width = std::atoi(argv[2]);
        numSpheres = std::atoi(argv[3]);
	    // Allocate memory for Sphere arrays
        numBytes = numSpheres * sphereSize;
        cpuSpheres = (Sphere*)malloc(numBytes);
        cudaMalloc((void**)&gpuSpheres, numBytes);
        // Generate random spheres within frustum
        srand(time(NULL));
        float x,y,z;
        int count = 0;
	    while (count < numSpheres - 1) {
            rand_sphere(x, y, z, height, width, 2);
            cpuSpheres[count].center[0] = x;
            cpuSpheres[count].center[1] = y;
            cpuSpheres[count].center[2] = z;
            cpuSpheres[count].radius = 2;
            cpuSpheres[count].radius2 = 4;
            cpuSpheres[count].surfaceColor[0] = random_float(0,1);
            cpuSpheres[count].surfaceColor[1] = random_float(0,1);
            cpuSpheres[count].surfaceColor[2] = random_float(0,1);
            cpuSpheres[count].reflection = 0.5;
            cpuSpheres[count].transparency = 0.5;
            cpuSpheres[count].emissionColor[0] = 0.0;
            cpuSpheres[count].emissionColor[1] = 0.0;
            cpuSpheres[count].emissionColor[2] = 0.0;
            count++;
        }
    }   
    

	// light
    int last = numSpheres - 1;
	cpuSpheres[last].center[0] = 0.0;
	cpuSpheres[last].center[1] = 20;
	cpuSpheres[last].center[2] = -30;
	cpuSpheres[last].radius = 3;
	cpuSpheres[last].radius2 = 9;
	cpuSpheres[last].surfaceColor[0] = 0.0;
	cpuSpheres[last].surfaceColor[1] = 0.0;
	cpuSpheres[last].surfaceColor[2] = 0.0;
	cpuSpheres[last].reflection = 0;
	cpuSpheres[last].transparency = 0;
	cpuSpheres[last].emissionColor[0] = 3.0;
	cpuSpheres[last].emissionColor[1] = 3.0;
	cpuSpheres[last].emissionColor[2] = 3.0;

	
	// Allocate memory for pixel arrays
	float *pixel_out_r = new float[width * height]; 
	float *pixel_out_g = new float[width * height]; 
	float *pixel_out_b = new float[width * height]; 
	float *gpu_pixel_r;	
	float *gpu_pixel_g;
	float *gpu_pixel_b;
	cudaMalloc((void**)&gpu_pixel_r, sizeof(float)*width*height); 
	cudaMalloc((void**)&gpu_pixel_g, sizeof(float)*width*height); 
	cudaMalloc((void**)&gpu_pixel_b, sizeof(float)*width*height); 

	// Derive display variables
	float invWidth = 1 / (float)width;
	float invHeight = 1 / (float)height;
    float fov = 30;
    float aspectratio = width / (float)height;
    float angle = tan(M_PI * 0.5 * fov / 180.);
	
	reset_and_start_timer();
	// Copy cpuSpheres to gpuSpheres
	cudaMemcpy(gpuSpheres, cpuSpheres, numBytes, cudaMemcpyHostToDevice); 
	 
	const dim3 blockSize(24,24,1);
	const dim3 gridSize(width/24+1, height/24+1, 1);
	render<<<gridSize,blockSize>>>(gpu_pixel_r, gpu_pixel_g, gpu_pixel_b, width, height, gpuSpheres, numSpheres, invWidth, invHeight, fov, aspectratio, angle);

	// Copy arrays back to cpu
	cudaMemcpy(pixel_out_r, gpu_pixel_r, sizeof(float)*width*height, cudaMemcpyDeviceToHost);
	cudaMemcpy(pixel_out_g, gpu_pixel_g, sizeof(float)*width*height, cudaMemcpyDeviceToHost);
	cudaMemcpy(pixel_out_b, gpu_pixel_b, sizeof(float)*width*height, cudaMemcpyDeviceToHost);
	// Stop timer
	double dt = get_elapsed_mcycles();
	printf("@time of CUDA run:\t\t\t[%.3f] million cycles\n", dt);

	
    // Save result to a PPM image (keep these flags if you compile under Windows)
    std::ofstream ofs("./CUDA_spheres.ppm", std::ios::out | std::ios::binary);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (unsigned i = 0; i < width * height; ++i) {
        ofs << (unsigned char)(std::min(float(1), pixel_out_r[i]) * 255) <<
               (unsigned char)(std::min(float(1), pixel_out_g[i]) * 255) <<
               (unsigned char)(std::min(float(1), pixel_out_b[i]) * 255);
    }   
    ofs.close();

	// Deallocate memory
    delete [] pixel_out_r;
    delete [] pixel_out_g;
    delete [] pixel_out_b;
	free(cpuSpheres);
	cudaFree(gpuSpheres);
	cudaFree(gpu_pixel_r);
	cudaFree(gpu_pixel_g);
	cudaFree(gpu_pixel_b);
	
	return 0;
}

	
	
