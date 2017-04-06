#include <math.h>

#define PI 3.14159265


float random_float(const float min, const float max) {
	return (max - min) * (static_cast <float> (rand()) / RAND_MAX) + min;
}

void rand_sphere(float &x, float &y, float &z,
                 const int height, const int width, const int radius) {
	float aspect = static_cast <float> (width)/height;
	z = random_float(20, 70);
	float frustum_h = (z * tan(30*0.5*(PI/180))) - static_cast <float> (radius);
	float frustum_w = frustum_h * aspect;
	x = random_float(-(frustum_w), frustum_w);
	y = random_float(-(frustum_h), frustum_h);
	z = -z;
}	
	
