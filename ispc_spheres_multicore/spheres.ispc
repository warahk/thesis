// ispc equivelent of scratchpixel example raytacers

typedef float<3> float3;

static inline void Normalize(float3 &v1) {
	float nor2 = v1.x * v1.x + v1.y * v1.y + v1.z * v1.z;
	if (nor2 > 0) {
		float invNor = 1 / sqrt(nor2);
		v1.x *= invNor, v1.y *= invNor, v1.z *= invNor;
	}
}

static inline float3 Cross(const float3 v1, const float3 v2) {
    float v1x = v1.x, v1y = v1.y, v1z = v1.z;
    float v2x = v2.x, v2y = v2.y, v2z = v2.z;
    float3 ret;
    ret.x = (v1y * v2z) - (v1z * v2y);
    ret.y = (v1z * v2x) - (v1x * v2z);
    ret.z = (v1x * v2y) - (v1y * v2x);
    return ret;
}

static inline float Dot(const float3 a, const float3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}


struct Sphere {
	float center[3];
	float radius, radius2;
	float transparency, reflection;
    float emissionColor[3], surfaceColor[3];
};

		
static bool SphereIntersect(const uniform Sphere &sphere, 
							const float3 &rayorig, 
							const float3 &raydir, 
							float &t0, float &t1)  {
	float3 center = {sphere.center[0], sphere.center[1], sphere.center[2]};
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

#define MAX_RAY_DEPTH 5
#define INFINITY 999999999
static inline float mix(const float &a, const uniform float &b, const uniform float &mix) {
	return b * mix + a * (1 - mix);
}


float3 trace(const float3 &rayorig, const float3 &raydir, const uniform Sphere spheres[], const uniform int depth, uniform int size) {
	float tnear = INFINITY;
	Sphere sphere;
	bool intersect_found = false;
	for (uniform int i = 0; i < size; ++i) {
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
		float3 ret = {2, 2, 2};
		return ret;
	}
	float3 surfaceColor = {0, 0, 0};
	float3 phit = rayorig + raydir * tnear; // point of intersection
	float3 center = {sphere.center[0], sphere.center[1], sphere.center[2]};
	float3 nhit = phit - center; // normal at the insection point

	Normalize(nhit); // normalize the direction
	// If the normal and the view direction are not opposite to each other
    // reverse the normal direction. That also means we are inside the sphere so set
    // the inside bool to true. Finally reverse the sign of IdotN which we want
    // positive.
	float bias = 1e-4;
	bool inside = false;

	if (Dot(raydir, nhit) > 0) nhit = -nhit, inside = true;
	if ((sphere.transparency > 0 || sphere.reflection > 0) && depth < MAX_RAY_DEPTH) {
		float facingratio = -(Dot(raydir, nhit));
		// change the mix value to the tweak effect
		const float fresneleffect = mix(pow(1 - facingratio, 3), 1.0, 0.1);
		// compute reflection direction (no need to normalize because all vectors
        // are already normalized)
		float3 refldir = raydir - nhit * 2 * Dot(raydir, nhit);
		Normalize(refldir);
		float3 reflection = trace(phit + nhit * bias, refldir, spheres, depth + 1, size);
		float3 refraction = 0;
		// if sphere is also transparent, compute refraction ray
		if (sphere.transparency) {
			float ior = 1.1;
			float eta;
			if (inside) eta = ior;
			else eta = 1 / ior; 
			float cosi = -(Dot(nhit, raydir));
			float k = 1 - eta * eta * (1 - cosi * cosi);
			float3 refrdir = raydir * eta + nhit * (eta * cosi - sqrt(k));
			Normalize(refrdir);
			refraction = trace(phit - nhit * bias, refrdir, spheres, depth + 1, size);
		}
		// the result is a mix of reflection and refraction	
		float3 sphereSurfaceColor = {sphere.surfaceColor[0],
								     sphere.surfaceColor[1],
									 sphere.surfaceColor[2]}; 
	    surfaceColor = (
		    reflection * fresneleffect + 
		    refraction * (1 - fresneleffect) * sphere.transparency) * sphereSurfaceColor;
	}
	else {
		// it's a diffuse object, no need to raytace any further
		for (int i = 0; i < size; ++i) {
			if (spheres[i].emissionColor[0] > 0) {
				// this is a light
				float3 transmission = 1;
				float3 centeri = {spheres[i].center[0],
								  spheres[i].center[1],
								  spheres[i].center[2]};
				float3 lightDirection = centeri - phit;
				Normalize(lightDirection);
				for (uniform int j = 0; j < size; ++j) {
					if (i != j) {
						float t0, t1;
						if (SphereIntersect(spheres[j], phit + nhit * bias, 
													lightDirection, t0, t1)) {
							transmission = 0;
							break;
						}
					}
				}
				float3 sphereSurfaceColor = {sphere.surfaceColor[0],
											 sphere.surfaceColor[1],
											 sphere.surfaceColor[2]}; 
				
				float3 sphereEmissionColor = {
											 spheres[i].emissionColor[0],
											 spheres[i].emissionColor[1],
											 spheres[i].emissionColor[2]}; 
				surfaceColor += sphereSurfaceColor * transmission *
				max(0., Dot(nhit, lightDirection)) * sphereEmissionColor;
			}
		}
	}
    float3 sphereEmissionColor = {
							 sphere.emissionColor[0],
							 sphere.emissionColor[1],
							 sphere.emissionColor[2]}; 
	return surfaceColor + sphereEmissionColor;
}

static void render_tile(uniform int x0, uniform int x1, 
                        uniform int y0, uniform int y1,
                        uniform int width, uniform int height,
                        uniform float * pixels_r, uniform float * pixels_g,
                        uniform float * pixels_b, uniform Sphere * uniform spheres,
                        uniform int sphere_count, uniform float invWidth,
                        uniform float invHeight, uniform float fov,
                        uniform float aspectratio, uniform float angle){
	// Trace rays
	foreach_tiled(y = y0 ... y1, x = x0 ... x1){
		float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio;
		float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
		float3 raydir = {xx, yy, -1};
		Normalize(raydir);
		int pixel_index = width * y + x;
		float3 ray_orig = { 0, 0, 0 };
		float3 pixel_color = trace(ray_orig, raydir, spheres, 0, sphere_count);
		pixels_r[pixel_index] = pixel_color.x;
		pixels_g[pixel_index] = pixel_color.y;
		pixels_b[pixel_index] = pixel_color.z;

	}
}

task void render_tile_task(uniform float * uniform pixels_r,
                           uniform float * uniform pixels_g,
                           uniform float * uniform pixels_b,
                           uniform int width, 
                           uniform int height,
                           uniform Sphere * uniform spheres,
                           uniform int sphere_count,
                           uniform float invWidth,
                           uniform float invHeight,
                           uniform float fov,
                           uniform float aspectratio,
                           uniform float angle){
    uniform int dx = 16, dy = 16;
    uniform int xBuckets = (width + (dx-1)) / dx;
    uniform int x0 = (taskIndex % xBuckets) * dx; 
    uniform int x1 = min(x0 + dx, width);
    uniform int y0 = (taskIndex / xBuckets) * dy; 
    uniform int y1 = min(y0 + dy, height);
    
    render_tile(x0, x1, y0, y1, width, height, pixels_r, pixels_g, pixels_b,
                spheres, sphere_count, invWidth, invHeight, fov, aspectratio,
                angle);
}


export void render_tasks(uniform float * uniform pixels_r, 
                         uniform float * uniform pixels_g,
                         uniform float * uniform pixels_b,
                         uniform int width, 
                         uniform int height,
                         uniform Sphere * uniform spheres,
                         uniform int sphere_count){
	uniform float invWidth = 1 / (float)width, invHeight = 1 / (float)height;
	uniform float fov = 30, aspectratio = width / (float)height;
	uniform float angle = tan(PI * 0.5 * fov / 180.);
    uniform int dx = 16, dy = 16; 
    uniform int xBuckets = (width + (dx-1)) / dx; 
    uniform int yBuckets = (height + (dy-1)) / dy; 
    uniform int nTasks = xBuckets * yBuckets;
    launch[nTasks] render_tile_task(pixels_r, pixels_g, pixels_b, width, height,
                                    spheres, sphere_count, invWidth, invHeight,
                                    fov, aspectratio, angle);
}


