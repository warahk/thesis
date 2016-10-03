default: spheres


dirs:
	mkdir -p objs

clean:
	rm -rf objs/sphere*

spheres: dirs objs/spheres.o objs/spheres_ispc.o
	g++ -o objs/$@ objs/spheres.o objs/spheres_ispc.o

objs/spheres.o: objs/spheres_ispc.h
	g++ -o $@ -c spheres.cpp -I./objs

objs/spheres_ispc.h objs/spheres_ispc.o: spheres.ispc
	ispc -o objs/spheres_ispc.o -h objs/spheres_ispc.h spheres.ispc --target=avx2-i32x16 
