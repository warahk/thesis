default: spheres


dirs:
	mkdir -p objs

clean:
	rm -rf objs/sphere*

spheres: dirs objs/spheres.o objs/spheres_ispc.o objs/tasksys.o
	g++ -o objs/$@ objs/spheres.o objs/spheres_ispc.o objs/tasksys.o -pthread

objs/spheres.o: objs/spheres_ispc.h
	g++ -o $@ -c spheres.cpp -I./objs

objs/tasksys.o:
	g++ -c tasksys.cpp -o $@ -DISPC_USE_PTHREADS

objs/spheres_ispc.h objs/spheres_ispc.o: spheres.ispc
	ispc -o objs/spheres_ispc.o -h objs/spheres_ispc.h spheres.ispc --target=avx2-i32x16 
