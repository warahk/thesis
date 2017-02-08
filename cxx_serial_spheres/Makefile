default: spheres


dirs:
	mkdir -p objs

clean:
	rm -rf objs

spheres: dirs 
	g++ -o objs/$@ --std=c++11 raytracer.cpp

