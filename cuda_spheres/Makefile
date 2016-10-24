default: spheres


dirs:
	mkdir -p objs

clean:
	rm -rf objs

spheres: dirs 
	nvcc -o objs/$@ spheres.cu

