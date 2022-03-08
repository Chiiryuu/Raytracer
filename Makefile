program: Vertex.h Face.h Object.h Object.cc Ray.h Sphere.h Light.h Camera.h
	g++ raytrace.cc Object.cc Ray.h Sphere.h Light.h Camera.h -o raytracer -O3 -fopenmp

tar:
	tar -cv $(MAKEFILE_LIST) *.cc $(wildcard *.h) README.txt Eigen > Brandon-Percin-P5.tar