Brandon Percin
eid percinbp
CS410 Assignment 2

As this assignment builds on the previous, some files from the previous assignment have been 
included even though not explicitly used in this one. This is because both pieces will
likely be needed for the next assignment.

Program arrives with a version of the Eigen library. Of course, nothing in the Eigen directory
was made by me. The program consists of 5 c++ files:
* Vertex.h 	(essentially an instantiable vector of doubles of length 3)
* Face.h 	(an instantiable vector of verticies and of face normals, both length 3)
* Object.h	(all header information for dealing with objects)
* Object.cc 	(actual code specific to each instance of an object as created from driver file)
* Camera.h 	(holds data relating to the camera)
* Light.h 	(holds data relating to lights)
* Ray.h		(holds all data required for a ray, including hit point)
* Sphere.h 	(holds information relating to a sphere object)
* raytrace.cc	(not instantiable; has a main method. Performs all calculations and writes to file.)

The default make target will produce a program 'raytracer', so simply use 'make' to compile.

After this, program runs as specified in the assignment, i.e.
	$./raytracer driver00.txt driver00.ppm
will run the application, taking in driver00.txt in the local directory as input.

Output files are placed in current directory under the assigned name.

