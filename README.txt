Brandon Percin
eid percinbp
CS410 Assignment 5

** Attributions for 3D Models and Textures Downloaded from the Internet **
	*TEXTURES*
		metal.ppm - https://www.123rf.com/photo_33951250_scratched-and-grunge-metal-sheet-texture-and-seamless-background.html
		earth.ppm, 
		moon.ppm - https://www.solarsystemscope.com/textures/
		stars.ppm - http://www.everystockphoto.com/photo.php?imageId=9268253
		
	*MODELS*
		mothership.obj - https://sketchfab.com/3d-models/mothership-896f36009f2d41878d01a855ae09eb12
		ship.obj - https://sketchfab.com/3d-models/gemini-8fdac465e8804fedaa1249f81484abac#download

New additions this version:
* Texture Mapping on Objects
* Texture tiling on Objects
* Texture mapping on spheres
* Texture tiling on spheres

NOTE: 	multithreaded rendering was implemented quickly, and because of this, extremely large images
	(typically > 512x512) can use up too much memory, resulting in a segfault, depending on the
	amount of objects in the scene and amount of textures in memory. To remedy this, if you 
	need a high resolution rendering, run the program in a single thread by calling the
	program with the flag -f.


Program arrives with a version of the Eigen library. Of course, nothing in the Eigen directory
was made by me. The program consists of 9 c++ files:
* Vertex.h 	(essentially an instantiable vector of doubles of length 3)
* Face.h 	(an instantiable vector of verticies and of face normals, both length 3)
* Object.h	(all header information for dealing with objects)
* Object.cc 	(actual code specific to each instance of an object as created from driver file)
* Camera.h 	(holds data relating to the camera)
* Light.h 	(holds data relating to lights)
* Ray.h		(holds all data required for a ray, including hit point)
* Sphere.h 	(holds information relating to a sphere object)
* Texture.h	(holds a texture file that can be accessed with u,v coordinates)
* raytrace.cc	(not instantiable; has a main method. Performs all calculations and writes to file.)

The default make target will produce a program 'raytracer', so simply use 'make' to compile.

After this, program runs as specified in the assignment, i.e.
	$./raytracer driver00.txt driver00.ppm
will run the application, taking in driver00.txt in the local directory as input.

Optionally, run single-threaded rendering by running with a flag, i.e.
	$./raytracer driver00.txt driver00.ppm -f
Only do this if you run into memory issues!

Output files are placed in current directory under the assigned name.

In addition, when run, this assignment will keep track of how close it is to completing 
its task, and will print a percentage to console. This was implemented to make sure the 
program isn't hung up on something.

Description of Drivers

driver00.txt: Composed of only 2 spheres, a box, and a light model, driver00.txt produces a render of the Moon casting a shadow
		over Earth (with both the Earth and Moon not utilizing texture tiling) with a giant box containing the entire scene,
		having a texture of stars tiling 200 times in both dimensions in order to make a starry background.
		This scene is a simple example of both texture mapping and tiling.
		

driver01.txt: Using Free-to-use models (attributed above) from sketchfab, driver01.txt is a camera looking directly at a giant metal sphere.
		This sphere, representing an alien spacecraft, reflects an image of Earth and a space fleet ready to engage in battle.
		The moon can be seen an acurate distance away in the reflection.
	      
