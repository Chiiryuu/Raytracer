#ifndef FACE_INCLUDED
#define FACE_INCLUDED
#include <string>


class Face{
  public:
  
  //actual data container
    int verticies[3];
    int normals[3];
  
	//Explicit constructor
	Face(int a, int b, int c) {
        verticies[0]=a;
        verticies[1]=b;
        verticies[2]=c;
        normals[0]=0;
        normals[1]=0;
        normals[2]=0;
    }
  //Explicit constructor With Normals
	Face(int a, int normA, int b, int normB, int c, int normC) {
        verticies[0]=a;
        verticies[1]=b;
        verticies[2]=c;
        normals[0]=normA;
        normals[1]=normB;
        normals[2]=normC;
    }
	//Default constructor
	Face()=delete;
	//Destructor (needs nothing)
	~Face(){
    }

	//Print vertex verticies for debug purposes
    std::string getPrettyString() {
        return "["+std::to_string(verticies[0])+"//"+std::to_string(normals[0]) +", "+std::to_string(verticies[1])+"//"+std::to_string(normals[1])+", "+std::to_string(verticies[2])+"//"+std::to_string(normals[2])+"]";
    }

	//Print vertex verticies as in .obj file
    std::string getString() {
        return "f "+std::to_string(verticies[0])+"//"+std::to_string(normals[0]) +" "+std::to_string(verticies[1])+"//"+std::to_string(normals[1])+" "+std::to_string(verticies[2])+"//"+std::to_string(normals[2]);
    }
	
  private:
	
    };
#endif 
