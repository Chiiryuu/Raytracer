#ifndef VERTEX_INCLUDED
#define VERTEX_INCLUDED
#include <string>


class Vertex{
  public:
  
  //actual data container
    double point[3];
  
	//Explicit constructor
	Vertex(double x, double y, double z) {
        point[0]=x;
        point[1]=y;
        point[2]=z;
    
    }
	//Default constructor
	Vertex() {
        point[0]=0;
        point[1]=0;
        point[2]=0;
    }
	//Destructor (needs nothing)
	~Vertex(){
    }
  
  void update(double x, double y, double z) {
     point[0]=x;
     point[1]=y;
     point[2]=z; 
  }
	

	//Print vertex point for debug purposes
    std::string getPrettyString() {
        return "["+std::to_string(point[0])+", "+std::to_string(point[1])+", "+std::to_string(point[2])+"]";
    }

	//Print vertex point as in .obj file
    std::string getString() {
        return "v "+std::to_string(point[0])+" "+std::to_string(point[1])+" "+std::to_string(point[2]);
    }
	
  private:
	
    };
#endif 
