#ifndef VERTEX_INCLUDED
#define VERTEX_INCLUDED
#include <string>
#include <vector>


class Vertex{
  public:
  
  //actual data container
    std::vector<double> point;
  
    std::vector<int> faces;
  
	//Explicit constructor
	Vertex(double x, double y, double z) {
        point.push_back(x);
        point.push_back(y);
        point.push_back(z);
    
    }
	//Default constructor
	Vertex() {
        point.push_back(0);
        point.push_back(0);
        point.push_back(0);
    }
	//Destructor (needs nothing)
	~Vertex(){
    }
  
  void update(double x, double y, double z) {
     point.clear();
     point.push_back(x);
     point.push_back(y);
     point.push_back(z);
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
