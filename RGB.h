#ifndef RGB_INCLUDED
#define RGB_INCLUDED
#include <string>
#include <vector>


class RGB{
  public:
  
  //actual data container
    std::vector<double> point;
  
	//Explicit constructor
	RGB(double x, double y, double z) {
        point.push_back(x);
        point.push_back(y);
        point.push_back(z);
    
    }
	//Default constructor
	RGB() {
        point.push_back(0);
        point.push_back(0);
        point.push_back(0);
    }
	//Destructor (needs nothing)
	~RGB(){
    }
  
  void update(double x, double y, double z) {
     point.clear();
     point.push_back(x);
     point.push_back(y);
     point.push_back(z);
  }
	

	//Print RGB point for debug purposes
    std::string getPrettyString() {
        return "["+std::to_string(point[0])+", "+std::to_string(point[1])+", "+std::to_string(point[2])+"]";
    }

	//Print RGB point as in .obj file
    std::string getString() {
        return "v "+std::to_string(point[0])+" "+std::to_string(point[1])+" "+std::to_string(point[2]);
    }
	
  private:
	
    };
#endif 
