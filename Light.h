#ifndef Light_INCLUDED
#define Light_INCLUDED
#include <string>
#include <vector>


class Light{
  public:
  
  //actual data container
    std::vector<double> point;
    std::vector<double> color;
  
	//Explicit constructor
	Light(double x, double y, double z, double w) {
    
        point.push_back(x);
        point.push_back(y);
        point.push_back(z);
        point.push_back(w);
    
        color.push_back(1);
        color.push_back(1);
        color.push_back(1);
    
    
    }
  
  Light(double x, double y, double z, double w, double colorR, double colorG, double colorB) {
    
        point.push_back(x);
        point.push_back(y);
        point.push_back(z);
        point.push_back(w);
    
        color.push_back(colorR);
        color.push_back(colorG);
        color.push_back(colorB);
    
    }
  
  
	//Default constructor
	Light() {
    
        point.push_back(0);
        point.push_back(0);
        point.push_back(0);
        point.push_back(0);
    
        color.push_back(1);
        color.push_back(1);
        color.push_back(1);
    
        
    }
	//Destructor (needs nothing)
	~Light(){
    }
  
	

	//Print Light point for debug purposes
    std::string getPrettyString() {
        return "w = "+std::to_string(point[3])+", ["+std::to_string(point[0])+", "+std::to_string(point[1])+", "+std::to_string(point[2])+"]";
    }

	//Print Light point as in .obj file
    std::string getString() {
        return "Light "+std::to_string(point[0])+" "+std::to_string(point[1])+" "+std::to_string(point[2])+" "+std::to_string(point[3]);
    }
	
  private:
	
    };
#endif 
