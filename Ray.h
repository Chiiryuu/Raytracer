#ifndef Ray_INCLUDED
#define Ray_INCLUDED
#include <string>
#include <vector>


class Ray{
  public:
  
  //actual data container
    std::vector<double> origin;
    double size = 0;
    std::vector<double> direction;
    std::vector<double> hitPoint;
    bool hit = false;
  
	//Explicit constructor
	Ray(double orgX, double orgY, double orgZ, double dirX, double dirY, double dirZ) {
        origin.push_back(orgX);
        origin.push_back(orgY);
        origin.push_back(orgZ);
    
        
    
        direction.push_back(dirX);
        direction.push_back(dirY);
        direction.push_back(dirZ);
    
    }
	//Default constructor
	Ray() {
        origin.push_back(0);
        origin.push_back(0);
        origin.push_back(0);
    
        
    
        direction.push_back(1);
        direction.push_back(0);
        direction.push_back(0);
    }
	//Destructor (needs nothing)
	~Ray(){
    }
  /*
  void update(double orgX, double orgY, double orgZ, double dirX, double dirY, double dirZ) {
        size = 0;
        origin[0] = orgX;
        origin[1] = orgY;
        origin[2] = orgZ;
    
        
    
        direction[0]=dirX;
        direction[1]=dirY;
        direction[2]=dirZ;
  }*/
  
  void setHitPosition(double x, double y, double z) {
    hitPoint.clear();
    hitPoint.push_back(x);
    hitPoint.push_back(y);
    hitPoint.push_back(z);
  }
  
  double tx(double t) {
    return origin[0] + t * direction[0];
  }
  
  double ty(double t) {
    return origin[1] + t * direction[1];
  }
  
  double tz(double t) {
    return origin[2] + t * direction[2];
  }
  
  std::vector<double> rt(double t) {
    return std::vector<double> {tx(t), ty(t), tz(t)};
  }
  
  bool checkHit() {
   return hit;
  }
  
  void resetHit() {
   hit = false;
  }
  
  void setHit() {
   hit = true;
  }
	

	//Print Ray direction for debug purposes
    std::string getPrettyString() {
        return "Origin: ["+std::to_string(origin[0])+", "+std::to_string(origin[1])+", "+std::to_string(origin[2])+"]\n"+"Direction: ["+std::to_string(direction[0])+", "+std::to_string(direction[1])+", "+std::to_string(direction[2])+"]";
    }

	//Print Ray direction as in .obj file
    std::string getString() {
        return "["+std::to_string(direction[0])+" "+std::to_string(direction[1])+" "+std::to_string(direction[2])+"]";
    }
	
  private:
	
    };
#endif 
