#ifndef Camera_INCLUDED
#define Camera_INCLUDED
#include <string>
#include <vector>


class Camera{
  public:
  
  //actual data container
    
    
  
    std::vector<double> eye;
    std::vector<double> look;
    std::vector<double> up;
    std::vector<double> w;
    std::vector<double> u;
    std::vector<double> v;
  
    double d;
  
    std::vector<double> bounds;
    std::vector<int> res;
  
 
  
  Camera(double x, double y, double z, double lookX, double lookY, double lookZ, double upX, double upY, double upZ, double dIn, 
         double bounds0, double bounds1, double bounds2, double bounds3, int resX, int resY,
        double wX, double wY, double wZ, double uX, double uY, double uZ, double vX, double vY, double vZ ) {
    
    eye.push_back(x);
    eye.push_back(y);
    eye.push_back(z);
    
    look.push_back(lookX);
    look.push_back(lookY);
    look.push_back(lookZ);
    
    up.push_back(upX);
    up.push_back(upY);
    up.push_back(upZ);
    
    d = dIn;
    
    bounds.push_back(bounds0);
    bounds.push_back(bounds1);
    bounds.push_back(bounds2);
    bounds.push_back(bounds3);
    
    res.push_back(resX);
    res.push_back(resY);
    
    w.push_back(wX);
    w.push_back(wY);
    w.push_back(wZ);
    
    u.push_back(uX);
    u.push_back(uY);
    u.push_back(uZ);
    
    v.push_back(vX);
    v.push_back(vY);
    v.push_back(vZ);
    
    
    }
  
  
	//Default constructor
	Camera() {
    
    
    eye.push_back(0);
    eye.push_back(0);
    eye.push_back(0);
    
    look.push_back(1);
    look.push_back(0);
    look.push_back(0);
    
    up.push_back(0);
    up.push_back(0);
    up.push_back(1);
    
    d = 1;
    
    bounds.push_back(-1);
    bounds.push_back(1);
    bounds.push_back(-1);
    bounds.push_back(1);
    
    res.push_back(10);
    res.push_back(10);
    
    
    }
	//Destructor (needs nothing)
	~Camera(){
    }
  
	

	//Print Camera eye for debug purposes
    std::string getPrettyString() {
    return "Camera: d = "+std::to_string(d)+"\neye = ["+std::to_string(eye[0])+", "+std::to_string(eye[1])+", "+std::to_string(eye[2])+"]"
      +"\nlook = ["+std::to_string(look[0])+", "+std::to_string(look[1])+", "+std::to_string(look[2])+"]"
        +"\nup = ["+std::to_string(up[0])+", "+std::to_string(up[1])+", "+std::to_string(up[2])+"]"
          +"\nw = ["+std::to_string(w[0])+", "+std::to_string(w[1])+", "+std::to_string(w[2])+"]"
            +"\nu = ["+std::to_string(u[0])+", "+std::to_string(u[1])+", "+std::to_string(u[2])+"]"
             +"\nv = ["+std::to_string(v[0])+", "+std::to_string(v[1])+", "+std::to_string(v[2])+"]";
    }

	//Print Camera eye as in .obj file
    std::string getString() {
    return "Camera "+std::to_string(eye[0])+" "+std::to_string(eye[1])+" "+std::to_string(eye[2]);
    }
	
  private:
	
    };
#endif 
