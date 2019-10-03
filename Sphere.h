#ifndef Sphere_INCLUDED
#define Sphere_INCLUDED
#include <string>
#include <vector>


class Sphere{
  public:
  
  //actual data container
    double radius;
    std::vector<double> point;
    std::vector<double> ambient;
    std::vector<double> diffuse;
    std::vector<double> specular;
    std::vector<double> attenuation;
  
	//Explicit constructor
	Sphere(double x, double y, double z, double rad) {
        radius = rad;
    
        point.push_back(x);
        point.push_back(y);
        point.push_back(z);
    
        ambient.push_back(0);
        ambient.push_back(0);
        ambient.push_back(0);
    
        diffuse.push_back(0);
        diffuse.push_back(0);
        diffuse.push_back(0);
    
        specular.push_back(0);
        specular.push_back(0);
        specular.push_back(0);
    
        attenuation.push_back(0);
        attenuation.push_back(0);
        attenuation.push_back(0);
    
    }
  
  Sphere(double x, double y, double z, double rad, double ambientR, double ambientG, double ambientB, double diffuseR, double diffuseG, double diffuseB, double specularR, double specularG, double specularB, double attenuationR, double attenuationG, double attenuationB) {
        radius = rad;
    
        point.push_back(x);
        point.push_back(y);
        point.push_back(z);
    
        ambient.push_back(ambientR);
        ambient.push_back(ambientG);
        ambient.push_back(ambientB);
    
        diffuse.push_back(diffuseR);
        diffuse.push_back(diffuseG);
        diffuse.push_back(diffuseB);
    
        specular.push_back(specularR);
        specular.push_back(specularG);
        specular.push_back(specularB);
    
        attenuation.push_back(attenuationR);
        attenuation.push_back(attenuationG);
        attenuation.push_back(attenuationB);
    
    }
  
  
	//Default constructor
	Sphere() {
        radius=1;
    
        point.push_back(0);
        point.push_back(0);
        point.push_back(0);
    
        ambient.push_back(0);
        ambient.push_back(0);
        ambient.push_back(0);
    
        diffuse.push_back(0);
        diffuse.push_back(0);
        diffuse.push_back(0);
    
        specular.push_back(0);
        specular.push_back(0);
        specular.push_back(0);
    
        attenuation.push_back(0);
        attenuation.push_back(0);
        attenuation.push_back(0);
    
        
    }
	//Destructor (needs nothing)
	~Sphere(){
    }
  
  void update(double x, double y, double z, double rad) {
     radius = rad;
     point[0]=x;
     point[1]=y;
     point[2]=z; 
  }
	

	//Print Sphere point for debug purposes
    std::string getPrettyString() {
        return "r = "+std::to_string(radius)+", ["+std::to_string(point[0])+", "+std::to_string(point[1])+", "+std::to_string(point[2])+"]";
    }

	//Print Sphere point as in .obj file
    std::string getString() {
        return "sphere "+std::to_string(radius)+" "+std::to_string(point[0])+" "+std::to_string(point[1])+" "+std::to_string(point[2]);
    }
	
  private:
	
    };
#endif 
