#ifndef Sphere_INCLUDED
#define Sphere_INCLUDED
#include <string>
#include <vector>
#include "Texture.h"


class Sphere{
  public:
  
  //actual data container
    double radius;
    std::vector<double> point;
    std::vector<double> ambient;
    std::vector<double> diffuse;
    std::vector<double> specular;
    std::vector<double> attenuation;
    bool textured;
    bool bumped;
    Texture texture;
    Texture bump;
    double phiOffset;
    double thetaOffset;
    int tile = 0;
    double Ni=0;
  
    static constexpr double twoPI = (2*M_PI);
  
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
  
  Sphere(double x, double y, double z, double rad, double ambientR, double ambientG, double ambientB, double diffuseR, double diffuseG, double diffuseB, double specularR, double specularG, double specularB, double attenuationR, double attenuationG, double attenuationB, double NiIn, std::string textureFile = "", double phiIn = 0.0, double thetaIn = 0.0, int tiling = 0) {
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
      
        Ni = NiIn;
        
        tile = tiling;
      
        //cout<<"sphere texture: '"<<textureFile<<"'\n";
        if (textureFile != "") {
          textured=true;
          //cout<<"making textured sphere\n";
          texture = Texture(textureFile, tile);
        }
    
    /*
      if (bumpFile != "") {
          bumped=true;
          //cout<<"making textured sphere\n";
          bump = Texture(bumpFile);
        }*/
    
    
    
      phiOffset = phiIn;
      thetaOffset = thetaIn;
      
    
          
    
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
  
  std::vector<double> pointToTexture(double x, double y, double z) {
    if (!textured)
      return vector<double> {1,1,1};
    
    x = (x - point[0]);
    y = (y - point[1]);
    z = (z - point[2]);
    
    //THETA: -Pi --> Pi
    //PHI: 0 --> Pi
    
    if (abs(x) < 0.0001)
      x = 0.0001;
    if (abs(z) < 0.0001)
      z = 0.0001;
    
    //cout<<"x: "<<x<<", y: "<<y<<", z: "<<z<<'\n';

    
      
    double theta = (atan2(y,x)) + thetaOffset;
    
    while (theta > (M_PI))
      theta = theta - (twoPI);
    
    while (theta < -(M_PI) )
       theta = theta + ( twoPI);
    
    double phi = (acos(z/radius));

    while (phi > (M_PI)) 
      phi = phi - (M_PI);
    while (phi < - 0)
       phi = phi + (M_PI);
    
    theta = (theta + (M_PI)) / ( twoPI);
    phi = ( phi ) / (M_PI);
    

    //double theta = (atan(y/x) + (0)) / (1);
    //double phi = (atan(pow((x*x)+(y*y),0.5)/z) + (0))  / (1);
    
    //cout<<radius<<", "<<theta<<", "<<phi<<'\n';
    //cout<<(theta*360 -180)<<", "<<(phi*180 -90)<<'\n';
   
    return texture.get(theta, (1-phi));
    
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
