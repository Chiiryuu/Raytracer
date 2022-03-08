#ifndef OBJECT_INCLUDED
#define OBJECT_INCLUDED
#include "Eigen/Dense"
#include "Vertex.h"
#include "Texture.h"
#include "Face.h"
#include <vector>



class Object{
  public:
	Object(std::string fileName);
    Object(std::string fileName, int tiling, double smoothingIn);
	Object();
	~Object();
  
  std::string name;
  std::string file;
    
  int tile=0;
  
  double sumOriginalToTransform=0;
  
  double sumOriginalToOriginal=0;
  
  Eigen::MatrixXd vertMatrix;
  
  Eigen::MatrixXd transformMatrix;
  Eigen::MatrixXd inverseTransformMatrix;
  
  std::string copyrightInfo;
  
  double smoothing=0; 
  bool smoothingGroup=false;
    
  std::vector<std::string> textureNames;  
  std::vector<Texture> textureFiles; 
  std::vector<int> textureIndex;
  
  std::vector<Vertex> verticies;
  std::vector<double> u;
  std::vector<double> v; 
  std::vector<Face> faces;
    
  std::vector<std::string> materialNames;
  std::vector<double> specularExponents;
//  std::vector<double[3]> ambient; //Ka
  std::vector<double> aR;
  std::vector<double> aG;
  std::vector<double> aB;
    
//  std::vector<double[3]> diffuse; //Kd
  std::vector<double> dR;
  std::vector<double> dG;
  std::vector<double> dB;
    
//  std::vector<double[3]> specular; //Ks
  std::vector<double> sR;
  std::vector<double> sG;
  std::vector<double> sB;
    
//  std::vector<double[3]> specular; //Tr
  std::vector<double> TrR;
  std::vector<double> TrG;
  std::vector<double> TrB;
    
  std::vector<double> Ni;
    
  std::vector<int> illuminationModels;
  
    
  
  void loadFromFile();

  void makeMatrix();
    
  void makeMaterial(std::string matPath);
  
  void updateVerts();
  
  
  void updateFaces();
  
  void makeSmoothing();
  
  void addVertex(double x, double y, double z);
  void addFace(int a, int b, int c, int id);
  void addFace(int a, int b, int c, int id, int e, int f, int g);
  //void addFace(int a, int normA, int b, int normB, int c, int normC);
  
  std::string printTransformMatrix();
  
  std::string printInverseTransformMatrix();
  
  std::string printSumOriginalToTransform();
  
  std::string printSumOriginalToOriginal();
  
  std::string printFile();
  
  std::string printTransformFile();
  
  void setName(std::string newName) {
    name = newName;
  }
  
  std::string getFileName() {
    return file;
  }
  
  std::string getName() {
    return name;
  }

	
  private:
	
    };
#endif 