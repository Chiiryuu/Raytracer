#ifndef OBJECT_INCLUDED
#define OBJECT_INCLUDED
#include "Eigen/Dense"
#include "Vertex.h"
#include "Face.h"
#include <vector>



class Object{
  public:
	Object(std::string fileName);
	Object();
	~Object();
  
  std::string name;
  std::string file;
  
  double sumOriginalToTransform=0;
  
  double sumOriginalToOriginal=0;
  
  Eigen::MatrixXd vertMatrix;
  
  Eigen::MatrixXd transformMatrix;
  Eigen::MatrixXd inverseTransformMatrix;
  
  std::string copyrightInfo;
  
  bool smoothingGroup=false;
  
  std::vector<Vertex> verticies;
  std::vector<Face> faces;
  
  void loadFromFile();

  void makeMatrix();
  
  void updateVerts();
  
  void addVertex(double x, double y, double z);
  void addFace(int a, int b, int c);
  void addFace(int a, int normA, int b, int normB, int c, int normC);
  
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