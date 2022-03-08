#ifndef FACE_INCLUDED
#define FACE_INCLUDED
#include <string>
#include <vector>


class Face{
  public:
  
  //actual data container
    std::vector<int> verticies;
    std::vector<int> texture;
    std::vector<double> N;
    std::vector<double> NInv;
    int matIndex;
    
    bool textured = false;
  
    std::vector<double> Na;
    std::vector<double> Nb;
    std::vector<double> Nc;
    
  
	//Explicit constructor
	Face(int a, int b, int c, int matIn) {
        verticies.push_back(a);
        verticies.push_back(b);
        verticies.push_back(c);
        matIndex = matIn;
    }
    //Explicit constructor
	Face(int a, int b, int c, int matIn, int e, int f, int g) {
        verticies.push_back(a);
        verticies.push_back(b);
        verticies.push_back(c);
        matIndex = matIn;
        texture.push_back(e);
        texture.push_back(f);
        texture.push_back(g);
        textured = true;
    }
	//Default constructor
	Face()=delete;
	//Destructor (needs nothing)
	~Face(){
    }

	//Print vertex verticies for debug purposes
    std::string getPrettyString() {
        return "["+std::to_string(verticies[0])+", "+std::to_string(verticies[1])+", "+std::to_string(verticies[2])+"]";
    }

	//Print vertex verticies as in .obj file
    std::string getString() {
        return "f "+std::to_string(verticies[0])+" "+std::to_string(verticies[1])+" "+std::to_string(verticies[2]);
    }
	
  private:
	
    };
#endif 
