#include "Object.h"
#include <iostream>
#include <fstream>
#include <string>
#include <regex> 
#include <cmath>
#include <sstream>
#include <iomanip>

using namespace std;


Object::Object(string fileName, int tiling, double smoothingIn): file(fileName){
    name = file;
    //Remove .obj
    name.erase(name.length()-4,5);
    
    tile = tiling;
    
    materialNames.push_back("None");
    specularExponents.push_back(16.0);
    aR.push_back(0.2);
    aG.push_back(0.2);
    aB.push_back(0.2);
    dR.push_back(0.8);
    dG.push_back(0.8);
    dB.push_back(0.8);
    sR.push_back(0.2);
    sG.push_back(0.2);
    sB.push_back(0.2);
    illuminationModels.push_back(2);
    Ni.push_back(1.0);
    TrR.push_back(1.0);
    TrG.push_back(1.0);
    TrB.push_back(1.0);
    textureFiles.push_back(Texture(""));
    textureNames.push_back("");
    textureIndex.push_back(0);
    
    
    
    loadFromFile();
    
    smoothing = smoothingIn;
    if (smoothing > 0.0001)
        smoothingGroup=true;
    
    makeMatrix();
    
    //cout<<"Texture 0: "<<textureFiles[0]<<'\n';
}


Object::Object(string fileName): file(fileName){
    name = file;
    //Remove .obj
    name.erase (name.length()-4,5);
    
    materialNames.push_back("None");
    specularExponents.push_back(16.0);
    aR.push_back(0.2);
    aG.push_back(0.2);
    aB.push_back(0.2);
    dR.push_back(0.8);
    dG.push_back(0.8);
    dB.push_back(0.8);
    sR.push_back(0.2);
    sG.push_back(0.2);
    sB.push_back(0.2);
    illuminationModels.push_back(2);
    Ni.push_back(1.0);
    TrR.push_back(1.0);
    TrG.push_back(1.0);
    TrB.push_back(1.0);
    textureFiles.push_back(Texture(""));
    textureNames.push_back("");
    textureIndex.push_back(0);
    
    
    loadFromFile();
    makeMatrix();
}

Object::Object() : file("None"s){
    name = file;
}


Object::~Object() {
    //delete self;
}


double dot(vector<double>& a, vector<double>& b) {
    return (a[0]*b[0]) + (a[1]*b[1]) + (a[2]*b[2]);
}

void Object::makeSmoothing() {
    double cutoff = cos(smoothing);
    for (int i=0;i<faces.size();i++) {
        Face& face = faces[i];
        vector<double>& N = face.N; 
        //cout<<"Face "<<i<<" neighbors: ";
        for (int j=0;j<face.verticies.size();j++) {
            Vertex& vertex = verticies[face.verticies[j]];
            //verticies[j].faces.push_back(i);
            int smoothingLength = 1;
            vector<double> newNormal = N;
            for (int k=0;k<vertex.faces.size();k++) {
                if (vertex.faces[k] != i) {                   
                    Face& face2 = faces[vertex.faces[k]];
                    vector<double>& N2 = face2.N; 
                    double cosTheta = dot(N, N2);
                    //cout<<"Expect "<<cutoff<<" <= "<<cosTheta<<'\n';
                    if (cosTheta >= cutoff) {
                        smoothingLength++;
                        newNormal[0] = newNormal[0] + N2[0];
                        newNormal[1] = newNormal[1] + N2[1];
                        newNormal[2] = newNormal[2] + N2[2];
                    }
                //cout<< vertex.faces[k] <<' ';
                }
            }
            //cout<<"Raw: [ "<<newNormal[0]<<' '<<newNormal[1]<<' '<<newNormal[2]<<", "<<smoothingLength<<" ]\n";
                
            
                newNormal[0] = newNormal[0] / smoothingLength;
                newNormal[1] = newNormal[1] / smoothingLength;
                newNormal[2] = newNormal[2] / smoothingLength;
            
            
            
                
                double mag = 0;
                for (double d : newNormal)
                    mag = mag + (d*d);
                mag = pow(mag,0.5);
            
                //cout<<"Raw: [ "<<newNormal[0]<<' '<<newNormal[1]<<' '<<newNormal[2]<<", "<<mag<<" ]\n";
        
        
                for (int k=0; k<3;k++) {
                    newNormal[k] = newNormal[k] / mag;
                }
                
            
                //cout<<"Norm [ "<<newNormal[0]<<' '<<newNormal[1]<<' '<<newNormal[2]<<", "<<1<<" ]\n";
                
                
                if (j==0)
                    face.Na=newNormal;
                else if (j==1)
                    face.Nb = newNormal;
                else
                    face.Nc = newNormal;
           
        }
         //cout<<'\n';
        
    }
}

void Object::updateVerts() {
    for (int i=0;i<verticies.size();i++) {
        verticies[i].update(vertMatrix(0,i),vertMatrix(1,i),vertMatrix(2,i));
    }
    vector<double> AvBv =vector<double>{0,0,0};
    vector<double> AvCv =vector<double>{0,0,0};
    vector<double> newN =vector<double>{0,0,0};
    vector<double> newNInv =vector<double>{0,0,0};
    for (int i=0;i<faces.size();i++) {
        Face& face = faces[i];
        
       // for (int j=0;j<face.verticies.size();j++) {
         //   verticies[face.verticies[j]].faces.push_back(i);
        //}
        
        vector<double> a = verticies[face.verticies[0]].point;
        vector<double> b = verticies[face.verticies[1]].point;
        vector<double> c = verticies[face.verticies[2]].point;
        
        //cout<<a[0]<<' '<<a[1]<<' '<<a[2]<<'\n';
        //cout<<b[0]<<' '<<b[1]<<' '<<b[2]<<'\n';
        //cout<<c[0]<<' '<<c[1]<<' '<<c[2]<<'\n';
        
        AvBv[0]=b[0]-a[0];
        AvBv[1]=b[1]-a[1];
        AvBv[2]=b[2]-a[2];
            
        AvCv[0]=c[0]-a[0];
        AvCv[1]=c[1]-a[1];
        AvCv[2]=c[2]-a[2];
        
        //cout<<AvBv[0]<<' '<<AvBv[1]<<' '<<AvBv[2]<<'\n';
        //cout<<AvCv[0]<<' '<<AvCv[1]<<' '<<AvCv[2]<<'\n';
        
        newN[0] = AvBv[1]*AvCv[2]  -   AvBv[2]*AvCv[1];
        newN[1] = AvBv[2]*AvCv[0]  -   AvBv[0]*AvCv[2];
        newN[2] = AvBv[0]*AvCv[1]  -   AvBv[1]*AvCv[0];
        
        double mag = 0;
        for (double d : newN)
            mag = mag + (d*d);
        mag = pow(mag,0.5);
        
        
        for (int j=0; j<3;j++) {
            newN[j] = newN[j] / mag;
        }
        
        newNInv[0] = newN[0] * -1;
        newNInv[1] = newN[1] * -1;
        newNInv[2] = newN[2] * -1;
        
        //cout<<newN[0]<<' '<<newN[1]<<' '<<newN[2]<<'\n'<<'\n';
            
        face.N = newN;
        face.NInv = newNInv;
        
    }
    if (smoothingGroup)
            makeSmoothing();
}

void Object::makeMatrix() {
    //Iterate over verticies, filling in columns of the matrix.
    vertMatrix.resize(4,verticies.size());
    int col = 0;
    for (Vertex v : verticies) {
        vertMatrix(0,col)=v.point[0];
        vertMatrix(1,col)=v.point[1];
        vertMatrix(2,col)=v.point[2];
        vertMatrix(3,col)=1;
        col++;
    }
    //cout<<vertMatrix<<'\n';
}

void Object::loadFromFile() {
    file = "models/"+file;
    //Open file stream
    ifstream rawFile(file);
    string line;
    //Check if file was opened
    if (!rawFile)
        cout << "File '"<<file<<"' not found\n";
    
    
    int matIndex = 0;
    string matName = "None";
    
    //Iterate through the file
    while (getline(rawFile,line)) {
        //Open string stream
        istringstream lineStream(line);
        //Skip up to 3 characters, until space
        lineStream.ignore(3,' ');
        if (line.at(0)=='v' && line.at(1)=='n'  && line.at(2)==' ');
            //Pass, ignore vn
        else if (line.at(0)=='v' && line.at(1)==' ') {
            //Create vertex from line in file
            double coords[3];
            lineStream>>coords[0];
            lineStream>>coords[1];
            lineStream>>coords[2];
            addVertex(coords[0],coords[1],coords[2]);
        }
        else if (line.at(0)=='v' && line.at(1)=='t'  && line.at(2)==' ') {
            //Create vertex from line in file
            double coords[2];
            lineStream>>coords[0];
            lineStream>>coords[1];
            u.push_back(coords[0]);
            v.push_back(coords[1]);
            //addTextureVertex(coords[0],coords[1],coords[2]);
        }
        else if (line.at(0)=='s' && line.at(1)==' ') {
            //Set smoothing group from s line
            if (line.at(line.length()-1) == 'n')
                ;
            //Smoothing now handled in driver file.
        }
        else if (line.at(0)=='f' && line.at(1)==' ') {
            //Get face verts and face normals from line
            int verts[3];
            int texture[3] = {-1,-1,-1};
            
            string word;
            char tempChar;
            
            bool textured = false;
            
            
            lineStream>>word;
            //cout<<word<<'\n';
            istringstream wordStream1(word);
            wordStream1>>verts[0];
            wordStream1>>tempChar; 
            wordStream1>>texture[0];
            
            lineStream>>word;
            istringstream wordStream2(word);
            wordStream2>>verts[1];
            wordStream2>>tempChar; 
            wordStream2>>texture[1];
            
            lineStream>>word;
            istringstream wordStream3(word);
            wordStream3>>verts[2];
            wordStream3>>tempChar; 
            wordStream3>>texture[2];
            
            if (texture[0] != -1 && texture[1] != -1 && texture[2] != -1 ) {
              textured=true;
              //cout<<texture[0]<<' '<<texture[1]<<' '<<texture[2]<<'\n'; 
            }
            
            /*
            lineStream>>verts[0];
            lineStream.ignore(2);
            lineStream>>norms[0];
            
            lineStream>>verts[1];
            lineStream.ignore(2);
            lineStream>>norms[1];
            
            lineStream>>verts[2];
            lineStream.ignore(2);
            lineStream>>norms[2];*/ 
            
            if (!textured)
                addFace(verts[0],verts[1],verts[2], matIndex);
            else {
                addFace(verts[0],verts[1],verts[2], matIndex, texture[0], texture[1], texture[2]);
                //cout<<"Face is textured.\n";
            }
        }
        
        else if (line.at(0)=='m' && line.at(1)=='t'  && line.at(2)=='l' && line.at(3)=='l' && line.at(4)=='i'  && line.at(5)=='b') {
            string word = line.substr(7);
            makeMaterial(word);
            //cout<<"Mat file: "<<word<<'\n';
        }
        else if (line.at(0)=='u' && line.at(1)=='s'  && line.at(2)=='e' && line.at(3)=='m' && line.at(4)=='t'  && line.at(5)=='l') {
            matName = line.substr(7);
            for (int i=1;i<materialNames.size();i++) {
                //cout<<materialNames[i]<<" == "<<matName<<" ? "<<(materialNames[i] == matName) <<'\n';
                if (materialNames[i] == matName) {
                    matIndex = i;
                    break;
                }
                
            }
            //cout<<"Material changed to "     <<  matIndex<<'\n';
            //cout<<"Mat file: "<<word<<'\n';
        }
        else
            //Assume text is copyright info
            copyrightInfo = copyrightInfo+line+'\n';
    }
}

/*
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
    
  std::vector<int> illuminationModels;*/

void Object::makeMaterial(string matPath) {
    matPath = "materials/"+matPath;
    ifstream matFile(matPath);
    if (!matFile) {
        cerr<<"Cannot access material file '"s+matPath+"'\n";
    }
    else {
        string line;
        string word;
        int tempInt;
        double tempDouble;
        double TrRVal = 0.0;
        double TrGVal = 0.0;
        double TrBVal = 0.0;
        double NiVal = 0.0;
        string file = "";
        while (getline(matFile,line)) {
            istringstream lineStream(line);
            if (lineStream>>word) {
                if (word=="newmtl") {
                    lineStream>>word;
                    materialNames.push_back(word);
                    int num = materialNames.size() - 1;
                    //cout<<"Created material "<<word<<", i="<<num<<'\n';
                }
                else if (word=="Ni") {
                    lineStream>>NiVal;
                    Ni.push_back(NiVal);
                }
                else if (word=="Tr") {
                    lineStream>>TrRVal; 
                    lineStream>>TrGVal; 
                    lineStream>>TrBVal; 
                    
                    TrR.push_back(TrRVal);
                    TrG.push_back(TrGVal);
                    TrB.push_back(TrBVal);
                }
                else if (word=="Ns") {
                    lineStream>>tempDouble;
                    specularExponents.push_back(tempDouble);
                }
                else if (word=="Ka") {
                    lineStream>>tempDouble;
                    aR.push_back(tempDouble);
                    lineStream>>tempDouble;
                    aG.push_back(tempDouble);
                    lineStream>>tempDouble;
                    aB.push_back(tempDouble);
                }
                else if (word=="Kd") {
                    lineStream>>tempDouble;
                    dR.push_back(tempDouble);
                    lineStream>>tempDouble;
                    dG.push_back(tempDouble);
                    lineStream>>tempDouble;
                    dB.push_back(tempDouble);
                }
                else if (word=="Ks") {
                    lineStream>>tempDouble;
                    sR.push_back(tempDouble);
                    lineStream>>tempDouble;
                    sG.push_back(tempDouble);
                    lineStream>>tempDouble;
                    sB.push_back(tempDouble);
                }
                else if (word=="illum") {
                    lineStream>>tempInt;
                    illuminationModels.push_back(tempInt);
                }
                else if (word=="map_Kd") {
                    lineStream>>file;
                    int index=-1;
                    
                    for (int i=0; i<textureNames.size(); i++) {
                        //cout<<textureNames[i]<<" == "<<file<<" ? -> "<< (textureNames[i] == file)<<'\n';
                        if (textureNames[i] == file)
                            index = i;
                    }
                    
                    if (index == -1) {
                        index = textureNames.size();
                        cout<<"Loading texture file: "<<file<<'\n';
                        textureFiles.push_back(Texture(file, tile));
                        textureNames.push_back(file);
                    }
                   else {
                        cout<<"Re-using texture file: "<<file<<'\n';
                   }
                    
                    textureIndex.push_back(index);
                    
                    //cout<<"File name: "<<file<<'\n';
                    
                    //cout<<file<<'\n';
                }
                
            }
        }

    }
    
}

string Object::printTransformMatrix() {
    //Print transform matrix legibly
    stringstream textStream;
    textStream << fixed << setprecision(3) 
               << transformMatrix(0,0) << ' ' << transformMatrix(0,1) << ' ' << transformMatrix(0,2) << ' ' << transformMatrix(0,3) << ' ' << '\n'
               << transformMatrix(1,0) << ' ' << transformMatrix(1,1) << ' ' << transformMatrix(1,2) << ' ' << transformMatrix(1,3) << ' ' << '\n'
               << transformMatrix(2,0) << ' ' << transformMatrix(2,1) << ' ' << transformMatrix(2,2) << ' ' << transformMatrix(2,3) << ' ' << '\n'
               << transformMatrix(3,0) << ' ' << transformMatrix(3,1) << ' ' << transformMatrix(3,2) << ' ' << transformMatrix(3,3) << ' ' << '\n';
        return textStream.str();
}

string Object::printInverseTransformMatrix() {
    //Print  inverse transform matrix legibly
    stringstream textStream;
    textStream << fixed << setprecision(3) 
               << inverseTransformMatrix(0,0) << ' ' << inverseTransformMatrix(0,1) << ' ' << inverseTransformMatrix(0,2) << ' ' << inverseTransformMatrix(0,3) << ' ' << '\n'
               << inverseTransformMatrix(1,0) << ' ' << inverseTransformMatrix(1,1) << ' ' << inverseTransformMatrix(1,2) << ' ' << inverseTransformMatrix(1,3) << ' ' << '\n'
               << inverseTransformMatrix(2,0) << ' ' << inverseTransformMatrix(2,1) << ' ' << inverseTransformMatrix(2,2) << ' ' << inverseTransformMatrix(2,3) << ' ' << '\n'
               << inverseTransformMatrix(3,0) << ' ' << inverseTransformMatrix(3,1) << ' ' << inverseTransformMatrix(3,2) << ' ' << inverseTransformMatrix(3,3) << ' ' << '\n';
        return textStream.str();
}

string Object::printSumOriginalToTransform() {
    stringstream textStream;
    textStream << fixed << setprecision(10) << sumOriginalToTransform;
    return textStream.str();
}

string Object::printSumOriginalToOriginal() {
    stringstream textStream;
    textStream << fixed << setprecision(10)  << sumOriginalToOriginal;
    return textStream.str();
}

string Object::printTransformFile() {
    string text = "# Transformation matrix\n";
    text = text + printTransformMatrix()+'\n';
    text = text+"# Inverse transformation matrix\n";
    text = text + printInverseTransformMatrix()+'\n';
    text = text+"# Sum absolute translations from original to transformed\n";
    text = text + printSumOriginalToTransform()+'\n'+'\n';
    text = text+"# Sum absolute translations from original to transformed to \"original\"\n";
    text = text + printSumOriginalToOriginal();
    return text;
}

string Object::printFile() {
    //Print all information in the format of a .obj file
    string text = copyrightInfo;
    for (Vertex v : verticies) {
        text = text + v.getString()+'\n';
    }
    if (smoothingGroup)  
      text = text + "s on\n";
    else
        text = text + "s off\n";
    for (Face f : faces) {
        text = text + f.getString()+'\n';
    }
    return text;
    
}
void Object::addVertex(double x, double y, double z) {
    verticies.push_back(Vertex(x,y,z));
}

void Object::addFace(int a, int b, int c, int id) {
    //cout << "vertex size: "<<verticies.size() <<" a: "<<a<<"b: "<<b<<"c: "<<c<<'\n';
    //Face newFace = Face(a-1, b-1, c-1,id);
    int i = faces.size();
    faces.push_back(Face(a-1, b-1, c-1,id));
    verticies[a-1].faces.push_back(i);
    verticies[b-1].faces.push_back(i);
    verticies[c-1].faces.push_back(i);
    //cout << newFace.getPrettyString() << '\n';
    
}

void Object::addFace(int a, int b, int c, int id, int e, int f, int g) {
    //cout << "vertex size: "<<verticies.size() <<" a: "<<a<<"b: "<<b<<"c: "<<c<<'\n';
    //Face newFace = Face(a-1, b-1, c-1,id);
    int i = faces.size();
    faces.push_back(Face(a-1, b-1, c-1,id, e-1, f-1, g-1));
    verticies[a-1].faces.push_back(i);
    verticies[b-1].faces.push_back(i);
    verticies[c-1].faces.push_back(i);
    //cout << newFace.getPrettyString() << '\n';
    
}
/*
void Object::addFace(int a, int normA, int b, int normB, int c, int normC) {
    faces.push_back(Face(a,normA,b,normB,c,normC));
}*/






