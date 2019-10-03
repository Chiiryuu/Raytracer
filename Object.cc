#include "Object.h"
#include <iostream>
#include <fstream>
#include <string>
#include <regex> 
#include <cmath>
#include <sstream>
#include <iomanip>

using namespace std;



Object::Object(string fileName): file(fileName){
    name = file;
    //Remove .obj
    name.erase (name.length()-4,5);
    loadFromFile();
    makeMatrix();
}

Object::Object() : file("None"s){
    name = file;
}


Object::~Object() {
    //delete self;
}

void Object::updateVerts() {
    for (int i=0;i<verticies.size();i++) {
        verticies[i].update(vertMatrix(0,i),vertMatrix(1,i),vertMatrix(2,i));
    }
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
    //Open file stream
    ifstream rawFile(file);
    string line;
    //Check if file was opened
    if (!rawFile)
        cout << "File '"<<file<<"' not found\n";
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
        else if (line.at(0)=='s' && line.at(1)==' ') {
            //Set smoothing group from s line
            if (line.at(line.length()-1) == 'n')
                smoothingGroup = true;
        }
        else if (line.at(0)=='f' && line.at(1)==' ') {
            //Get face verts and face normals from line
            int verts[3];
            int norms[3];
            lineStream>>verts[0];
            lineStream.ignore(2);
            lineStream>>norms[0];
            
            lineStream>>verts[1];
            lineStream.ignore(2);
            lineStream>>norms[1];
            
            lineStream>>verts[2];
            lineStream.ignore(2);
            lineStream>>norms[2];
            
            addFace(verts[0],norms[0],verts[1],norms[1],verts[2],norms[2]);
        }
        else
            //Assume text is copyright info
            copyrightInfo = copyrightInfo+line+'\n';
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

void Object::addFace(int a, int b, int c) {
    faces.push_back(Face(a,b,c));
}

void Object::addFace(int a, int normA, int b, int normB, int c, int normC) {
    faces.push_back(Face(a,normA,b,normB,c,normC));
}








