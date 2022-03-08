#include "Object.h"
#include "Ray.h"
#include "Camera.h"
#include "Sphere.h"
#include "Light.h"
#include "Vertex.h"
#include "Face.h"
#include "Eigen/Dense"
#include <cmath>
#include <utility>
#include <map>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <ctime>

using namespace std;

vector<Object> objects;

vector<Sphere> spheres;

vector<Light> lights;

vector<Camera> cameras;

map<string,int> objectsProcessed;

vector<double> ambient;

const double sphereAlpha = 16.0;

const double minT = 0.00001;

string arrToString(vector<double> vect) {
    string arr = "[ ";
    for (double d : vect)
        arr = arr + to_string(d)+' ';
    return arr+']';
}


double degreeToRad(double degrees) {
    return degrees * M_PI / 180.0;
}

double dotProduct(vector<double> a, vector<double> b) {
    return (a[0]*b[0]) + (a[1]*b[1]) + (a[2]*b[2]);
}

vector<double> pairwiseProduct(vector<double> a, vector<double> b) {
    vector<double> result;
    result.push_back(a[0]*b[0]);
    result.push_back(a[1]*b[1]);
    result.push_back(a[2]*b[2]);
    return result;
}

vector<double> vectorAddition(vector<double> a, vector<double> b) {
    vector<double> result;
    result.push_back(a[0]+b[0]);
    result.push_back(a[1]+b[1]);
    result.push_back(a[2]+b[2]);
    return result;
}

vector<double> vectorSubtraction(vector<double> a, vector<double> b) {
    vector<double> result;
    result.push_back(a[0]-b[0]);
    result.push_back(a[1]-b[1]);
    result.push_back(a[2]-b[2]);
    return result;
}

vector<double> crossProduct(vector<double> a, vector<double> b) {
    vector<double> result;
    result.push_back(   a[1]*b[2]  -   a[2]*b[1]  );
    result.push_back(   a[2]*b[0]  -   a[0]*b[2]  );
    result.push_back(   a[0]*b[1]  -   a[1]*b[0]  );
    return result;
}

double magnitude(vector<double>& vect){
    double mag = 0;
    for (double d : vect)
        mag = mag + (d*d);
    
    return pow(mag,0.5);
}

vector<double> scalarMult(vector<double> vect, double scalar) {
    double mag = magnitude(vect);
    
    vector<double> normVect;
    for (double d : vect)
        normVect.push_back(d * scalar);
    
    return normVect;
}

vector<double> normVect(vector<double> vect) {
    double mag = magnitude(vect);
    
    vector<double> normVect;
    for (double d : vect)
        normVect.push_back(d/mag);
    
    return normVect;
}

vector<double> nearestAxis(vector<double> vect) {
    double min = 99999999;
    int minElement = 0;
    for (int i=0; i<vect.size();i++)
        if (vect[i] < min) {
            min = vect[i];
            minElement=i;
        }
    vector<double> closest(vect);
    closest[minElement]=1;
    closest = normVect(closest);
    
    return closest;
}

Eigen::MatrixXd zRotMatrix(double rads) {
    Eigen::MatrixXd rotMatrix;
    rotMatrix.resize(4,4);
    rotMatrix(0,0)=cos(rads);
    rotMatrix(0,1)=-1*sin(rads);
    rotMatrix(0,2)=0;
    rotMatrix(0,3)=0;
    
    rotMatrix(1,0)=sin(rads);
    rotMatrix(1,1)=cos(rads);
    rotMatrix(1,2)=0;
    rotMatrix(1,3)=0;
    
    rotMatrix(2,0)=0;
    rotMatrix(2,1)=0;
    rotMatrix(2,2)=1;
    rotMatrix(2,3)=0;
    
    rotMatrix(3,0)=0;
    rotMatrix(3,1)=0;
    rotMatrix(3,2)=0;
    rotMatrix(3,3)=1;
    
    return rotMatrix;
}

Eigen::MatrixXd ScaleToMatrix(vector<double> scale) {
    Eigen::MatrixXd scaleMatrix;
    scaleMatrix.resize(4,4);
    
    scaleMatrix(0,0)=scale[0];
    scaleMatrix(0,1)=0;
    scaleMatrix(0,2)=0;
    scaleMatrix(0,3)=0;
    
    scaleMatrix(1,0)=0;
    scaleMatrix(1,1)=scale[1];
    scaleMatrix(1,2)=0;
    scaleMatrix(1,3)=0;
    
    scaleMatrix(2,0)=0;
    scaleMatrix(2,1)=0;
    scaleMatrix(2,2)=scale[2];
    scaleMatrix(2,3)=0;
    
    scaleMatrix(3,0)=0;
    scaleMatrix(3,1)=0;
    scaleMatrix(3,2)=0;
    scaleMatrix(3,3)=1;
    
    return scaleMatrix;
}

Eigen::MatrixXd TranslationToMatrix(vector<double> translate) {
    Eigen::MatrixXd transMatrix;
    transMatrix.resize(4,4);
    
    transMatrix(0,0)=1;
    transMatrix(0,1)=0;
    transMatrix(0,2)=0;
    transMatrix(0,3)=translate[0];
    
    transMatrix(1,0)=0;
    transMatrix(1,1)=1;
    transMatrix(1,2)=0;
    transMatrix(1,3)=translate[1];
    
    transMatrix(2,0)=0;
    transMatrix(2,1)=0;
    transMatrix(2,2)=1;
    transMatrix(2,3)=translate[2];
    
    transMatrix(3,0)=0;
    transMatrix(3,1)=0;
    transMatrix(3,2)=0;
    transMatrix(3,3)=1;
    
    return transMatrix;
}

Eigen::MatrixXd MxMatrix(vector<double> Av, vector<double> Bv, vector<double> Cv, vector<double> Dv) {
    Eigen::MatrixXd transMatrix;
    transMatrix.resize(3,3);
    
    transMatrix(0,0)=Av[0]-Bv[0];std::vector<double> crossProduct(std::vector<double> a, std::vector<double> b);
    transMatrix(0,1)=Av[0]-Cv[0];
    transMatrix(0,2)=Dv[0];
    
    transMatrix(1,0)=Av[1]-Bv[1];
    transMatrix(1,1)=Av[1]-Cv[1];
    transMatrix(1,2)=Dv[1];
    
    transMatrix(2,0)=Av[2]-Bv[2];
    transMatrix(2,1)=Av[2]-Cv[2];
    transMatrix(2,2)=Dv[2];
    
    return transMatrix;
}

Eigen::MatrixXd MyMatrix(vector<double> Av, vector<double> Lv) {
    Eigen::MatrixXd transMatrix;
    transMatrix.resize(3,1);
    
    transMatrix(0,0)=Av[0]-Lv[0];
    
    transMatrix(1,0)=Av[1]-Lv[1];
    
    transMatrix(2,0)=Av[2]-Lv[2];
    
    return transMatrix;
}


Eigen::MatrixXd AxisAngleToMatrix(vector<double> rotAxis, double theta) {
    Eigen::MatrixXd rotMatrix;
    
    if (magnitude(rotAxis) == 0)
        rotAxis[0]=1;
    
    
    //Set up known parts of matrix
    rotMatrix.resize(4,4);
    rotMatrix(3,0)=0;
    rotMatrix(3,1)=0;
    rotMatrix(3,2)=0;
    rotMatrix(3,3)=1;
    rotMatrix(0,3)=0;
    rotMatrix(1,3)=0;
    rotMatrix(2,3)=0;
    
    //Normalize axis
    vector<double> w = normVect(rotAxis);
    
    //Add w to matrix
    rotMatrix(2,0)=w[0];
    rotMatrix(2,1)=w[1];
    rotMatrix(2,2)=w[2];
    
    vector<double> m = nearestAxis(rotAxis);
    
    vector<double> u = crossProduct(w,m);
    rotMatrix(0,0)=u[0];
    rotMatrix(0,1)=u[1];
    rotMatrix(0,2)=u[2];
    
    vector<double> v = crossProduct(w,u);
    rotMatrix(1,0)=v[0];
    rotMatrix(1,1)=v[1];
    rotMatrix(1,2)=v[2];
    
    //cout<<"Theta:\t"<<theta<<'\n';
    //cout<<"Axis:\t"<<arrToString(rotAxis)<<'\n';
    //cout<<"Norm:\t"<<arrToString(normVect(rotAxis))<<'\n';
    
    Eigen::MatrixXd zRot = zRotMatrix(degreeToRad(theta));

    Eigen::MatrixXd invRotMatrix = rotMatrix.inverse();
    
    
    return invRotMatrix * zRot * rotMatrix;
}

void setAmbient(string input) {
    istringstream lineStream(input);
    double temp;
    lineStream>>temp;
    ambient.push_back(temp);
    lineStream>>temp;
    ambient.push_back(temp);
    lineStream>>temp;
    ambient.push_back(temp);
    

}



void makeCamera(string input) {
    istringstream lineStream(input);
    
    vector<double> eye;
    
    double x;
    lineStream>>x;
    double y;
    lineStream>>y;
    double z;
    lineStream>>z;
    
    double lookX;
    lineStream>>lookX;
    double lookY;
    lineStream>>lookY;
    double lookZ;
    lineStream>>lookZ;
    
    double upX;
    lineStream>>upX;
    double upY;
    lineStream>>upY;
    double upZ;
    lineStream>>upZ;
    
    vector<double> up{upX, upY, upZ};
    up = normVect(up);
    
    double d;
    lineStream>>d;
    
    double bounds0;
    lineStream>>bounds0;
    double bounds1;
    lineStream>>bounds1;
    double bounds2;
    lineStream>>bounds2;
    double bounds3;
    lineStream>>bounds3;
    
    double resX;
    lineStream>>resX;
    double resY;
    lineStream>>resY;
    
    vector<double> w{x-lookX, y-lookY, z-lookZ};
    w = normVect(w);
    
    vector<double> u = crossProduct(up, w);
    u = normVect(u);
    
    vector<double> v = crossProduct(w, u);
    
    Camera camera(x, y, z, lookX, lookY, lookZ, up[0], up[1], up[2], d, bounds0, bounds1, bounds2, bounds3, resX, resY, w[0], w[1], w[2], u[0], u[1], u[2], v[0], v[1], v[2]);
    cameras.push_back(camera);
    

}

void makeLight(string input) {
    istringstream lineStream(input);
    
    vector<double> pos;
    
    double temp;
    lineStream>>temp;
    pos.push_back(temp);
    lineStream>>temp;
    pos.push_back(temp);
    lineStream>>temp;
    pos.push_back(temp);
    lineStream>>temp;
    pos.push_back(temp);
    
    double colorR;
    lineStream>>colorR;
    double colorG;
    lineStream>>colorG;
    double colorB;
    lineStream>>colorB;
    
    //Temp = w for the light source at this time
    if (temp==0) {
        pos = scalarMult(normVect(pos), 99999999.0);
    }
    
    Light newLight(pos[0], pos[1], pos[2], pos[3], colorR, colorG, colorB);
    lights.push_back(newLight);
    

}

void makeSphere(string input) {
    istringstream lineStream(input);
    
    double x;
    lineStream>>x;
    double y;
    lineStream>>y;
    double z;
    lineStream>>z;
    
    double r;
    lineStream>>r;
    
    double ambientR;
    lineStream>>ambientR;
    double ambientG;
    lineStream>>ambientG;
    double ambientB;
    lineStream>>ambientB;
    
    double diffuseR;
    lineStream>>diffuseR;
    double diffuseG;
    lineStream>>diffuseG;
    double diffuseB;
    lineStream>>diffuseB;
    
    double specularR;
    lineStream>>specularR;
    double specularG;
    lineStream>>specularG;
    double specularB;
    lineStream>>specularB;
    
    double attenuationR;
    lineStream>>attenuationR;
    double attenuationG;
    lineStream>>attenuationG;
    double attenuationB;
    lineStream>>attenuationB;
    
    double Ni=0.0;
    lineStream>>Ni;
    
    string textureFile="";
    lineStream>>textureFile;
    
    int tiling=0;
    lineStream>>tiling;
    
    double phi = 0.0;
    lineStream>>phi;
    phi = degreeToRad(phi);
    
    double theta = 0.0;
    lineStream>>theta;
    theta = degreeToRad(theta);
    
    //string bumpFile="";
    //lineStream>>bumpFile;
    
    
    
    //cout<<a<<'\n';
    
    //Sphere newSphere(x, y, z, r, ambientR, ambientG, ambientB, diffuseR, diffuseG, diffuseB, specularR, specularG, specularB, attenuationR, attenuationG, attenuationB, a);
    Sphere newSphere(x, y, z, r, ambientR, ambientG, ambientB, diffuseR, diffuseG, diffuseB, specularR, specularG, specularB, attenuationR, attenuationG, attenuationB, Ni, textureFile, phi, theta, tiling);
    spheres.push_back(newSphere);
    

}


void transform(string input) {
    
    //cout<<input<<'\n';
    
    //cout<<input<<'\n';
    istringstream lineStream(input);
    
    //Temp double
    double tempD;
    
    //Other useful declarations
    vector<double> rotAxis;
    double theta;
    double scale;
    vector<double> trans;
    vector<double> scaleVector;
    double smoothing;
    string objectName;
    
    //Build rotation axis vector
    lineStream>>tempD;
    rotAxis.push_back(tempD);
    lineStream>>tempD;
    rotAxis.push_back(tempD);
    lineStream>>tempD;
    rotAxis.push_back(tempD);
    
    //Add theta
    lineStream>>theta;
    
    //Build scale vector (uniform now, can easily be expanded)
    lineStream>>scale;
    scaleVector.push_back(scale);
    scaleVector.push_back(scale);
    scaleVector.push_back(scale);
    
    //Build translation vector
    lineStream>>tempD;
    trans.push_back(tempD);
    lineStream>>tempD;
    trans.push_back(tempD);
    lineStream>>tempD;
    trans.push_back(tempD);
    
    //Set smoothing cutoff
    lineStream>>smoothing;
    smoothing = degreeToRad(smoothing);
    
    //Set object name
    lineStream>>objectName;
    
    int tiling = 0;
    lineStream>>tiling;
    
    //cout << objectName<<'\n';
    ifstream checkFileExists("models/"+objectName);
    if (!checkFileExists.good()){
        cout<<"Unable to find model file '"<<objectName<<"'\n";
        return;
    }
        
    
    //Make object
    Object object(objectName, tiling, smoothing);
    
    //Update name following naming conventions in assignment, with help of a map
    string newName = object.getName();
    if(objectsProcessed.insert(make_pair(newName , 0)).second == false)
    {
        objectsProcessed[newName] =   objectsProcessed[newName]+1;
    }
    int numObject = objectsProcessed[newName];
     newName = newName + "_mw";
    if (numObject<10)
        newName = newName+"0";
    newName = newName+to_string(numObject);
    object.setName(newName);
    
    
    cout<<"Transforming '"<<newName<<"'\n";
    
    
    //the order of transformation is rotate, scale, and then translate
        //SO REVERSE ORDER IN MATRIX COMPUTATION
    Eigen::MatrixXd R;
    if (theta == 0)
        R = Eigen::Matrix<double, 4, 4>::Identity();
    else
        R = AxisAngleToMatrix(rotAxis,theta);
    
    Eigen::MatrixXd T = TranslationToMatrix(trans);
    
    Eigen::MatrixXd S = ScaleToMatrix(scaleVector);
    
    object.transformMatrix = T * S * R;
    
    //cout<<"Rounded Transform\n"<<object.printTransformMatrix()<<'\n';
    
    //cout<<"Transform Matrix:\n"<<object.transformMatrix<<'\n';
    
    object.inverseTransformMatrix = object.transformMatrix.inverse();
    
    //cout<<"Rounded Inverse Transform\n"<<object.printInverseTransformMatrix()<<'\n';
    
    //cout<<"Inverse Transform Matrix:\n"<<object.inverseTransformMatrix<<'\n';
    
    
    
    object.vertMatrix = object.transformMatrix * object.vertMatrix;
    
    //Calculate absolute transform, original to transform
    object.sumOriginalToTransform=0;
    for (int i=0; i<object.verticies.size(); i++) {
        object.sumOriginalToTransform += ( abs(object.verticies[i].point[0] - object.vertMatrix(0,i)) + abs(object.verticies[i].point[1]  - object.vertMatrix(1,i)) + abs(object.verticies[i].point[2]  - object.vertMatrix(2,i)) );
    }
    
    //cout<<"Sum of Original To Transform:\n"<<object.printSumOriginalToTransform();
    
    //Calculate the same, but then back again. (Not using the object's actual matrix, as we'd need to reset it later)
    Eigen::MatrixXd theoreticalTransform = object.inverseTransformMatrix * object.vertMatrix;
    object.sumOriginalToOriginal=0;
    for (int i=0; i<object.verticies.size(); i++) {
        object.sumOriginalToOriginal += ( abs(object.verticies[i].point[0] - theoreticalTransform(0,i)) + abs(object.verticies[i].point[1]  - theoreticalTransform(1,i)) + abs(object.verticies[i].point[2]  - theoreticalTransform(2,i)) );
    }
    //Round down small values to 0
    if (object.sumOriginalToOriginal <= 0.00000000001)
        object.sumOriginalToOriginal=0.0;
    //cout<<"Sum of Original To Original:\n"<<object.printSumOriginalToOriginal();
    
    
    //This updates the vertex vector in object from its corresponding vertex matrix. Call this after we do calculations.
    object.updateVerts();
    
    objects.push_back(object);
    
    
    //Create Files
    
    
    
}

void printPercent(double percent) {
    cout << fixed << setprecision(2) <<percent  << '%' << '\r' << flush << setprecision(9) ;
}

bool pathBlocked(Ray &ray, double maxDistance=999999999999999999.0) {
    for (int i=0; i<spheres.size();i++) {
        vector<double> rayToSphere = vectorSubtraction(spheres[i].point, ray.origin);
        
        double v = dotProduct(rayToSphere, ray.direction);
        
        double cc = dotProduct(rayToSphere, rayToSphere);
        
        double disc = (pow(spheres[i].radius,2)) - (cc - (pow(v,2)));
        if (disc > 0 ) {
            double d = pow(disc,0.5);
            double distance = v - d;
            if (distance > minT  && distance < maxDistance) {
                return true;
            }
            else {
                distance = v + d; 
                if (distance > minT && distance < maxDistance) {
                    return true;
                }
            }
        }
    }
    
    if (objects.size() == 0)
        return false;
    
    
    vector<double> Lv(ray.origin);
            
    vector<double> Dv(ray.direction);
    
    Eigen::MatrixXd Mx;
    Mx.resize(3,3);
    Mx(0,2)=Dv[0];
    
    Mx(1,2)=Dv[1];
    
    Mx(2,2)=Dv[2];
    
    Eigen::MatrixXd My;
    My.resize(3,1);
    
    
    for (int i=0; i<objects.size();i++) {
        Object& object = objects[i];
        for (int j=0; j<object.faces.size();j++) {
            Face& face = object.faces[j];
            vector<double> Av = object.verticies[face.verticies[0]].point;
            vector<double> Bv = object.verticies[face.verticies[1]].point;
            vector<double> Cv = object.verticies[face.verticies[2]].point;
            

            Mx(0,0)=Av[0]-Bv[0];
            Mx(0,1)=Av[0]-Cv[0];
    
            Mx(1,0)=Av[1]-Bv[1];
            Mx(1,1)=Av[1]-Cv[1];
    
            Mx(2,0)=Av[2]-Bv[2];
            Mx(2,1)=Av[2]-Cv[2];
            
            double detMx = Mx.determinant();
            if (abs(detMx) > 0) {
                My(0,0)=Av[0]-Lv[0];  
                My(1,0)=Av[1]-Lv[1];
                My(2,0)=Av[2]-Lv[2];
                Eigen::MatrixXd x = Mx.fullPivLu().solve(My);
                double B = x(0,0);
                double Y = x(1,0);
                double t = x(2,0);
                if (t>minT && B >= 0 && Y >= 0 && (B+Y) <= 1  && t < maxDistance) {
                    return true;
                }
            }
        }
    }
    
    
    return false;
    
}


int checkRaySphere(Ray &ray, double &minDistance, double maxDistance=999999999999999999.0) {
    if (spheres.size() == 0)
        return -1;
    int index = -1;
    for (int i=0; i<spheres.size();i++) {
        vector<double> rayToSphere = vectorSubtraction(spheres[i].point, ray.origin);
        double v = dotProduct(rayToSphere, ray.direction);
        double cc = dotProduct(rayToSphere, rayToSphere);
        double disc = (pow(spheres[i].radius,2)) - (cc - (pow(v,2)));
        //cout<<disc<<'\n';
        if (disc > 0 ) {
            double d = pow(disc,0.5);
            double distance = v - d;
            if (distance > minT && distance < minDistance && distance < maxDistance) {
                minDistance = distance;
                index = i;
            }
            else {
                distance = v + d; 
                if (distance > minT && distance < minDistance && distance < maxDistance) {
                    minDistance = distance;
                    index = i;
                }
            }
        }
    }
    return index;
    
}

vector<double> smoothedN(Face& face, double A, double B, double Y) {
    vector<double> N{0.0,0.0,0.0};
    N[0] = (A*face.Na[0]) + (B*face.Nb[0]) + (Y*face.Nc[0]);
    N[1] = (A*face.Na[1]) + (B*face.Nb[1]) + (Y*face.Nc[1]);
    N[2] = (A*face.Na[2]) + (B*face.Nb[2]) + (Y*face.Nc[2]);
    return N;
}

void getRefractPoint(Ray &ray, double& minDistance, bool &hitFace, vector<double> &N, int& faceIndex, int i) {
    if (objects.size() == 0)
        return;
    int hitIndex = -1;
    vector<double> Lv(ray.origin); 
    vector<double> Dv(ray.direction);
    
    Eigen::MatrixXd Mx;
    Mx.resize(3,3);
    Mx(0,2)=Dv[0];
    
    Mx(1,2)=Dv[1];
    
    Mx(2,2)=Dv[2];
    
    Eigen::MatrixXd My;
    My.resize(3,1);
    
    Object& object = objects[i];
        for (int j=0; j<object.faces.size();j++) {
            if (j != faceIndex) {
            //cout<<"checking face "<<j<<'\n';
            Face& face = object.faces[j];
            
            vector<double> Av = object.verticies[face.verticies[0]].point;
            vector<double> Bv = object.verticies[face.verticies[1]].point;
            vector<double> Cv = object.verticies[face.verticies[2]].point;
            

            Mx(0,0)=Av[0]-Bv[0];
            Mx(0,1)=Av[0]-Cv[0];
    
            Mx(1,0)=Av[1]-Bv[1];
            Mx(1,1)=Av[1]-Cv[1];
    
            Mx(2,0)=Av[2]-Bv[2];
            Mx(2,1)=Av[2]-Cv[2];
            
            double detMx = Mx.determinant();
            if (abs(detMx) > 0) {
                My(0,0)=Av[0]-Lv[0];  
                My(1,0)=Av[1]-Lv[1];
                My(2,0)=Av[2]-Lv[2];
                Eigen::MatrixXd x = Mx.fullPivLu().solve(My);
                double B = x(0,0);
                double Y = x(1,0);
                double t = x(2,0);
                if (t>minT && B >= 0 && Y >= 0 && (B+Y) <= 1 && t<minDistance) {
                    faceIndex = j;
                    hitFace = true;
                    minDistance = t;
                    hitIndex = i;
                    
                    //cout<<"A: "<<(1-(B+Y))<<"B: "<<B<<"Y: "<<Y<<'\n';
                    
                    if (object.smoothingGroup) {
                        N = smoothedN(face, (1-(B+Y)), B, Y);
                        if (dotProduct(N, Dv) > 0)
                            N = face.NInv;
                    }
                    else {
                        N = face.N;
                        if (dotProduct(N, Dv) > 0)
                            N = face.NInv;
                    }
                    
                    
                    //cout<<N[0]<<' '<<N[1]<<' '<<N[2]<<'\n';
                    
                }
           
            }
            }
        
        }
    
}


int checkRayFace(Ray &ray, double &minDistance, bool &hitFace, vector<double> &N, int &faceIndex, double &finalB, double &finalY, double maxDistance=999999999999999999.0) {
    if (objects.size() == 0)
        return -1;
    int hitIndex = -1;
    vector<double> Lv(ray.origin); 
    vector<double> Dv(ray.direction);
    
    Eigen::MatrixXd Mx;
    Mx.resize(3,3);
    Mx(0,2)=Dv[0];
    
    Mx(1,2)=Dv[1];
    
    Mx(2,2)=Dv[2];
    
    Eigen::MatrixXd My;
    My.resize(3,1);
    
    
    for (int i=0; i<objects.size();i++) {
        Object& object = objects[i];
        for (int j=0; j<object.faces.size();j++) {
            //cout<<"checking face "<<j<<'\n';
            Face& face = object.faces[j];
            
            vector<double> Av = object.verticies[face.verticies[0]].point;
            vector<double> Bv = object.verticies[face.verticies[1]].point;
            vector<double> Cv = object.verticies[face.verticies[2]].point;
            

            Mx(0,0)=Av[0]-Bv[0];
            Mx(0,1)=Av[0]-Cv[0];
    
            Mx(1,0)=Av[1]-Bv[1];
            Mx(1,1)=Av[1]-Cv[1];
    
            Mx(2,0)=Av[2]-Bv[2];
            Mx(2,1)=Av[2]-Cv[2];
            
            double detMx = Mx.determinant();
            if (abs(detMx) > 0) {
                My(0,0)=Av[0]-Lv[0];  
                My(1,0)=Av[1]-Lv[1];
                My(2,0)=Av[2]-Lv[2];
                Eigen::MatrixXd x = Mx.fullPivLu().solve(My);
                double B = x(0,0);
                double Y = x(1,0);
                double t = x(2,0);
                if (t>minT && B >= 0 && Y >= 0 && (B+Y) <= 1 && t<minDistance  && t < maxDistance) {
                    finalB = B;
                    finalY = Y;
                    faceIndex = j;
                    hitFace = true;
                    minDistance = t;
                    hitIndex = i;
                    
                    //cout<<"A: "<<(1-(B+Y))<<"B: "<<B<<"Y: "<<Y<<'\n';
                    
                    if (object.smoothingGroup) {
                        N = smoothedN(face, (1-(B+Y)), B, Y);
                        if (dotProduct(N, Dv) > 0)
                            N = face.NInv;
                    }
                    else {
                        N = face.N;
                        if (dotProduct(N, Dv) > 0)
                            N = face.NInv;
                    }
                    
                    
                    //cout<<N[0]<<' '<<N[1]<<' '<<N[2]<<'\n';
                    
                }
           
            }
        
        }
    }
    
    return hitIndex;
    
}

/*
def refract_tray(self, W, pt, N, eta1, eta2) :
        etar  = eta1 / eta2
        a     = - etar
        wn    = np.dot(W,N)
        radsq = etar**2 * (wn**2 - 1) + 1
        if (radsq < 0.0) :
           T = np.array([0.0,0.0,0.0])
        else :
           b = (etar * wn) - sqrt(radsq)
           T = a * W + b * N 
        return(T)
        
    def refract_exit(self, W, pt, eta_in, eta_out) :
        T1 = self.refract_tray(W, pt, make_unit(pt - self.C), eta_out, eta_in) 
        if (sum(T1) == 0.0) :
            return None
        else :
            exit = pt + 2 * np.dot((self.C - pt),T1) * T1
            Nin = make_unit(self.C - exit)
            T2 = self.refract_tray(-T1, exit, Nin, eta_in, eta_out)
            refR = Ray(exit, T2)
            return refR
            */

vector<double> refractRay(vector<double> W, vector<double> N, double etaIn, double etaOut) {
    
    
    
    double etaRatio = etaIn / etaOut;
    double a = 0 - etaRatio;
    double WN = dotProduct(W,N);
    double radSq = (pow(etaRatio, 2) * (pow(WN,2) -1 )) + 1;
    if (radSq < 0.0) {
        return vector<double>{0.0, 0.0, 0.0};
    }
    double b = (etaRatio * WN) - sqrt(radSq);
    vector<double> T = vectorAddition(scalarMult(W,a),scalarMult(N,b));
    return T;
    
}

Ray refractExit(Sphere& sphere, vector<double> W, vector<double> pt, vector<double> N, double etaIn, double etaOut) {
   //cout<<"W: "<<W[0]<<' '<<W[1]<<' '<<W[2]<<" pt: "<<pt[0]<<' '<<pt[1]<<' '<<pt[2]<<'\n';
   vector<double> T1 = refractRay(W, N, etaOut, etaIn);
   if (magnitude(T1) < 0.0001)
       return Ray();
    //cout<<"T1: "<<T1[0]<<' '<<T1[1]<<' '<<T1[2]<<'\n';
    vector<double> exit = vectorAddition(pt, scalarMult(T1, 2.0 * dotProduct(vectorSubtraction(sphere.point, pt),T1)));
    vector<double> Nin = normVect(vectorSubtraction(sphere.point, exit));
    vector<double> T2 = refractRay(scalarMult(T1,-1), Nin, etaIn, etaOut);
    //cout<<"T2: "<<T2[0]<<' '<<T2[1]<<' '<<T2[2]<<'\n';
    return Ray(exit[0],exit[1],exit[2],T2[0],T2[1],T2[2]);
    
}



//Check for sphere collision, calculate color
vector<double> checkCollision(Ray& ray, vector<double> totalColor, vector<double> totalAttenuation, int recursions=0, double etaOut=1.0) {
    double minDistance = 9999999;
    int hitIndex = -1;
    int faceIndex = -1;
    bool hitObject = false;
    bool reflect = true;
    bool refract = false;
    vector<double> hitPoint;
    vector<double> sphereCenter;
    vector<double> N;
    double etaIn;
    double B;
    double Y;
    //vector<double> face;
    
    vector<double>& Lv(ray.origin);
            
    vector<double>& Dv(ray.direction);
    
    //cout<<"[ "<<Dv[0]<<' '<<Dv[1]<<' '<<Dv[2]<<" ]\n";
    //cout<<"[ "<<totalColor[0]<<' '<<totalColor[1]<<' '<<totalColor[2]<<" ]\n";
    //cout<<"[ "<<totalAttenuation[0]<<' '<<totalAttenuation[1]<<' '<<totalAttenuation[2]<<" ]\n";
    
    int sphereIndex = checkRaySphere(ray, minDistance);
    
    int objectIndex = checkRayFace(ray, minDistance, hitObject, N, faceIndex, B, Y);
    
    if (hitObject)
        hitIndex = objectIndex;
    else
        hitIndex = sphereIndex;

    
    
    /*std::vector<std::string> materialNames;
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
    
    
    
    if (hitIndex > -1) {
        //cout<<hitIndex<<'\n'<<'\n';
        vector<double> color;
        //vector<double> ambient;
        vector<double> specular;
        vector<double> attenuation;
        vector<double> diffuse;
        double alpha=0;
        double etaIn = 0;
        hitPoint = vectorAddition(Lv, scalarMult(Dv, minDistance));
        ray.setHit();
        ray.setHitPosition(hitPoint[0],hitPoint[1],hitPoint[2]);
        
        //cout << "min: " << minDistance << '\n';
        if (!hitObject) {
            Sphere& hitSphere = spheres[hitIndex];
            color = pairwiseProduct(ambient,hitSphere.ambient);
            diffuse = hitSphere.diffuse;
            specular = hitSphere.specular;
            attenuation = hitSphere.attenuation;
            
            if (hitSphere.textured) {
                vector<double> faceTexture = hitSphere.pointToTexture(hitPoint[0], hitPoint[1], hitPoint[2]);
            
                color = pairwiseProduct(color, faceTexture);
                
                diffuse = pairwiseProduct(faceTexture ,diffuse);
                
                specular = pairwiseProduct(faceTexture ,diffuse);
                
            }
            
            
            if (magnitude(attenuation) < minT)
                reflect = false;
            sphereCenter = hitSphere.point;
            N = vectorSubtraction(hitPoint, hitSphere.point);
            if (dotProduct(N, Dv) > 0)
                N = scalarMult(N,-1);
            N = normVect(N);
            alpha=sphereAlpha;
            etaIn = hitSphere.Ni;
            if (etaIn > 0.000001)
                refract = true;
        }
        else {
            //color = pairwiseProduct(ambient,vector<double>{0.2, 0.2, 0.2});
            Object& object = objects[hitIndex];
            Face& face = object.faces[faceIndex];
            int matIndex = face.matIndex;
            
            vector<double> faceAmbient;
            faceAmbient.push_back(object.aR[matIndex]);
            faceAmbient.push_back(object.aG[matIndex]);
            faceAmbient.push_back(object.aB[matIndex]);
            
            color = pairwiseProduct(ambient,faceAmbient);
            
            diffuse.push_back(object.dR[matIndex]);
            diffuse.push_back(object.dG[matIndex]);
            diffuse.push_back(object.dB[matIndex]);
            
            specular.push_back(object.sR[matIndex]);
            specular.push_back(object.sG[matIndex]);
            specular.push_back(object.sB[matIndex]);
            
            
            
            if (face.textured) {
                
                
                //cout<<"Matindex: "<<matIndex<<'\n';
                
                //cout<<"Object textures: "<<object.textureFiles.size()<<'\n';
                
                double A = 1-(B+Y);
                //B, Y
                Texture& texture = object.textureFiles[object.textureIndex[matIndex]];
                
                //cout<<"found texture\n"<<flush;
                
                double Au = object.u[face.texture[0]];
                
                //cout<<"found Au\n"<<flush;
                
                double Av = object.v[face.texture[0]];
                
                //cout<<"found Av\n"<<flush;
                
                double Bu = object.u[face.texture[1]];
                double Bv = object.v[face.texture[1]];
                double Yu = object.u[face.texture[2]];
                double Yv = object.v[face.texture[2]];
                
                double u = (Au*A) + (Bu * B) + (Yu * Y);
                double v = (Av*A) + (Bv * B) + (Yv * Y);
                
                //cout<<"Au: "<<Au<<", Av: "<<Av<<"\n"<<flush;
                
                vector<double> faceTexture = texture.get(u, v);
            
                color = pairwiseProduct(color, faceTexture);
                
                diffuse = pairwiseProduct(faceTexture ,diffuse);
                
                specular = pairwiseProduct(faceTexture ,diffuse);
            }
            
            //cout<<"B: "<<B<<", Y: "<<Y<<'\n';
            
            
            if (object.illuminationModels[matIndex] == 2) {
                attenuation  = vector<double>{0,0,0};
                reflect = false;
            }
            else if (object.illuminationModels[matIndex] == 6) {
                attenuation  = specular;
                refract = true;
                etaIn = object.Ni[matIndex];
                //etaIn = 1.5;
                //cout<<"mat index: "<<matIndex<<'\n';
                //cout<<"read eta: "<<object.Ni[matIndex]<<'\n';
            }
            else {
                //attenuation  = vector<double>{1, 1, 1};
                attenuation  = specular;
            }
            alpha = object.specularExponents[matIndex];
        }
        
       for (Light l:lights) {
           vector<double> toL = vectorSubtraction(l.point, hitPoint);
           double distanceToL= magnitude(toL);
           toL = normVect(toL);
           double dotProd = dotProduct(N, toL);  
           
           Ray testRay(hitPoint[0],hitPoint[1],hitPoint[2], toL[0], toL[1], toL[2]);
           
           
           if (dotProd > 0.0 && !pathBlocked(testRay, distanceToL)) {
               //Diffuse
               color = vectorAddition(color, scalarMult(pairwiseProduct(diffuse, l.color), dotProd)); 
               
               //Specular
               vector<double> toC = normVect(vectorSubtraction(Lv,hitPoint));
               vector<double> spR = normVect(vectorSubtraction(scalarMult(N,dotProd*2),toL));
               double CdR = dotProduct(toC, spR);
               if (CdR > 0.0)
                   color = vectorAddition(color, scalarMult(pairwiseProduct(specular, l.color), pow(CdR,alpha))); 
           }
           
           
       }
        
        totalColor = vectorAddition(totalColor, pairwiseProduct(totalAttenuation, color));
        
        if (recursions>0 && (reflect || refract) ) {
           vector<double> Uinv= scalarMult(Dv,-1.0);
            if (reflect) {
                vector<double> R = normVect(vectorSubtraction((scalarMult(N,2*dotProduct(N,Uinv))),Uinv));
                if (dotProduct(N,R) > 0.001) {
                    Ray newRay(hitPoint[0],hitPoint[1],hitPoint[2],R[0],R[1],R[2]);
                    totalColor = checkCollision(newRay, totalColor, pairwiseProduct(totalAttenuation, attenuation), recursions-1);
                }
            }
            
            if (refract) {
            //fraR = ray.best_sph.refract_exit(-1 * ray.D, ray.best_pt, mat.eta, eta_outside)
                if (!hitObject) {
                //cout<<"hit: "<<hitPoint[0]<<' '<<hitPoint[1]<<' '<<hitPoint[2]<<'\n';
                Ray refRay = refractExit(spheres[hitIndex], Uinv, hitPoint, N, etaIn, etaOut);
                if (!refRay.invalid) {
                    //cout<<"exit: "<<refRay.origin[0]<<' '<<refRay.origin[1]<<' '<<refRay.origin[2]<<" dv: "<<refRay.direction[0]<<' '<<refRay.direction[1]<<' '<<refRay.direction[2]<<'\n'<<'\n';
                    vector<double> refractVal = vectorSubtraction(vector<double>{1.0,1.0,1.0}, attenuation);
                    totalColor = checkCollision(refRay, totalColor, pairwiseProduct(totalAttenuation, refractVal), recursions-1); 
                    }
                }
                else {
                    //cout<<"Refracting on object...\n";
                    
                    //cout<<"etaIn: "<<etaIn<<", etaOut: "<<etaOut<<'\n';
                    
                    vector<double> T1 = refractRay(Uinv, N, etaOut, etaIn);
                    if (magnitude(T1) < 0.0001)
                        return totalColor;
                    
                    double minDistance = 999999999999999999.999999;
                    bool hitFace = false;
                    Ray insideRay(hitPoint[0], hitPoint[1], hitPoint[2], T1[0], T1[1], T1[2]);
                    getRefractPoint(insideRay, minDistance, hitFace, N, faceIndex, hitIndex);
                    if (!hitFace)
                        return totalColor;
                    hitPoint = vectorAddition(hitPoint, scalarMult(T1, minDistance));
                    vector<double> Tinv = scalarMult(T1,-1.0);
                    vector<double> T2 = refractRay(Tinv, N, etaIn, etaOut);
                    
                    Ray refRay(hitPoint[0], hitPoint[1], hitPoint[2], T2[0], T2[1], T2[2]);
                    
                    vector<double> refractVal = vectorSubtraction(vector<double>{1.0,1.0,1.0}, attenuation);
                    
                    totalColor = checkCollision(refRay, totalColor, pairwiseProduct(totalAttenuation, refractVal), recursions-1); 
                    
                }
            }
        }
    }
            
    return totalColor;
}



void calculateRayPixelColor(int* colorIn, Ray& ray, int recursionLevel=0) {
    
    vector<double> color{0, 0, 0};
    
    color = vectorAddition(color, checkCollision(ray,color, vector<double>{1.0,1.0,1.0}, recursionLevel));
    //color = vectorAddition(color, checkObjectCollision(ray));
    
    /*
    for (int i=0; i<color.size();i++) {
        /*
            if (color[i] < 0)
                color[i]=0;
            else if (color[i]>1)
                color[i]=1;
        color[i] = (int)
        }*/
    
    
    
    colorIn[0] = (int)(max(0,min(255,(int)ceil(255.0*color[0]))));
    colorIn[1] = (int)(max(0,min(255,(int)ceil(255.0*color[1]))));
    colorIn[2] = (int)(max(0,min(255,(int)ceil(255.0*color[2]))));
    
    
}

void makeProgressBar() {
    cout<<"[--------------------------------------------------]"<<'\r'<<'['<<flush;
}
void finishProgressBar() {
    cout<<'\r'<<"[==================================================]                     "<<'\n';
}


void rayTrace(string fileName, int recursionLevel=0, bool multithread=false) {
    Camera camera = cameras[0];
    
    double left = camera.bounds[0];
    double right = camera.bounds[1];
    double bottom = camera.bounds[2];
    double top = camera.bounds[3];
    int width = camera.res[0];
    int height = camera.res[1];
    int total = width*height;
    double near = camera.d;
    
    if (multithread && total < 50) {
        cout<<"Multithread disabled; resolution too low to benefit."<<'\n';
        multithread=false;
    }
    
    if (multithread && total > 697650) {
        cout<<"Multithread disabled; memory overhead too high. Continue or lower resolution."<<'\n';
        multithread=false;
    }
    
    cout<<"Simulating rays...\n";
    
    string imageString = "P3\n"s +  to_string(width)+" "s+to_string(height)+"\n255\n"s;
    
    //SAVE MEMORY BY ONLY SPLITTING ROWS UP FOR PARALLEL
    
    if (multithread) {
        makeProgressBar();
        int divNum = total / 50;
        int numDone = 0;
        
        int red[total]; 
        int blue[total]; 
        int green[total]; 
        
        #pragma omp parallel for
        for (int t=0;t<total;++t) {
            double orgX;
            double orgY;
            double orgZ;
            double dirX;
            double dirY;
            double dirZ;
            
            int i = t / width;
            int j = t % width;
            
            double px = (j*1.0)/(width-1)*(right-left)+left;
            double py = (i*1.0)/(height-1)*(bottom-top)+top;
            
            vector<double> origin = vectorAddition( vectorAddition(camera.eye, scalarMult(camera.w, near)), vectorAddition(scalarMult(camera.u, px), scalarMult(camera.v, py))   );
            
            
            vector<double> dir = vectorSubtraction(origin, camera.eye);
            dir = normVect(dir);
            
            Ray ray(origin[0],origin[1],origin[2],dir[0],dir[1],dir[2]);
            
            int color[3];
            calculateRayPixelColor(color, ray, recursionLevel);
            red[t] = (color[0]);
            blue[t] = (color[1]);
            green[t] = (color[2]);
            
            
            ++numDone;
            if (numDone%divNum == 0)
                cout<<'='<<flush;
        }
        finishProgressBar();
        
        cout<<"Building image file '"+fileName+"'...\n";
        makeProgressBar();
        numDone = 0;
        for (int i=0;i<total;i++) {
        //cout << rays[i].getString() << '\t';
        //double hit = checkSphereCollision(rays[i]);
        imageString = imageString + to_string(red[i]) + " "s + to_string(blue[i]) + " "s + to_string(green[i])+ " "s ;
            
        if ((i+1)%(width) == 0)
            imageString += '\n';
        ++numDone;
        if (numDone%divNum == 0)
                cout<<'='<<flush;
        }
        
        finishProgressBar();
    }
    
    else {
        printPercent(0);
        for (int t=0;t<total;++t) {
            double orgX;
            double orgY;
            double orgZ;
            double dirX;
            double dirY;
            double dirZ;
            
            int i = t / width;
            int j = t % width;
            
            double px = (j*1.0)/(width-1)*(right-left)+left;
            double py = (i*1.0)/(height-1)*(bottom-top)+top;
            
            vector<double> origin = vectorAddition( vectorAddition(camera.eye, scalarMult(camera.w, near)), vectorAddition(scalarMult(camera.u, px), scalarMult(camera.v, py))   );
            
            
            vector<double> dir = vectorSubtraction(origin, camera.eye);
            dir = normVect(dir);
            
            Ray ray(origin[0],origin[1],origin[2],dir[0],dir[1],dir[2]);
            
            int color[3];
            calculateRayPixelColor(color, ray, recursionLevel);
            
            imageString = imageString + to_string((color[0])) + " "s + to_string((color[1])) + " "s + to_string((color[2]))+ " "s ;
            
            printPercent(((t*1.0)/total)*100);
            
            }
    }
    
    
    
    
    
    
    cout<<"Image file '"+fileName+"' successfully generated.\n";
    
    /*
    def pixelRay(i,j):
    px = i/(width-1)*(right-left)+left;
    #py = j/(height-1)*(top-bottom)+bottom;
    py = j/(height-1)*(bottom-top)+top;
    pixpt = EV + (near * WV) + (px * UV) + (py * VV);
    shoot = pixpt - EV; shoot = shoot / shoot.norm();
    raypt = pixpt + shoot * abs(far-near);
    return arrow3d(pixpt, raypt, width=16, color=pixcolor(i,j));
    */
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    //cout<<imageString<<'\n';
    
    ofstream writeObject;
    writeObject.open(fileName, ofstream::out | ofstream::trunc);
    writeObject << imageString;
    writeObject.close();
    
    if(!writeObject)
      std::cout<<"error: failed to write\n";
    
    
}




/*

Main method for assignment 1; may be useful

int main(int argc, char *argv[]) {
    
    if (argc<2) {
        cerr<<"No driver file specified. Usage: "s+argv[0]+" driver_file\n";
        return -1;
    }
    
    string driverPath = argv[1];
    
    ifstream driverFile(driverPath);
    if (!driverFile) {
        cerr<<"Cannot access driver file '"s+argv[1]+"'\n";
        return -1;
    }
    
    cout<<"Opening "<<driverPath<<'\n';
    
    string line;
    string word;
    string driverName = driverPath.substr(0,driverPath.length()-4);
    
    while (getline(driverFile,line)) {
        //Open string stream
        istringstream lineStream(line);
        //Check if line is a transform
        if (lineStream>>word) {
            //Ignore comment lines
            if (word=="#");
            //Transform models with model flag
            if (word=="model") {
                transform(line.substr(6));
            }
        }
    }
    
    mkdir(driverName.c_str(),0777);
    
    for (Object o:objects) {
        cout<<"Writing '"<<o.name<<"'  to file '"<<driverName<<'/'<<o.name<<".obj'\n";
        //cout<<o.name<<'\n';
        //cout<<o.printFile()<<'\n';
        //cout<<o.printTransformFile()<<'\n';
        
        ofstream writeObject;
        writeObject.open(driverName+"/"+o.name+".obj");
        writeObject << o.printFile();
        writeObject.close();
        string transformName = o.name.substr(0,o.name.length()-4);
        transformName+="transform_"+o.name.substr(o.name.length()-4);
        
        cout<<"Writing '"<<o.name<<"' transformation data to file '"<<driverName<<'/'<<transformName<<".txt'\n";
        
        writeObject.open(driverName+"/"+transformName+".txt");
        writeObject << o.printTransformFile();
        writeObject.close();
        
        
    }
    
}
*/


int main(int argc, char *argv[]) {
    
    if (argc<2) {
        cerr<<"No driver file specified. Usage: "s+argv[0]+" driver_file output_file\n";
        return -1;
    }
    
    if (argc<3) {
        cerr<<"No output file specified. Usage: "s+argv[0]+" driver_file output_file\n";
        return -1;
    }
    
    bool multithread = true;
    
    if (argc==4) {
        multithread = false;
        cout<<"Program executing as a single thread. To use multithreading, run without an additional flag.\n";
    }
    
    
    
    string driverPath = argv[1];
    
    ifstream driverFile(driverPath);
    if (!driverFile) {
        cerr<<"Cannot access driver file '"s+argv[1]+"'\n";
        return -1;
    }
    
    cout<<"Opening "<<driverPath<<'\n';
    
    string line;
    string word;
    string driverName = driverPath.substr(0,driverPath.length()-4);
    
    string eye;
    string look;
    string up;
    string d;
    string bounds;
    string res;
    
    int recursionLevel = 0;
    
    clock_t initial = clock();
    
    
    while (getline(driverFile,line)) {
        //Open string stream
        istringstream lineStream(line);
        //Check if line is a transform
        if (lineStream>>word) {
            //Ignore comment lines
            if (word=="#");
            //Transform models with model flag
            else if (word=="model") {
                transform(line.substr(6));
            }
            else if (word=="sphere") {
                makeSphere(line.substr(7));
            }
            else if (word=="light") {
                makeLight(line.substr(6));
            }
            else if (word=="ambient") {
                setAmbient(line.substr(8));
            }
            else if (word=="eye") {
                eye = line.substr(4);
            }
            else if (word=="look") {
                look = line.substr(5);
            }
            else if (word=="up") {
                up = line.substr(3);
            }
            else if (word=="d") {
                d = line.substr(2);
            }
            else if (word=="bounds") {
                bounds = line.substr(7);
            }
            else if (word=="res") {
                res = line.substr(4);
            }
            else if (word=="model") {
                transform(line.substr(6));
            }
            else if (word=="recursionlevel") {
                istringstream lineStream(line.substr(15));
                lineStream>>recursionLevel;
            }
        }
    }
    
    string camSpecs = eye + ' ' + look + ' ' + up + ' ' + d + ' '+ bounds + ' ' + res;
    
    makeCamera(camSpecs);
    
    /*
    
    for (Camera c:cameras) {
        cout << c.getPrettyString() << '\n';
    }
    
    cout << "Ambient: ["<< ambient[0]<<' '<< ambient[1]<<' '<< ambient[2]<<"]\n";
    
    for (Sphere s:spheres) {
        cout << s.getPrettyString() << '\n';
    }
    
    for (Light l:lights) {
        cout << l.getPrettyString() << '\n';
    }
    */
    
    rayTrace(argv[2], recursionLevel, multithread);
    
    cout<<"Took "<<((clock() - initial) / CLOCKS_PER_SEC)<<" seconds to finish.\n";
    return 0;
}
