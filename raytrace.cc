#include "Object.h"
#include "Ray.h"
#include "Camera.h"
#include "Sphere.h"
#include "Light.h"
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

using namespace std;

vector<Object> objects;

vector<Sphere> spheres;

vector<Light> lights;

vector<Camera> cameras;

map<string,int> objectsProcessed;

vector<double> ambient;

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

double magnitude(vector<double> vect){
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
    
    Sphere newSphere(x, y, z, r, ambientR, ambientG, ambientB, diffuseR, diffuseG, diffuseB, specularR, specularG, specularB, attenuationR, attenuationG, attenuationB);
    spheres.push_back(newSphere);
    

}


void transform(string input) {
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
    
    //Set object name
    lineStream>>objectName;
    
    ifstream checkFileExists(objectName);
    if (!checkFileExists.good()){
        cout<<"Unable to find model file '"<<objectName<<"'\n";
        return;
    }
        
    
    //Make object
    Object object(objectName);
    
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
    cout << fixed << setprecision(2) <<percent  << '%' << '\r' << flush;
}



vector<double> checkSphereCollision(Ray ray) {
    double minDistance = 9999999;
    int hitSphereIndex = -1;
    vector<double> hitPoint;
    for (int i=0; i<spheres.size();i++) {
        vector<double> rayToSphere = vectorSubtraction(spheres[i].point, ray.origin);
        
        double v = dotProduct(rayToSphere, ray.direction);
        double cc = dotProduct(rayToSphere, rayToSphere);
        
        double d = (spheres[i].radius*spheres[i].radius) - (cc - (v*v));
        d = pow(d,0.5);
        double distance = v - d;
        
        if (d > 0 ) {
            //cout << distance << ' ';
            if (distance < minDistance) {
                minDistance = distance;
                hitSphereIndex = i;
                hitPoint = vectorAddition(ray.origin, scalarMult(ray.direction, (v-d)));
            }
        }
        //return 0;
        
        
    }
    if (hitSphereIndex > -1) {
        //cout << "min: " << minDistance << '\n';
        ray.setHit();
        ray.setHitPosition(hitPoint[0],hitPoint[1],hitPoint[2]);
        Sphere hitSphere = spheres[hitSphereIndex];
        
        vector<double> sphereNormal = vectorSubtraction(hitPoint, hitSphere.point);
        sphereNormal = normVect(sphereNormal);
        
            
        vector<double> color = pairwiseProduct(ambient,hitSphere.ambient);
        
        
        /*
        for lt in lights :
            ptL = lt['p']
            emL = lt['e']
            toL = ptL - ptos; toL = toL / toL.norm()
            if (snrm.dot_product(toL) > 0.0) :
                color += mat1['kd'].pairwise_product(emL) * snrm.dot_product(toL)
                */
       for (Light l:lights) {
           vector<double> toL = vectorSubtraction(l.point, hitPoint);
           toL = normVect(toL);
           double dotProd = dotProduct(sphereNormal, toL);
           if (dotProd > 0.0) {
                color = vectorAddition(color, scalarMult(pairwiseProduct(hitSphere.diffuse, l.color), dotProd));   
           }
           
       }
        
        
        
        
        return color;
    }
    return vector<double> {0, 0, 0};
}

string calculateRayPixelColor(Ray ray) {
    
    vector<double> color{0, 0, 0};
    
    color = vectorAddition(color, checkSphereCollision(ray));
    
    
    for (int i=0; i<color.size();i++) {
            if (color[i] < 0)
                color[i]=0;
            else if (color[i]>1)
                color[i]=1;
        }
    
    return to_string((int)(color[0]*255)) + " " + to_string((int)(color[1]*255)) + " " + to_string((int)(color[2]*255)) + " ";
    
}


void rayTrace(string fileName) {
    Camera camera = cameras[0];
    
    double left = camera.bounds[0];
    double right = camera.bounds[1];
    double bottom = camera.bounds[2];
    double top = camera.bounds[3];
    int width = camera.res[0];
    int height = camera.res[1];
    double near = camera.d;
    
    cout<<"Generating rays...\n";
    
    vector<Ray> rays;
    printPercent(0);
    for (int i=0; i<height;i++) {
        for (int j=0; j<width;j++) {
            double orgX;
            double orgY;
            double orgZ;
            double dirX;
            double dirY;
            double dirZ;
            
            
            double px = (j*1.0)/(width-1)*(right-left)+left;
            double py = (i*1.0)/(height-1)*(bottom-top)+top;
            
            vector<double> origin = vectorAddition( vectorAddition(camera.eye, scalarMult(camera.w, near)), vectorAddition(scalarMult(camera.u, px), scalarMult(camera.v, py))   );
            
            
            vector<double> dir = vectorSubtraction(origin, camera.eye);
            dir = normVect(dir);
            
            Ray ray(origin[0],origin[1],origin[2],dir[0],dir[1],dir[2]);
            rays.push_back(ray);
            printPercent(((((i*1.0*width)+j)/(width*height)) * 100.0));
            
        }
        
    }
    cout<<"Ray generation complete.\n";
    /*
    for (int i=0;i<rays.size();i++) {
        cout << rays[i].getString() << '\t';
        if ((i+1)%(width) == 0)
            cout<<'\n';
    }*/
    
    cout<<"Building image file '"+fileName+"'...\n";
    string imageString = "P3\n"s +  to_string(width)+" "s+to_string(height)+"\n255\n"s;
    
    printPercent(0);
    for (int i=0;i<rays.size();i++) {
        //cout << rays[i].getString() << '\t';
        //double hit = checkSphereCollision(rays[i]);
        imageString += calculateRayPixelColor(rays[i]);
            
        if ((i+1)%(width) == 0)
            imageString += '\n';
        printPercent(((i*1.0)/rays.size())*100);
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
    rayTrace(argv[2]);
    
}
