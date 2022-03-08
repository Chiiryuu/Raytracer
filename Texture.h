#ifndef Texture_INCLUDED
#define Texture_INCLUDED
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;

class Texture{
  public:
  
  //actual data container
    vector<double> values;
    
    int width=0;
    int height=0;
    int tiling=1;
    bool bump;
    
    Texture() {
    }
  
	//Explicit constructor
	Texture(string fileName, int tiles=1, bool isBump=false) {
        
        /*
        P3
# CREATOR: GIMP PNM Filter Version 1.1
60 52
255
226
*/      bump = isBump;    
        
    
        if (fileName == "") {
            //cout<<"empty file"<<'\n';   
        }
        
        else {
            //cout<<"Raw tiling: "<<tiles<<'\n';
            if (tiles==0)
              tiles = 1;
            tiling = tiles;
            fileName = "materials/"+fileName;
            string line;
            ifstream file (fileName);
            if (!file) {
              cout<<"Unable to find texture file "<<fileName<<'\n';
            }
            else if (file.is_open())
            {
                bool foundWidth = false;
                bool foundMax = false;
                int i = 0;
                while ( getline (file,line) )
                    {
                        if (line.length() == 0);
                        else if (line == "P3");
                        else if (line.at(0) == '#');
                        else if (!foundMax && line=="255") {
                            foundMax = true;
                        }
                        else {
                          
                          istringstream lineStream(line);
                          
                          
                          if (!foundWidth) {
                            foundWidth = true;
                            
                            lineStream>>width;
                            lineStream>>height;
                              
                          }
                          else {
                            if (bump && i==0) {
                              double tempDouble;
                              lineStream>>tempDouble;
                              tempDouble = tempDouble / 255.0;
                              values.push_back(tempDouble);
                              i++;
                              if (i==3)
                                i=0;
                              
                            }
                            else {
                              double tempDouble;
                              lineStream>>tempDouble;
                              tempDouble = tempDouble / 255.0;
                              values.push_back(tempDouble);
                            }
                            //cout << line << '\n';
                          }
                        }
                        
                        
                    }
                file.close();
            }

      
    
        }
      //cout<<"Width: "<<width<<", Height: "<<height<<'\n';
    
      //cout<<"height: "<<height<<", width: "<<width<<'\n';
    
      /*
      for (int i=0;i<width*height;i++) {
        //cout << tempDouble << '\n';
        int row = i % width;
        int col = i / width;
        vector<double> color = get(row,col);
        cout << color[0] << ' ' << color[1] << ' ' << color[2] << ' ' << '\n';
      }
      */
    
    }
  
	

	//Print Light point for debug purposes
  
    vector<double> get(double u, double v) {
        
        if (u<0)
          u=0;
        if (v<0)
          v=0;
        //cout<<"tiling: "<<tiling<<'\n';
      
        if (tiling > 1) {
          u = u * tiling;
          while (u>1)
            u = u - 1;
          
          v = v * tiling;
          while (v>1)
            v = v - 1;
        }
        
      
        int row = round(u * width);
        int col = round((1-v) * height);
        
        if (row > width)
          row = width;
        if (col > height)
          col = height;
      
        //cout<<"Width: "<<width<<", Height: "<<height<<'\n';  
      
        //cout<<"Row: "<<row<<", Col: "<<col<<'\n';
      
        int i = row + col*width;
      
        if (bump) {
          if (i>= values.size())
            return vector<double> {0.5,0.5,0.5};
          return vector<double> {values[i]};
        }
        else {
          i = i * 3;
          if (i>= values.size())
            return vector<double> {1,1,1};
          return vector<double> {values[i], values[i+1], values[i+2]};
        }
        
    }
  
  
    vector<double> get(int row, int col) {
        int i = row + col*width;
        i = i * 3;
        if (i>= values.size())
          return vector<double> {1,1,1};
        return vector<double> {values[i], values[i+1], values[i+2]};
        
    }
	
  private:
	
    };
#endif 
