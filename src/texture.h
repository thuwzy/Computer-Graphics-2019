#pragma once
#include "vec.h"
#include <iostream>
#include <stdlib.h>

#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"

class Texture{
public:
	unsigned char* img;
	std::string filename;
	unsigned short xi[3];
	int w,h,n;//n = number of channel
	Texture(){
		img=NULL;
		xi[0]=1,xi[1]=123,xi[2]=6234;
	}
	Texture(std::string _filename):filename(_filename){
		xi[0]=1,xi[1]=123,xi[2]=6234;
		if (filename!=""){
			img=stbi_load(filename.c_str(),&w,&h,&n,0);
			std::cout<<"Successively load "<<filename<<"!\n";
		}
		else{
			img=NULL;
		}
	}
	Vec get_color(double a,double b){
		if (img==NULL){
			std::cout<<"Using non-exist texture!\n";
			return Vec();
		}
		int _w=(int(a*w)%w+w)%w,_h=(int(b*h)%h+h)%h;
		int idx=_h*w*n+_w*n;
		return Vec(double((int)img[idx+0]/255.),double((int)img[idx+1]/255.),double((int)img[idx+2]/255.));
	}

	Vec get_texture(Vec pos,Vec centre=Vec()){
		if (filename=="texture/wood.jpg"){
			double x,z;
			double scaling=1.2;
			if (fmod(pos.x,8*scaling)<4*scaling)
				x=fmod(pos.x,4*scaling)/(4.0*scaling),z=fmod(pos.z,30*scaling)/(30.0*scaling);
			else
				x=fmod(pos.x,4*scaling)/(4.0*scaling),z=fmod(pos.z+15*scaling,30.0*scaling)/(30.0*scaling);
			return get_color(x,z);
		}
		else if (filename=="texture/vase.png"){
			return get_color(pos.x,pos.y);
		}
		return get_color(pos.x,1-pos.y);
		//return get_color(pos.x,pos.y);
		std::cout<<"Using non-exist texture!\n";
		return Vec();
	}

};





