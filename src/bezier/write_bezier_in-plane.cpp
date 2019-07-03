#include <math.h>
#include <stdlib.h>
#include <stdio.h>
//#include "rand.hpp"
#include <iostream>
#include <fstream>
#include <stdlib.h>
using namespace std;

const int max_size=100;


class Vec
{
public:
	double x,y,z;
	Vec(){x=y=z=0;}
	Vec(double _x,double _y=0,double _z=0){x=_x;y=_y;z=_z;}
	Vec operator+ (const Vec &b) const{return Vec(x+b.x,y+b.y,z+b.z);}
	Vec operator- (const Vec &b) const{return Vec(x-b.x,y-b.y,z-b.z);}
	//void operator=(const Vec &b) const{x=b.x,y=b.y,z=b.z;}
	Vec operator* (const double b) const{return Vec(x*b,y*b,z*b);}
	Vec mult(const Vec &b) const{return Vec(x*b.x,y*b.y,z*b.z);}
	Vec& norm(){return *this=*this * (1/sqrt(x*x+y*y+z*z));}
	double dot(const Vec& b) const{return x*b.x+y*b.y+z*b.z;}
	Vec operator%(Vec b)const{return Vec(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);} 
	friend ostream& operator<<(ostream& out,const Vec s){
    out << "Vec: x="<<s.x<<" y="<<s.y<<" z="<<s.z<<endl;
    return out; }
};
const int w=400,h=400;
Vec image[w][h];

Vec P[10][10];
void draw_bezier(int n)
{
	int fine=10000;
	double t=0;
	for (int i=0;i<=fine;i++)
	{
		t=double(i)/double(fine);
		for (int j=0;j<n;j++)
		{
			for (int k=0;k<n-j;k++)
			{
				P[j+1][k]=P[j][k]*t+P[j][k+1]*(1-t);
			}
		}
		image[int(P[n][0].y+0.5)][int(P[n][0].x+0.5)]=Vec(255,255,255);
		//cout<<int(P[n][0].x+0.5)<<" | "<<int(P[n][0].y+0.5)<<endl;
	}
}


void drawline(Vec a,Vec b){
	int x1=a.x,x2=b.x,y1=a.y,y2=b.y;
	double t=double(y2-y1+.0)/double(x2-x1+.0);
	double now=y1;
	for (int i=x1;i<=x2;i++){
		image[int(now+0.5)][i]=Vec(255,255,255);
		now+=t;
	}

}
void draw_circle(Vec p,int rad)
{
	for (int x=0;x<=rad;x++)
	{
		int y=int(sqrt(rad*rad-x*x)+0.5);
		image[int(0.5+p.y+y)][int(0.5+p.x+x)]=Vec(255,255,255);
		image[int(0.5+p.y+y)][int(0.5+p.x-x)]=Vec(255,255,255);
		image[int(0.5+p.y-y)][int(0.5+p.x+x)]=Vec(255,255,255);
		image[int(0.5+p.y-y)][int(0.5+p.x-x)]=Vec(255,255,255);
	}
}
int main(int argc, char *argv[])
{ 

  for (int i=0; i<h; i++)
 	for (int j=0;j<w;j++)
 		image[i][j]=Vec(0,0,0);

 	P[0][0]=Vec(0,0)+Vec(100,100);
	P[0][1]=Vec(27,37)+Vec(100,100);
	P[0][2]=Vec(80,75)+Vec(100,100);
	P[0][3]=Vec(80.1,0)+Vec(100,100);
	draw_bezier(3);
	//draw_circle(Vec(100+70,150),82);
	//drawline(Vec(10,150,0),Vec(390,150,0));
	FILE *f = fopen("out_bezier.ppm", "w");         // Write image to PPM file. 
	fprintf(f, "P3\n%d %d\n%d\n", w, h, 255); 
	for (int i=h-1; i>=0; i--)
		for (int j=0;j<w;j++)
			fprintf(f,"%d %d %d ", 255-int(image[i][j].x+0.5), 255-int(image[i][j].y+0.5), 255-int(image[i][j].z+0.5)); 
 
}



