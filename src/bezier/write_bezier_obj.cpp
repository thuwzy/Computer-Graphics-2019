#include <math.h>
#include <stdlib.h>
#include <stdio.h>
//#include "rand.hpp"
#include <iostream>
#include <fstream>
#include <stdlib.h>
using namespace std;

const int max_size=100;
int n=3;

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
int theta_fine=90;
int t_fine=110;


Vec P[100][100];
void draw_bezier()
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


Vec get_pos(double t,double theta){
	for (int j=0;j<n;j++){
		for (int k=0;k<n-j;k++){
			P[j+1][k]=P[j][k]*t+P[j][k+1]*(1-t);
		}
	}
	return Vec(P[n][0].x*cos(double(theta)*2*M_PI),P[n][0].y,P[n][0].x*sin(double(theta)*2*M_PI));
}


Vec point[200][200];
int f[100000][3];
int id(int x,int y){
	return x*(t_fine+1)+y+1;
}

void print_point(){
	FILE *f = fopen("my_bezier.obj", "w");         // Write image to PPM file. 
	for (int i=0;i<=theta_fine;i++){
		for (int j=0;j<=t_fine;j++){
			point[i][j]=get_pos(double(j)/double(t_fine),double(i)/double(theta_fine));
			fprintf(f, "v %lf %lf %lf\n", point[i][j].x, point[i][j].y, point[i][j].z);
		}
	}
	for (int i=0;i<=theta_fine;i++){
		for (int j=0;j<=t_fine;j++){
			point[i][j]=get_pos(double(j)/double(t_fine),double(i)/double(theta_fine));
			fprintf(f, "vt %lf %lf\n", double(i)/double(theta_fine),double(j)/double(t_fine));
		}
	}
	for (int i=0;i<=theta_fine-1;i++){
		for (int j=1;j<=t_fine-1;j++){
			fprintf(f, "f %d %d %d\n", id(i,j+1), id(i,j), id(i+1,j));
			fprintf(f, "f %d %d %d\n", id(i,j), id(i+1,j-1), id(i+1,j));
		}
	}
}

int main(int argc, char *argv[])
{ 

  for (int i=0; i<h; i++)
 	for (int j=0;j<w;j++)
 		image[i][j]=Vec(0,0,0);
/*
 	P[0][16]=Vec(0,80.1);
	P[0][15]=Vec(30,80);
	P[0][14]=Vec(30,80);
	P[0][13]=Vec(30,80);
	P[0][12]=Vec(30,80);
	P[0][11]=Vec(10,70.1);
	P[0][10]=Vec(10,70.1);
	P[0][9]=Vec(20,60.1);
	P[0][8]=Vec(25,40.1);
	P[0][7]=Vec(30,30.1);
	P[0][6]=Vec(30,30.1);
	P[0][5]=Vec(30,20.1);
	P[0][4]=Vec(30,10.1);
	P[0][3]=Vec(30,0.1);
	P[0][2]=Vec(30,0.1);
	P[0][1]=Vec(30,0.1);
	P[0][0]=Vec(0,0);
*/		
	P[0][3]=Vec(0,80.1);
	P[0][2]=Vec(75,80);
	P[0][1]=Vec(37,27);
	P[0][0]=Vec(0,0);


	draw_bezier();
	print_point();
	FILE *f = fopen("out_bezier.ppm", "w");         // Write image to PPM file. 
	fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
	for (int i=h-1; i>=0; i--)
		for (int j=0;j<w;j++)
			fprintf(f,"%d %d %d ", 255-int(image[i][j].x+0.5), 255-int(image[i][j].y+0.5), 255-int(image[i][j].z+0.5)); 
 
}


