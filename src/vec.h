#pragma once
#include "utils.h"

struct Vec {        // Usage: time ./smallpt  5000 && xv image.ppm 
    double x, y, z;                  // position, also color (r,g,b) 
    Vec(double x_=0, double y_=0, double z_=0){ x=x_; y=y_; z=z_; } 
    Vec operator-()const {return Vec(-x,-y,-z);}
    Vec operator+(const Vec &b) const { return Vec(x+b.x,y+b.y,z+b.z); } 
    Vec operator-(const Vec &b) const { return Vec(x-b.x,y-b.y,z-b.z); } 
    Vec operator*(double b) const { return Vec(x*b,y*b,z*b); }
    Vec operator+(double b) const { return Vec(x+b,y+b,z+b); }
    Vec operator/(double b) const { return Vec(x/b,y/b,z/b); }
    Vec operator-(double b) const { return Vec(x-b,y-b,z-b); }
	bool operator==(const Vec&a) const {return x==a.x && y==a.y && z==a.z;}
	bool operator!=(const Vec&a) const {return x!=a.x || y!=a.y || z!=a.z;}
	Vec&operator+=(const Vec&a) {return *this = *this + a;}
	Vec&operator-=(const Vec&a) {return *this = *this - a;}
	Vec&operator+=(double p) {return *this = *this + p;}
	Vec&operator-=(double p) {return *this = *this - p;}
	Vec&operator*=(double p) {return *this = *this * p;}
	Vec&operator/=(double p) {return *this = *this / p;} 
    Vec mult(const Vec &b) const { return Vec(x*b.x,y*b.y,z*b.z); } 
    Vec& norm() {return *this = *this * (1/sqrt(x*x+y*y+z*z));} 
	Vec max(const Vec&a) const {return Vec(std::max(x,a.x), std::max(y,a.y), std::max(z,a.z));}
	Vec min(const Vec&a) const {return Vec(std::min(x,a.x), std::min(y,a.y), std::min(z,a.z));}
    double dot(const Vec &b) const { return x*b.x+y*b.y+z*b.z; } // cross: 
    Vec operator%(Vec b)const{return Vec(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);} 
	friend std::ostream& operator<<(std::ostream& out,const Vec s){
    out << "Vec: x="<<s.x<<" y="<<s.y<<" z="<<s.z<<std::endl;
    return out; }
}; 