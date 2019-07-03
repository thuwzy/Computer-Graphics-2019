#pragma once
#include <bits/stdc++.h>

#define M_PI 3.14159265358979323846
enum Refl_t { DIFF, SPEC, REFR };  // material types, used in radiance() 
enum Shape {EMPTY,SPHERE,KDTREE,BEZIER};
const double eps=1e-6;
const double INF=1<<20;



inline double clamp(double x){ return x<0 ? 0 : x>1 ? 1 : x; } 
inline int toInt(double x){ return int(pow(clamp(x),1/2.2)*255+.5); }