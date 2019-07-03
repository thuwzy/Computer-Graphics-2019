#pragma once
#include "vec.h"

struct Ray {
    Vec o, d;
    Vec get_point(double t){
        return o+d*t;
    } 
    Ray(Vec o_, Vec d_) : o(o_), d(d_) {} 
}; 