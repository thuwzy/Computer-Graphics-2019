#pragma once
#include "utils.h"
#include "vec.h"
#include "ray.h"
#include "texture.h"
class Object{
public:
};


class Sphere:public Object {
public: 
    double rad;       // radius 
    Vec p, e, c;      // position, emission, color 
    Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive) 
	Texture t;
    Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_,std::string filename=""): 
    rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {
		t=Texture(filename);
	} 
    double intersect(const Ray &r){ // returns distance, 0 if nohit 
        Vec op = p-r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0 
        double t, eps=1e-4, b=op.dot(r.d), det=b*b-op.dot(op)+rad*rad; 
        if (det<0) return 0; else det=sqrt(det); 
        return (t=b-det)>eps ? t : ((t=b+det)>eps ? t : 0); 
    } 
};

