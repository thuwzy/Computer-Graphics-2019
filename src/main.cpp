#include <math.h>   
#include <stdlib.h> 
#include <stdio.h>  
#include "rand.hpp"
#include <iostream>
#include "vec.h"
#include "ray.h"
#include "objects.h"
#include "KDTree.h"
#include "utils.h"
#include "scene.h"
#include "ppm.cpp"
using std::cout;

unsigned short Xx[3]={12,44,88};

inline bool intersect(const Ray &r,double &t,int &id,Shape& shape,Vec &ref_color ,Vec& _n_ans){ 
	shape=EMPTY;
	Vec __n_ans;
	double n=sizeof(spheres)/sizeof(Sphere),d,inf=t=1e20;
	for (int i=int(n);i--;)
		if ((d=spheres[i].intersect(r))&&d<t)
		{
			t=d;
			id=i;
			shape=SPHERE;
			ref_color=spheres[i].c;
			if (spheres[i].t.img!=NULL){
				ref_color=spheres[i].t.get_texture(r.o+r.d*t);
			}
		}
	n=sizeof(kdtrees)/sizeof(KDTree);
	for (int i=n;i--;){
		Ray _r=r;
		double ord_x,ord_y;
		int group;
		if ((d=kdtrees[i].intersect(_r,__n_ans,ord_x,ord_y,group))&&d<t)
		{
			Vec _ref_color=Vec(0.99,0.99,0.99);
			if (group!=-1&&kdtrees[i].texture[group].img!=NULL){
				_ref_color=kdtrees[i].texture[group].get_texture(Vec(ord_x,ord_y));
			}
			t=d;
			id=i;
			shape=KDTREE;
			ref_color=_ref_color;
			_n_ans=__n_ans;
		}
	}
	n=sizeof(beziers)/sizeof(Bezier);
	for (int i=n;i--;){
		Ray _r=r;
		double ord_x,ord_y;
		if ((d=beziers[i].b_intersect(_r,__n_ans,ord_x,ord_y))&&d<t)
		{
			t=d;
			id=i;
			shape=BEZIER;
			ref_color=Vec(0.99,0.99,0.99);
			_n_ans=-__n_ans;
			if (beziers[i].texture.img!=NULL){
				ref_color=beziers[i].texture.get_texture(Vec(ord_x,ord_y));
			} 
		}
	}
	return t<inf;
} 

int main(int argc, char *argv[]){
    int w=2880, h=2400, samps = argc==2 ? atoi(argv[1])/4 : 1; // # samples
    w=1080*2,h=768*2;
    samps=1; 

    ProgressivePhotonMapping ppm(w,h,samps);
    ppm.run();
    Vec *c = ppm.result_img;

    FILE *f = fopen("image.ppm", "w");         // Write image to PPM file. 
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255); 
    for (int i=w*h-1; i>=0; i--) 
        fprintf(f,"%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z)); 
}