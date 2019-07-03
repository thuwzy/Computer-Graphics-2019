#pragma once
#include "ppm.h"
#include "scene.h"
#include <math.h>
#include <cmath>
#include <cstring>
#include <omp.h>
#define _Depth
void ProgressivePhotonMapping::RayTrace(const Ray &r, int depth, unsigned short *Xi,Vec color,int p_x,int p_y){

	double t;
	int id=0;
	Shape shape;
	Vec x,n,nl,ref_color,emit;
	Refl_t refl;
	intersect(r,t,id,shape,ref_color,n);
	x=r.o+r.d*t;//光线击中点
	if (shape==EMPTY) return ;
	if (shape==SPHERE)
	{
		Sphere &obj=spheres[id];
		n=(x-obj.p).norm();//光线击中点的外法线
		nl=n.dot(r.d)<0?n:n*-1;//nl与光线方向相反
		emit=obj.e;
		refl=obj.refl;
	}
    if (shape==KDTREE){
 		KDTree &obj=kdtrees[id];
		//n=obj.n_ans;//光线击中点外法线
		nl=n.dot(r.d)<0?n:n*-1;//nl与光线方向相反

		//n=nl;

		emit=Vec();
		refl=obj.refl;       
    }
    if (shape==BEZIER){
 		Bezier &obj=beziers[id];
		//n=obj.n_ans;//光线击中点外法线
		nl=n.dot(r.d)<0?n:n*-1;//nl与光线方向相反
		emit=Vec();
		refl=obj.refl; 
		if (erand48(Xi)<obj.refl_rate) refl=SPEC;   
    }

	double p;
	if (ref_color.x>ref_color.y&&ref_color.x>ref_color.z) p=ref_color.x;
	else if (ref_color.y>ref_color.z) p=ref_color.y;
	else p=ref_color.z;

	if (++depth>5)
	{
		if (erand48(Xi)<p) ref_color=ref_color*(1/p);
		else return;
	}
	if (refl==DIFF)
	{

        tmp_vector_vp[omp_get_thread_num()].push_back(ViewPoint(x,nl,r.d,ref_color.mult(color),all_energy,p_x,p_y));
        return; 
	}
	else if (refl==SPEC)
	{
		RayTrace(Ray(x,r.d-n*2*n.dot(r.d)),depth,Xi,color.mult(ref_color),p_x,p_y);
        return;
	}
	Ray reflRay(x,r.d-n*2*n.dot(r.d));
	bool into=n.dot(nl)>0;
	double nc=1,nt=1.5,nnt=into?nc/nt:nt/nc,ddn=r.d.dot(nl),cos2t;
	cos2t=1-nnt*nnt*(1-ddn*ddn);
	if (cos2t<0){
		RayTrace(reflRay,depth,Xi,color.mult(ref_color),p_x,p_y);
		return;
	}
		
	Vec tdir=(r.d*nnt-n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm();
	double a=nt-nc,b=nt+nc,R0=a*a/(b*b),c=1-(into?-ddn:tdir.dot(n));
	double Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=0.25+0.5*Re,RP=Re/P,TP=Tr/(1-P);
    if (depth>2){
        if (erand48(Xi)<P){
            RayTrace(reflRay,depth,Xi,color.mult(ref_color)*RP,p_x,p_y);
        }
        else{
            RayTrace(Ray(x,tdir),depth,Xi,color.mult(ref_color)*TP,p_x,p_y);
        }
    }
    else{
        RayTrace(reflRay,depth,Xi,color.mult(ref_color)*Re,p_x,p_y);
        RayTrace(Ray(x,tdir),depth,Xi,color.mult(ref_color)*Tr,p_x,p_y);
        
    }    
}

void ProgressivePhotonMapping::creating_viewpoint(){
	for (int __i=0;__i<cpu_num;__i++){
		tmp_vector_vp[__i].clear();
	}
    printf("Creating ViewPoint...\n");
    vector_vp.clear();
    Ray cam(Vec(50,52,295.6), Vec(0,-0.042612,-1).norm()); // cam pos, dir 
    Vec cx=Vec(width*.5135/height), cy=(cx%cam.d).norm()*.5135, r;
 /* 
#pragma omp parallel for schedule(dynamic, 1) private(r) 
    for (int y=0; y<height; y++){                       // Loop over image rows 
		fprintf(stderr,"\rCreating (%d spp) %5.2f%%\n",samps*4,100.*y/(height-1)); 
		for (unsigned short x=0; x<width; x++)   // Loop cols 
			//for (int sy=0; sy<2; sy++)     // 2x2 subpixel rows 
				//for (int sx=0; sx<2; sx++)

				{        // 2x2 subpixel cols 
					for (int s=0; s<samps; s++){
						int sy=erand48(Xi)<0.5?0:1;
						int sx=erand48(Xi)<0.5?0:1; 
						double r1=2*erand48(Xi), dx=r1<1 ? sqrt(r1)-1: 1-sqrt(2-r1); 
						double r2=2*erand48(Xi), dy=r2<1 ? sqrt(r2)-1: 1-sqrt(2-r2); 
						Vec d = cx*( ( (sx+.5 + dx)/2 + x)/width - .5) + 
								cy*( ( (sy+.5 + dy)/2 + y)/height - .5) + cam.d; 
						//if (erand48(Xi)<0.25)
						RayTrace(Ray(cam.o+d*140,d.norm()),0,Xi,Vec(1,1,1),x,y); 
					} // Camera rays are pushed ^^^^^ forward to start in interior 

				} 
    }
*/
	double f1=2,aperture=1,f2=220;
#pragma omp parallel for schedule(dynamic, 1) private(r) 
//#pragma omp parallel for num_threads(3)
    for (int y=0; y<height; y++){                       // Loop over image rows 
		fprintf(stderr,"\rCreating (%d spp) %5.2f%%\n",samps*4,100.*y/(height-1)); 
		for (unsigned short x=0; x<width; x++)   // Loop cols 
			//for (int sy=0, i=(height-y-1)*width+x; sy<2; sy++)     // 2x2 subpixel rows 
				//for (int sx=0; sx<2; sx++, r=Vec())
				{        // 2x2 subpixel cols 
					for (int s=0; s<samps; s++){ 
						int sy=erand48(Xi)<0.5?0:1;
						int sx=erand48(Xi)<0.5?0:1;
						double r1=2*erand48(Xi), dx=r1<1 ? sqrt(r1)-1: 1-sqrt(2-r1); 
						double r2=2*erand48(Xi), dy=r2<1 ? sqrt(r2)-1: 1-sqrt(2-r2); 
						Vec d = cx*( ( (sx+.5 + dx)/2 + x)/width - .5) + 
								cy*( ( (sy+.5 + dy)/2 + y)/height - .5) + cam.d; 

						double theta=erand48(Xi)*2*M_PI;
						Vec point1=cam.o+cam.d*f1+(cx*cos(theta)*aperture+cy*sin(theta)*aperture)*(erand48(Xi)-0.5);
						RayTrace(Ray(point1+(cam.o+d*f2-point1).norm()*123,(cam.o+d*f2-point1).norm()),0,Xi,Vec(1,1,1),x,y); 
					} // Camera rays are pushed ^^^^^ forward to start in interior 

				} 
    }

	for (int __i=0;__i<cpu_num;__i++){
		for (int j=0;j<tmp_vector_vp[__i].size();j++){
			vector_vp.push_back(tmp_vector_vp[__i][j]);
		}
		tmp_vector_vp[__i].clear();
	}
    printf("Created %d veiw points.\n",(int)vector_vp.size()); 
}

void ProgressivePhotonMapping::PhotonTracing(const Ray &r, int depth, unsigned short *Xi,Vec color){
	double t;
	int id=0;
	Shape shape;
	Vec x,n,nl,ref_color,emit;
	Refl_t refl;
	intersect(r,t,id,shape,ref_color,n);
	x=r.o+r.d*t;//光线击中点
	if (shape==EMPTY) return ;
	if (shape==SPHERE)
	{
		Sphere &obj=spheres[id];
		n=(x-obj.p).norm();//光线击中点的外法线
		nl=n.dot(r.d)<0?n:n*-1;//nl与光线方向相反
		emit=obj.e;
		refl=obj.refl;
	}
    if (shape==KDTREE){
 		KDTree &obj=kdtrees[id];
		//n=obj.n_ans;//光线击中点外法线
		nl=n.dot(r.d)<0?n:n*-1;//nl与光线方向相反
		//n=nl;
		emit=Vec();
		refl=obj.refl;       
    }
    if (shape==BEZIER){
 		Bezier &obj=beziers[id];
		//n=obj.n_ans;//光线击中点外法线
		nl=n.dot(r.d)<0?n:n*-1;//nl与光线方向相反
		emit=Vec();
		refl=obj.refl;
		if (erand48(Xi)<obj.refl_rate) refl=SPEC;         
    }


	double p;
	if (ref_color.x>ref_color.y&&ref_color.x>ref_color.z) p=ref_color.x;
	else if (ref_color.y>ref_color.z) p=ref_color.y;
	else p=ref_color.z;

	if (++depth>5)
	{
		if (erand48(Xi)<p) ref_color=ref_color*(1/p);
		else return;
	}
	if (refl==DIFF)
	{
		double r1=2*M_PI*erand48(Xi),r2=erand48(Xi),r2s=sqrt(r2);
		Vec w=nl,u=((fabs(w.x)>0.1?Vec(0,1,0):Vec(1,0,0))%w).norm(),v=w%u;
		Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm();
		PhotonTracing(Ray(x,d),depth,Xi,color.mult(ref_color));


        std::vector<ViewPoint> result;
        query(x,result);
        for (int i=0;i<result.size();i++){
            buf_img[omp_get_thread_num()][result[i].pix_x+result[i].pix_y*width]+=
                result[i].color.mult(color)*all_energy*(1./(double)samps);//*d.dot(nl);
        }


        return; 
	}
	else if (refl==SPEC)
	{
		PhotonTracing(Ray(x,r.d-n*2*n.dot(r.d)),depth,Xi,color.mult(ref_color));
        return;
	}
	Ray reflRay(x,r.d-n*2*n.dot(r.d));
	bool into=n.dot(nl)>0;
	double nc=1,nt=1.5,nnt=into?nc/nt:nt/nc,ddn=r.d.dot(nl),cos2t;
	cos2t=1-nnt*nnt*(1-ddn*ddn);
	if (cos2t<0){
		PhotonTracing(reflRay,depth,Xi,color.mult(ref_color));
		return;
	}
		
	Vec tdir=(r.d*nnt-n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm();
	double a=nt-nc,b=nt+nc,R0=a*a/(b*b),c=1-(into?-ddn:tdir.dot(n));
	double Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=0.25+0.5*Re,RP=Re/P,TP=Tr/(1-P);
    if (depth>2){
        if (erand48(Xi)<P){
            PhotonTracing(reflRay,depth,Xi,color.mult(ref_color)*RP);
        }
        else{
            PhotonTracing(Ray(x,tdir),depth,Xi,color.mult(ref_color)*TP);
        }
    }
    else{
        PhotonTracing(reflRay,depth,Xi,color.mult(ref_color)*Re);
        PhotonTracing(Ray(x,tdir),depth,Xi,color.mult(ref_color)*Tr);
        
    }     
}
void ProgressivePhotonMapping::EmitPhoton(){
    
    Vec p=Vec(50,81.5,40.8);
    //Vec d=Vec(erand48(r_seed)-0.5,-10,erand48(r_seed)-0.5).norm();
#pragma omp parallel for schedule(dynamic,int(photon_num/cpu_num))
//#pragma omp parallel for num_threads(3)
    for (int i=0;i<photon_num;i++){
        Vec d=Vec(erand48(r_seed)-0.5,erand48(r_seed)-1,erand48(r_seed)-0.5).norm();
        PhotonTracing(Ray(p+Vec((erand48(r_seed)-0.5),0,(erand48(r_seed)-0.5))*1.6,d),0,r_seed,Vec(1,1,1));
        if (i%10000==0) printf("Emitting %d th photon\n",i);
    }
}




void ProgressivePhotonMapping::run(){
    all_radius=1.58;
    all_energy=0.2/(all_radius*all_radius);
    radius_decay=0.9991;
    photon_num=100000;
    round_num=502;

	creating_viewpoint();
    build_tree();   
    result_img=new Vec[width*height+10];
    for (int i=0;i<=width*height;i++) result_img[i]=Vec();
    
    for (int i=0;i<round_num;i++){
        if (i>0){
            creating_viewpoint();
            delete[] tree;
            build_tree();           
        }
        printf("Emitting Photon...Round %d\n",i);
		for (int __i=0;__i<cpu_num;__i++){
        	buf_img[__i]=new Vec[width*height+10];
        	for (int _i=0;_i<=width*height;_i++) buf_img[__i][_i]=Vec();
		}


        EmitPhoton();

		for (int __i=0;__i<cpu_num;__i++){
			for (int _i=0;_i<=width*height;_i++){
				result_img[_i]+=buf_img[__i][_i]*(1./round_num);
			}
		}

		for (int __i=0;__i<cpu_num;__i++)
        	delete[] buf_img[__i];
        all_radius*=radius_decay;
        all_energy/=(radius_decay);
        if (i<30||i%50==0){
            std::string s1 = "result12/image";
            std::string s2 = ".ppm";
            std::ostringstream oss;
            oss << s1 << i << s2;
            FILE *f = fopen(oss.str().c_str(), "w");         // Write image to PPM file. 
            fprintf(f, "P3\n%d %d\n%d\n", width, height, 255);
			for (int _i=height-1;_i>=0;_i--){
				for (int j=0;j<width;j++){
					int pos=_i*width+j;
					fprintf(f,"%d %d %d ",
						toInt(result_img[pos].x*(round_num/(i+1))), toInt(result_img[pos].y*(round_num/(i+1))), toInt(result_img[pos].z*(round_num/(i+1))));
				}
			}         
        }

    }
}
