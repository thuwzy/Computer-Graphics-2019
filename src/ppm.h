#pragma once
#include <math.h>  
#include <stdlib.h>
#include <stdio.h>  
#include "rand.hpp"
#include <iostream>
#include "vec.h"
#include "ray.h"
#include "objects.h"
#include <omp.h>

inline bool intersect(const Ray &r,double &t,int &id,Shape& shape,Vec &ref_color ,Vec& _n_ans);

class ViewPoint{
public:
    Vec position,n,w;//w is the direction
    int pix_x,pix_y;
    int N;//number of photons
    double wgt;//weight to pixel
    Vec color;
    double energy;
    ViewPoint(){
        N=0;
    }
    ViewPoint(Vec _p,Vec _n,Vec _w,Vec _c,double _e,int _px,int _py)
        :position(_p),n(_n),w(_w),color(_c),energy(_e),pix_x(_px),pix_y(_py){
        N=0;
    }

};

class ProgressivePhotonMapping{
public:
    static int D;
    struct KDTreeNode{
        ViewPoint view_point;
        Vec m[2];//min and max
        int son[2];
		bool operator<(const KDTreeNode&a) const {
			if (D == 0)
				return view_point.position.x < a.view_point.position.x;
			else if (D == 1)
				return view_point.position.y < a.view_point.position.y;
			else
				return view_point.position.z < a.view_point.position.z;
		}
        KDTreeNode(){
            m[0]=Vec(INF,INF,INF);
            m[1]=Vec(-INF,-INF,-INF);
            son[0]=son[1]=-1;
        }
        void expand(){
            m[0]=view_point.position;
            m[1]=view_point.position;
        }
        bool intersect_box(Vec x){
            double dx=std::max(std::max(0.,m[0].x-x.x),x.x-m[1].x);
            double dy=std::max(std::max(0.,m[0].y-x.y),x.y-m[1].y);
            double dz=std::max(std::max(0.,m[0].z-x.z),x.z-m[1].z);
            if (dx*dx+dy*dy+dz*dz<all_radius*all_radius) return 1;
            return 0;            
        }
    };
    KDTreeNode* tree;
    int photon_num;
    int round_num;
	//omp_lock_t mylock;
    static double all_radius;
    static double all_energy;
    double radius_decay;
    int kdt_root;
    unsigned short Xi[3];
    Vec* buf_img[100];
    Vec* result_img;
	int cpu_num;

    int width,height;
    int samps;
    ProgressivePhotonMapping(int _w=400,int _h=400,int _s=32)
        :width(_w),height(_h),samps(_s){
        r_seed[0]=11;r_seed[1]=22;r_seed[2]=233;
        Xi[0]=0,Xi[1]=1,Xi[2]=332;
		cpu_num=omp_get_num_procs();
		std::cout<<cpu_num<<" cpus in all!\n";
    }
    unsigned short r_seed[3];
    std::vector<ViewPoint> vector_vp;
	std::vector<ViewPoint> tmp_vector_vp[100];
    void RayTrace(const Ray &r, int depth, unsigned short *Xi,Vec color,int p_x,int p_y);
    void run();
    void creating_viewpoint();
    void EmitPhoton();
    void PhotonTracing(const Ray &r, int depth, unsigned short *Xi,Vec color);
    double dist2(int a,int b){
        return (tree[a].view_point.position.x-tree[b].view_point.position.x)*(tree[a].view_point.position.x-tree[b].view_point.position.x)
              +(tree[a].view_point.position.y-tree[b].view_point.position.y)*(tree[a].view_point.position.y-tree[b].view_point.position.y)
              +(tree[a].view_point.position.z-tree[b].view_point.position.z)*(tree[a].view_point.position.z-tree[b].view_point.position.z);
    }
    double dist2(Vec a,Vec b){
        return (a.x-b.x)*(a.x-b.x)+
               (a.y-b.y)*(a.y-b.y)+
               (a.z-b.z)*(a.z-b.z);
    }
    void expand_bynode(int a,int b){
        tree[a].m[0]=tree[a].m[0].min(tree[b].m[0]);
        tree[a].m[1]=tree[a].m[1].max(tree[b].m[1]);
    }
    int _build(int l,int r,int d){
        D=d;
        int subroot=(l+r)>>1;
        std::nth_element(tree+l,tree+subroot,tree+r+1);
        tree[subroot].expand();
        if (l<subroot) tree[subroot].son[0]=_build(l,subroot-1,(d+1)%3),expand_bynode(subroot,tree[subroot].son[0]);
        if (subroot<r) tree[subroot].son[1]=_build(subroot+1,r,(d+1)%3),expand_bynode(subroot,tree[subroot].son[1]);
        return subroot;
    }
    void build_tree(){
        printf("Building KDTree...\n");
        int v_size=vector_vp.size();
        tree=new KDTreeNode[v_size+10];
        for (int i=1;i<=v_size;i++){
            tree[i].view_point=vector_vp[i-1];
        }
        kdt_root=_build(1,v_size,0);
		vector_vp.clear();
    }


    void _query(Vec x,int rt,std::vector<ViewPoint> &result){
        if (!tree[rt].intersect_box(x)) return;
        if (dist2(x,tree[rt].view_point.position)<all_radius*all_radius){
            result.push_back(tree[rt].view_point);
        }
        if (tree[rt].son[0]!=-1) _query(x,tree[rt].son[0],result);
        if (tree[rt].son[1]!=-1) _query(x,tree[rt].son[1],result);
    }
    void query(Vec x,std::vector<ViewPoint> &result){
        result.clear();
        _query(x,kdt_root,result);        
    }
};

int ProgressivePhotonMapping::D = 0;
double ProgressivePhotonMapping::all_radius=1;
double ProgressivePhotonMapping::all_energy=1;










