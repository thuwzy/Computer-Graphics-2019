#pragma once
#include "vec.h"
#include "ray.h"
#include "utils.h"
#include "objects.h"
#include "texture.h"
#include <fstream>
#include <iostream>
#include <algorithm>
class KDTree{
public:
    static int D;
    struct Triangle{//also KDTree Node
        Vec n;//法向
        Vec v[3];//顶点
		Vec vn[3];
		Vec vt[3];
        Vec m[2];//min and max;
        int son[2];//lson and rson
		int group;
        double intersect(const Ray& ray){
            Vec E1=v[0]-v[1];
		    Vec E2=v[0]-v[2];
		    Vec S=v[0]-ray.o;
		    double under=(ray.d%E1).dot(E2);
		    double t=(S%E1).dot(E2)/under;
		    double beta=(ray.d%S).dot(E2)/under;
		    double gamma=(ray.d%E1).dot(S)/under;
		    if (t<=0||beta<0||gamma<0||beta+gamma>1) return 0;
		    return t;
        }
		void get_ord(Vec point,Vec &my_n_ans,double&my_ord_x,double&my_ord_y,int&my_group){
			Vec E1=v[1]-v[0];
		    Vec E0=v[2]-v[0];
		    Vec E2=point-v[0];
			double lambda1=(E1.dot(E1)*E0.dot(E2)-E0.dot(E1)*E1.dot(E2))/(E0.dot(E0)*E1.dot(E1)-E0.dot(E1)*E0.dot(E1));
			double lambda2=(E0.dot(E0)*E1.dot(E2)-E0.dot(E1)*E0.dot(E2))/(E0.dot(E0)*E1.dot(E1)-E0.dot(E1)*E0.dot(E1));
			double lambda3=1-lambda1-lambda2;
			my_n_ans=(vn[0]*lambda3+vn[2]*lambda1+vn[1]*lambda2).norm();
			Vec tmp=vt[0]*lambda3+vt[2]*lambda1+vt[1]*lambda2;
			my_ord_x=tmp.x,my_ord_y=tmp.y;
			my_group=group;
		}
        bool intersect_box(const Ray& ray){
            Ray _r=ray;
            if (_r.d.x==0)_r.d.x+=1e-10;
            if (_r.d.y==0)_r.d.y+=1e-10;
            if (_r.d.z==0)_r.d.z+=1e-10;

            double t_max=INF,t_min=-INF;
            double t1=0,t2=0;

            t1=(m[0].x-_r.o.x)/_r.d.x;
            t2=(m[1].x-_r.o.x)/_r.d.x;
            if (t1>t2){
                double tmp=t2;
                t2=t1;
                t1=tmp;
            }
            if (t_min<t1) t_min=t1;
            if (t2<t_max) t_max=t2;
            
            t1=(m[0].y-_r.o.y)/_r.d.y;
            t2=(m[1].y-_r.o.y)/_r.d.y;
            if (t1>t2){
                double tmp=t2;
                t2=t1;
                t1=tmp;
            }
            if (t_min<t1) t_min=t1;
            if (t2<t_max) t_max=t2;
                    
            t1=(m[0].z-_r.o.z)/_r.d.z;
            t2=(m[1].z-_r.o.z)/_r.d.z;
            if (t1>t2){
                double tmp=t2;
                t2=t1;
                t1=tmp;
            }
            if (t_min<t1) t_min=t1;
            if (t2<t_max) t_max=t2;

            if (t_min<t_max) return 1;
            else return 0;            
        }
        Triangle(){
            m[0]=Vec(INF,INF,INF);
            m[1]=Vec(-INF,-INF,-INF);
            son[0]=son[1]=-1;
        }
        void expand(){
            for (int i=0;i<3;i++){
                m[0]=m[0].min(v[i]);
                m[1]=m[1].max(v[i]);
            }
        }

		bool operator<(const Triangle&a) const {
			if (D == 0)
				return v[0].x < a.v[0].x;
			else if (D == 1)
				return v[0].y < a.v[0].y;
			else
				return v[0].z < a.v[0].z;
		}
    };
    Triangle* meshs;
	Texture texture[15];
    Vec* v_buffer;
    Vec* vn_buffer;
	Vec* vt_buffer;
    int v_num,vt_num,mesh_num,vn_num;
    double query_ans;
    int root;
    static Vec n_ans;
	Refl_t refl;
	double scaling_size;
    KDTree(const char* file_in,Refl_t _refl=SPEC,double _size=300,Vec position=Vec(10,0,10)):refl(_refl),scaling_size(_size){
        query_ans=INF;
		 
		texture[0]=Texture("mesh/room/floor1.jpg");
		texture[1]=Texture("mesh/room/veshalka.jpg");
		texture[2]=Texture("mesh/room/doors1.jpg");
		texture[3]=Texture("mesh/room/Alyscamps_Arles_dos.JPG");
		texture[4]=Texture("mesh/room/belye.jpg");
		texture[5]=Texture("mesh/room/bad1.jpg");
		texture[6]=Texture("mesh/room/rest.jpg");
		texture[7]=Texture("mesh/room/wall_staff.jpg");
		texture[8]=Texture("mesh/room/table.jpg");
		texture[9]=Texture("mesh/room/chair1.jpg");
		texture[10]=Texture("mesh/room/wall1.jpg");
		texture[11]=Texture("mesh/room/win_right.jpg");
		texture[12]=Texture("mesh/room/win_left.jpg");
		texture[13]=Texture("mesh/room/glass_right.png");
		texture[14]=Texture("mesh/room/glass_left.png");
			
        if (meshs!=NULL) delete[] meshs;
        load(file_in);
        transfer(position,scaling_size);
        root=build_tree(1,mesh_num,0);
    }
    void load(const char* file_in){

        std::ifstream in;
        in.open(file_in);
        if (!in.is_open()){ fprintf(stderr, "Wrong when opening file");return;}
        in>>v_num>>vn_num>>vt_num>>mesh_num;

        v_buffer=new Vec[v_num+10];
        vn_buffer=new Vec[vn_num+10];
		vt_buffer=new Vec[vt_num+10];
        meshs=new Triangle[mesh_num+10];

		for (int i=0;i<v_num;i++){
			in>>v_buffer[i].x>>v_buffer[i].y>>v_buffer[i].z;
		}
		fprintf(stderr, "\rLoading %s  1/3 finished\n ",file_in);
		for (int i=0;i<vn_num;i++){
			in>>vn_buffer[i].x>>vn_buffer[i].y>>vn_buffer[i].z;
			vn_buffer[i]=vn_buffer[i].norm();
		}
		for (int i=0;i<vt_num;i++){
			in>>vt_buffer[i].x>>vt_buffer[i].y;

		}
		fprintf(stderr, "\rLoading %s  2/3 finished\n ",file_in);
		for (int i=1;i<=mesh_num;i++){
			int tmp0,tmp1,tmp2;
			in>>meshs[i].group>>tmp0>>tmp1>>tmp2;
			tmp0-=1,tmp1-=1,tmp2-=1;

            meshs[i].v[0]=v_buffer[tmp0];
            meshs[i].v[1]=v_buffer[tmp1];
            meshs[i].v[2]=v_buffer[tmp2];
            meshs[i].n=((meshs[i].v[0]-meshs[i].v[1])%(meshs[i].v[0]-meshs[i].v[2])).norm();
			
			in>>meshs[i].group>>tmp0>>tmp1>>tmp2;
			tmp0-=1,tmp1-=1,tmp2-=1;
            meshs[i].vn[0]=vn_buffer[tmp0];
            meshs[i].vn[1]=vn_buffer[tmp1];
            meshs[i].vn[2]=vn_buffer[tmp2];
			if (vt_num!=0){
				in>>meshs[i].group>>tmp0>>tmp1>>tmp2;
				tmp0-=1,tmp1-=1,tmp2-=1;
				meshs[i].vt[0]=vt_buffer[tmp0];
				meshs[i].vt[1]=vt_buffer[tmp1];
				meshs[i].vt[2]=vt_buffer[tmp2];
			}


		}
        fprintf(stderr, "\rLoading %s  3/3 finished\n ",file_in);
        delete[] v_buffer;
        delete[] vn_buffer;
		delete[] vt_buffer;         
    }

    void transfer(Vec position,double size){
        Vec min_position=Vec(INF,INF,INF);
        for (int i=1;i<mesh_num;i++){
            for (int j=0;j<3;j++){
				
				double tmpx=meshs[i].v[j].x*cos(-M_PI/60)-meshs[i].v[j].z*sin(-M_PI/60);
				double tmpz=meshs[i].v[j].x*sin(-M_PI/60)+meshs[i].v[j].z*cos(-M_PI/60);
				meshs[i].v[j].x=tmpx;
				meshs[i].v[j].z=tmpz;

				double tmpy=meshs[i].v[j].y*cos(-M_PI/60)-meshs[i].v[j].z*sin(-M_PI/60);
				tmpz=meshs[i].v[j].y*sin(-M_PI/60)+meshs[i].v[j].z*cos(-M_PI/60);
				meshs[i].v[j].y=tmpy;
				meshs[i].v[j].z=tmpz;	
					
            }

        }
        for (int i=1;i<mesh_num;i++){
            for (int j=0;j<3;j++){
                min_position=min_position.min(meshs[i].v[j]);
            }
        }
        for (int i=1;i<mesh_num;i++){
            for (int j=0;j<3;j++){
                meshs[i].v[j]=(meshs[i].v[j])*size+position;
            }
        }
    }

    void expand_bynode(int a,int b){
        meshs[a].m[0]=meshs[a].m[0].min(meshs[b].m[0]);
        meshs[a].m[1]=meshs[a].m[1].max(meshs[b].m[1]);
    }
    int build_tree(int l,int r,int d){
        D=d;
        int subroot=(l+r)>>1;
        std::nth_element(meshs+l,meshs+subroot,meshs+r+1);
        meshs[subroot].expand();
        if (l<subroot) meshs[subroot].son[0]=build_tree(l,subroot-1,(d+1)%3),expand_bynode(subroot,meshs[subroot].son[0]);
        if (subroot<r) meshs[subroot].son[1]=build_tree(subroot+1,r,(d+1)%3),expand_bynode(subroot,meshs[subroot].son[1]);
        return subroot;
    }
    void _query(const Ray& ray,int rt,double& _query_ans,Vec &my_n_ans,double&my_ord_x,double&my_ord_y,int&my_group){
        if (!meshs[rt].intersect_box(ray)) return;
        double t=meshs[rt].intersect(ray);
        if (t!=0&&t>eps&&t<_query_ans) _query_ans=t,meshs[rt].get_ord(ray.o+ray.d*t,my_n_ans,my_ord_x,my_ord_y,my_group);
        if (meshs[rt].son[0]!=-1) _query(ray,meshs[rt].son[0],_query_ans,my_n_ans,my_ord_x,my_ord_y,my_group);
        if (meshs[rt].son[1]!=-1) _query(ray,meshs[rt].son[1],_query_ans,my_n_ans,my_ord_x,my_ord_y,my_group);
    }
    double intersect(const Ray& ray,Vec& my_n_ans,double&my_ord_x,double&my_ord_y,int&my_group){
        double _query_ans=INF;
        _query(ray,root,_query_ans,my_n_ans,my_ord_x,my_ord_y,my_group);
        if (_query_ans<INF/2) return _query_ans;
        else return 0;
    }

};

int KDTree::D = 0;
Vec KDTree::n_ans=Vec();