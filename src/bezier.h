#pragma once
#include "vec.h"
#include "ray.h"
#include "utils.h"
#include "objects.h"
#include <fstream>
#include <iostream>
#include <algorithm>
class Bezier{
public:
    static int D;
	static double ord_x,ord_y;
	Texture texture;
	double refl_rate;
	static Vec n_ans;
    struct Triangle{//also Bezier Node
        Vec n;//法向
        Vec v[3];//顶点
		Vec vt[3];
		Vec vn[3];
        Vec m[2];//min and max;
        int son[2];//lson and rson
		double beta,gamma;
        double intersect(const Ray& ray){
            Vec E1=v[0]-v[1];
		    Vec E2=v[0]-v[2];
		    Vec S=v[0]-ray.o;
		    double under=(ray.d%E1).dot(E2);
		    double t=(S%E1).dot(E2)/under;
		    beta=(ray.d%S).dot(E2)/under;
		    gamma=(ray.d%E1).dot(S)/under;
		    if (t<=0||beta<0||gamma<0||beta+gamma>1) return 0;
		    return t;
        }
		void get_ord(Vec point,Vec& my_n_ans,double& my_ord_x,double& my_ord_y){
			Vec E1=v[1]-v[0];
		    Vec E0=v[2]-v[0];
		    Vec E2=point-v[0];
			double lambda1=(E1.dot(E1)*E0.dot(E2)-E0.dot(E1)*E1.dot(E2))/(E0.dot(E0)*E1.dot(E1)-E0.dot(E1)*E0.dot(E1));
			double lambda2=(E0.dot(E0)*E1.dot(E2)-E0.dot(E1)*E0.dot(E2))/(E0.dot(E0)*E1.dot(E1)-E0.dot(E1)*E0.dot(E1));
			double lambda3=1-lambda1-lambda2;
			Vec tmp=vt[0]*lambda3+vt[2]*lambda1+vt[1]*lambda2;
			my_ord_x=tmp.x,my_ord_y=tmp.y;
			my_n_ans=(vn[0]*lambda3+vn[2]*lambda1+vn[1]*lambda2).norm();
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
    Vec* v_buffer;
    Vec* vn_buffer;
	Vec* vt_buffer;
    int v_num,mesh_num;
    double query_ans;
    int root;
    
	Refl_t refl;
	double scaling_size;
	Vec P_bezier[100][100];//bezier point
	Vec Pd_bezier[100][100];
	Vec position;
	int n_bezier;

    Bezier(const char* file_in,Refl_t _refl=SPEC,double _size=300,Vec _position=Vec(10,0,10),std::string filename="")
		:refl(_refl),scaling_size(_size),position(_position),texture(filename){
		refl_rate=0;
		P_bezier[0][3]=Vec(0,80.1);
		P_bezier[0][2]=Vec(75,80);
		P_bezier[0][1]=Vec(37,27);
		P_bezier[0][0]=Vec(0,0);

		n_bezier=3;

		for (int i=0;i<=n_bezier;i++){
			P_bezier[0][i]=P_bezier[0][i]*scaling_size;
		}
		for (int i=0;i<n_bezier;i++){
			Pd_bezier[0][i]=(P_bezier[0][i+1]-P_bezier[0][i])*n_bezier;
		}
        query_ans=INF;
        if (meshs!=NULL) delete[] meshs;
        load(file_in);
        transfer(position,scaling_size);
        root=build_tree(1,mesh_num,0);
    }
    void load(const char* file_in){

        std::ifstream in;
        in.open(file_in);
        if (!in.is_open()){ fprintf(stderr, "Wrong when opening file");return;}
        in>>v_num>>v_num>>mesh_num;

        v_buffer=new Vec[v_num+10];
        vt_buffer=new Vec[v_num+10];
        meshs=new Triangle[mesh_num+10];

		for (int i=0;i<v_num;i++){
			in>>v_buffer[i].x>>v_buffer[i].y>>v_buffer[i].z;
		}
		fprintf(stderr, "\rLoading %s  1/3 finished\n ",file_in);
		for (int i=0;i<v_num;i++){
			in>>vt_buffer[i].x>>vt_buffer[i].y;
		}
		fprintf(stderr, "\rLoading %s  2/3 finished\n ",file_in);
		//int _cnt=0;
		for (int i=1;i<=mesh_num;i++){
			int tmp0,tmp1,tmp2;
			in>>tmp0>>tmp1>>tmp2;
			tmp0-=1,tmp1-=1,tmp2-=1;

            meshs[i].v[0]=v_buffer[tmp0];
            meshs[i].v[1]=v_buffer[tmp1];
            meshs[i].v[2]=v_buffer[tmp2];
            meshs[i].n=((meshs[i].v[0]-meshs[i].v[1])%(meshs[i].v[0]-meshs[i].v[2])).norm();

            meshs[i].vt[0]=vt_buffer[tmp0];
            meshs[i].vt[1]=vt_buffer[tmp1];
            meshs[i].vt[2]=vt_buffer[tmp2];

			for (int _i=0;_i<3;_i++){
				double theta=meshs[i].vt[_i].x*M_PI*2;
				double t=meshs[i].vt[_i].y;
				for (int j=0;j<n_bezier;j++){
					for (int k=0;k<n_bezier-j;k++){
						P_bezier[j+1][k]=P_bezier[j][k]*t+P_bezier[j][k+1]*(1-t);
					}
				}
				for (int j=0;j<n_bezier-1;j++){
					for (int k=0;k<n_bezier-1-j;k++){
						Pd_bezier[j+1][k]=Pd_bezier[j][k]*t+Pd_bezier[j][k+1]*(1-t);
					}
				}
				Vec dS_dtheta=Vec(-sin(theta)*P_bezier[n_bezier][0].x,0,cos(theta)*P_bezier[n_bezier][0].x);
				Vec dS_dt=Vec(cos(theta)*Pd_bezier[n_bezier-1][0].x,Pd_bezier[n_bezier-1][0].y,sin(theta)*Pd_bezier[n_bezier-1][0].x);
				meshs[i].vn[_i]=(dS_dtheta%dS_dt).norm();
			}

		}


        fprintf(stderr, "\rLoading %s  3/3 finished\n ",file_in);

        delete[] v_buffer;
        delete[] vt_buffer;     
    }

    void transfer(Vec position,double size){
        Vec min_position=Vec(INF,INF,INF);

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
    void _query(const Ray& ray,int rt,double& _query_ans,Vec &my_n_ans,double& my_ord_x,double& my_ord_y){
        if (!meshs[rt].intersect_box(ray)) return;
        double t=meshs[rt].intersect(ray);

        if (t!=0&&t>1e-1&&t<_query_ans) _query_ans=t,meshs[rt].get_ord(ray.o+ray.d*t,my_n_ans,my_ord_x,my_ord_y);
        if (meshs[rt].son[0]!=-1) _query(ray,meshs[rt].son[0],_query_ans,my_n_ans,my_ord_x,my_ord_y);
        if (meshs[rt].son[1]!=-1) _query(ray,meshs[rt].son[1],_query_ans,my_n_ans,my_ord_x,my_ord_y);
    }
    double intersect(const Ray& ray,Vec &my_n_ans,double& my_ord_x,double& my_ord_y){
        double _query_ans=INF;
        _query(ray,root,_query_ans,my_n_ans,my_ord_x,my_ord_y);
        if (_query_ans<INF/2) return _query_ans;
        else return 0;
    }
    double b_intersect(const Ray& ray,Vec& my_n_ans,double& my_ord_x,double& my_ord_y){
		double l=intersect(ray,my_n_ans,my_ord_x,my_ord_y);

		if (l==0) return 0;

		double theta=my_ord_x*M_PI*2;
		double t=my_ord_y;
		double _eps=1e-5;
		Vec f;
		Vec _P_bezier[100][100];
		Vec _Pd_bezier[100][100];
		for (int ii=0;ii<=n_bezier;ii++){
			_P_bezier[0][ii]=P_bezier[0][ii];
		}
		for (int ii=0;ii<n_bezier;ii++){
			_Pd_bezier[0][ii]=Pd_bezier[0][ii];
		}


		int cnt=0;
		do 
		{
			cnt++;
			for (int j=0;j<n_bezier;j++){
				for (int k=0;k<n_bezier-j;k++){
					_P_bezier[j+1][k]=_P_bezier[j][k]*t+_P_bezier[j][k+1]*(1-t);
				}
			}
			for (int j=0;j<n_bezier-1;j++){
				for (int k=0;k<n_bezier-1-j;k++){
					_Pd_bezier[j+1][k]=_Pd_bezier[j][k]*t+_Pd_bezier[j][k+1]*(1-t);
				}
			}
		
			f=ray.o+ray.d*l-(Vec(cos(theta)*_P_bezier[n_bezier][0].x,_P_bezier[n_bezier][0].y,
				sin(theta)*_P_bezier[n_bezier][0].x)+position);
			if (sqrt(f.dot(f))<_eps) break;
			Vec dS_dtheta=Vec(-sin(theta)*_P_bezier[n_bezier][0].x,0,cos(theta)*_P_bezier[n_bezier][0].x);
			Vec dS_dt=-Vec(cos(theta)*_Pd_bezier[n_bezier-1][0].x,_Pd_bezier[n_bezier-1][0].y,sin(theta)*_Pd_bezier[n_bezier-1][0].x);
			double D=ray.d.dot(dS_dtheta%dS_dt);
			double dl=((dS_dtheta.dot(dS_dt%f))/D);
			double dtheta=(ray.d.dot(dS_dt%f)/D);
			double dt=(ray.d.dot(dS_dtheta%f)/D);
			int cntt=0;
			while (cntt<10){
				double tmpt=t+dt;
				for (int j=0;j<n_bezier;j++){
					for (int k=0;k<n_bezier-j;k++){
						_P_bezier[j+1][k]=_P_bezier[j][k]*tmpt+_P_bezier[j][k+1]*(1-tmpt);
					}
				}
				for (int j=0;j<n_bezier-1;j++){
					for (int k=0;k<n_bezier-1-j;k++){
						_Pd_bezier[j+1][k]=_Pd_bezier[j][k]*tmpt+_Pd_bezier[j][k+1]*(1-tmpt);
					}
				}
				Vec tmpf=ray.o+ray.d*(l-dl)-(Vec(cos(theta-dtheta)*_P_bezier[n_bezier][0].x,_P_bezier[n_bezier][0].y,
					sin(theta-dtheta)*_P_bezier[n_bezier][0].x))-position;
				if (tmpf.dot(tmpf)<f.dot(f)){
					break;
				}
				dl*=0.5;
				dtheta*=0.5;
				dt*=0.5;
				cntt++;
			}
			if(cntt==10) break;
			l-=dl;
			theta-=dtheta;
			t+=dt;
			theta=std::fmod(theta,2*M_PI);
			if (cnt>15){
				return 0;
			}
		} while (1);

		for (int j=0;j<n_bezier;j++){
			for (int k=0;k<n_bezier-j;k++){
				_P_bezier[j+1][k]=_P_bezier[j][k]*t+_P_bezier[j][k+1]*(1-t);
			}
		}
		for (int j=0;j<n_bezier-1;j++){
			for (int k=0;k<n_bezier-1-j;k++){
				_Pd_bezier[j+1][k]=_Pd_bezier[j][k]*t+_Pd_bezier[j][k+1]*(1-t);
			}
		}

		Vec dS_dtheta=Vec(-sin(theta)*_P_bezier[n_bezier][0].x,0,cos(theta)*_P_bezier[n_bezier][0].x);
		Vec dS_dt=Vec(cos(theta)*_Pd_bezier[n_bezier-1][0].x,_Pd_bezier[n_bezier-1][0].y,sin(theta)*_Pd_bezier[n_bezier-1][0].x);
		my_n_ans=(dS_dtheta%dS_dt).norm();
		my_ord_x=theta/(2*M_PI);my_ord_y=t;
		return l;		


    }
};

int Bezier::D = 0;
double Bezier::ord_y = 0;
double Bezier::ord_x = 0;
Vec Bezier::n_ans=Vec();