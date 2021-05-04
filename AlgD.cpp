//	AlgD.cpp
#include <functional>
#include <numeric> 
#include <fstream>
#include <cmath>
#include <ctime>
#include <vector>
#include <string>
#include <random>
#include "AlgD.h"
#include "emo/Sel.h"
#include "../SelT.h"
#include "alg/AR.h"
#include "alg/Matrix.h"
#include "emo/GenMod.h"
#include "alg/normal.h"
#include "assert.h"
#include "../gd_generate.h"
#include <map>
#include "../WeighVector.h"
#include <algorithm>
#include <yvals.h>
#if (_MSC_VER >= 1600)
#define __STDC_UTF_16__
#endif
//#include"engine.h"
#include"memory.h"
//#pragma comment(lib,"libmat.lib") 
//#pragma comment(lib,"libmx.lib") 
//#pragma comment(lib,"libeng.lib") 
//#define SAVE_CEN 1

//!\brief	az namespace, the top namespace
namespace az
{
	//!\brief	mea namespace, the multiobjective evolutionary algorithm namespace
	namespace mea
	{
		//!\brief namespace of dynamic evolutionary algoirhtm
		namespace dea
		{
			//const double PI = 3.141592653589793;
			double		 T	= 0.0, T0;

			// FDAs are from M. Farina, K. Deb and P. Amato's paper
			//FDA1
			void FDA1(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
			{
				double gx = 1.0, G = sin(0.5*PI*T);
				for(unsigned int i=1; i<X.size(); i++)
					gx += (X[i]-G)*(X[i]-G);
				F[0] = X[0];
				F[1] = gx*(1-sqrt(F[0]/gx));
			}
			//FDA2
			void FDA2(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
			{
				double gx = 1.0, hx = 0.0, H = 0.75+0.7*sin(0.5*PI*T);
				unsigned int i,xii = (unsigned int)(X.size()/2);
				for(i=1; i<xii; i++) gx += X[i]*X[i];
				for(i=xii; i<X.size(); i++) hx += (X[i]-H)*(X[i]-H); hx += H;
				F[0] = X[0];
				F[1] = gx*(1-pow(F[0]/gx, 1.0/hx));
			}
			//FDA3
			void FDA3(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
			{
				double gx = 1.0, Ft = pow(10.0,2.0*sin(0.5*PI*T)), Gt = fabs(sin(0.5*PI*T));
				unsigned int i;
				for(i=1; i<X.size(); i++) gx += (X[i]-Gt)*(X[i]-Gt); gx += Gt;
				F[0] = pow(X[0],Ft);
				F[1] = gx*(1-pow(F[0]/gx, 0.5));
			}
			//FDA4
			void FDA4(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
			{
				double gx = 0.0, G = fabs(sin(0.5*PI*T));
				unsigned int i;
				for(i=2; i<X.size(); i++) gx += (X[i]-G)*(X[i]-G);
				F[0] = (1.0+gx)*cos(0.5*PI*X[0])*cos(0.5*PI*X[1]);
				F[1] = (1.0+gx)*cos(0.5*PI*X[0])*sin(0.5*PI*X[1]);
				F[2] = (1.0+gx)*sin(0.5*PI*X[0]);
			}
			// == DMOPs are from C.K Goh and K.C Tan 's paper
			//dMOP1
			void DMOP1(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
			{
				double gx = 0.0, H = 0.75*sin(0.5*PI*T)+1.25;
				unsigned int i;
				for(i=1; i<X.size(); i++) gx += X[i]*X[i]; gx = 1.0 + 9.0*gx;
				F[0] = X[0];
				F[1] = gx*(1-pow(F[0]/gx, H));
			}
			//dMOP2
			void DMOP2(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
			{
				double gx = 1.0, H = 0.75*sin(0.5*PI*T)+1.25, G = sin(0.5*PI*T);
				unsigned int i;
				for(i=1; i<X.size(); i++) gx += (X[i]-G)*(X[i]-G);
				F[0] = X[0];
				F[1] = gx*(1-pow(F[0]/gx, H));
			}
			//dMOP3
			void DMOP3(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
			{
				double gx = 1.0, H = 0.75*sin(0.5*PI*T)+1.25, G = sin(0.5*PI*T);
				unsigned int i;
				for(i=1; i<X.size(); i++) gx += (X[i]-G)*(X[i]-G);
				F[0] = X[0];
				F[1] = gx*(1-pow(F[0]/gx, 0.5));
			}

			// newly designed DMOPs problems
			//centre moving functions
			void LT(double t, unsigned int index, double& x1, double& x2)   
			{
				t *= 0.5;
				switch(index)
				{	
				case 0:	// Lissajous curve
					x1 = (cos(t*PI)+1.0)*2.0;
					x2 = (sin(2.0*t*PI)+1.0)*2.0;
					break;
				case 1:	// Rose curve
					x1 = (cos(1.5*t*PI)*sin(0.5*t*PI)+1.0)*2.0;
					x2 = (cos(1.5*t*PI)*cos(0.5*t*PI)+1.0)*2.0;
					break;
				case 2: // Heart curve
					x1 = ((1.0-sin(t*PI))*sin(t*PI)+2.0)*1.7;
					x2 = ((1.0-sin(t*PI))*cos(t*PI)+1.5)*1.4;
					break;
				case 3: // discontinus Lissajous curve
					t  = t-floor(t);
					x1 = (cos(t*PI)+1.0)*2.0;
					x2 = (sin(2.0*t*PI)+1.0)*2.0;
					break;
				default:
					x1 = x2 = 0.0;
					break;
				}
			}
			//shape moving functions
			double HT(double t)
			{
				return 1.25+0.75*sin(t*PI);
			}
			// problem framework
			void DMOP(std::vector< double >& F, std::vector< double >& X, unsigned int index, bool turn)
			{
				double a, b, Gi, gx1=0.0, gx2=0.0, ht=HT(T);

				bool old = (((unsigned int)(T*10.0+0.001)) % 2 == 1);

				LT(T, index, a, b);

				for(unsigned int i=1; i<X.size(); i++)
				{
					//if(turn && old)
					//	Gi = pow(fabs(2.0*X[0]-2.0*a-1.0), ht+double(i-1.0)/double(X.size()-2.0));
					//else	
					//	Gi = 1.0 - pow(fabs(2.0*X[0]-2.0*a-1.0), ht+double(i-1.0)/double(X.size()-2.0));

					if(turn && old)
						Gi = pow(fabs(X[0]-a), ht+double(i+1.0)/double(X.size()));
					else	
						Gi = 1.0 - pow(fabs(X[0]-a), ht+double(i+1.0)/double(X.size()));


					if(i % 2 == 1)
						gx1 += pow(X[i] - b - Gi ,2.0);
					else 
						gx2 += pow(X[i] - b - Gi ,2.0);
				}
				F[0] = pow(fabs(X[0]-a), ht)     + 0.5*gx1;
				F[1] = pow(fabs(X[0]-a-1.0), ht) + 0.5*gx2;

			}
			void DMOPA(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X) //F5
			{
				DMOP(F,X,0,false);
			}
			void DMOPB(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)//F6
			{
				DMOP(F,X,1,false);
			}
			void DMOPC(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)//F7
			{
				DMOP(F,X,2,false);
			}
			void DMOPD(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X) //F9
			{
				DMOP(F,X,3,false);
			}
			void DMOPE(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)//F10
			{
				DMOP(F,X,0,true);
			}
			void DMOPF(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)//F8
			{	
				double gx = 0.0, G = sin(0.5*PI*T);
				for(unsigned int i=2; i<X.size(); i++) gx += pow(X[i]-pow(0.5*(X[0]+X[1]), HT(T)+double(i+1.0)/double(2.0*X.size()))-G,2.0);
				F[0] = (1.0+gx)*cos(0.5*PI*X[0])*cos(0.5*PI*X[1]);
				F[1] = (1.0+gx)*cos(0.5*PI*X[0])*sin(0.5*PI*X[1]);
				F[2] = (1.0+gx)*sin(0.5*PI*X[0]);
			}
			//多项式偏转：P(t)
			double Pt()
			{
				return 1.0 + 16.0 * pow(cos(0.5 * PI * T) , 4);
			}

			//欺骗函数中的A(t)与多模函数中的C(t)
			double tempT()
			{
				return 0.5+0.35*sin(PI*T);
			}

			//欺骗函数  A(t),B=0.001,C=0.1
			double S_decept(double x)
			{
				double A = tempT(),B = 0.001,C = 0.1;
				assert(!(A-B));assert(!(1-A-B));
				double temp1 = 0,temp2 = 0;
				temp1 = floor(x-A+B)*(1-C+(A-B)/B)/(A-B);
				temp2 =  floor(A+B-x)*(1-C+(1-A-B)/B)/(1-A-B);
				return 1+(fabs(x-A)-B)*(temp1+temp2+1/B);
			}

			//多模函数  A=20,B=10,C(t)
			double S_multi(double x)
			{
				double A = 20,B = 10,C = tempT();
				double temp = fabs(x-C)/(2*(floor(C-x)+C));
				return (1+cos((4*A+2)*PI*(0.5-temp))+4*B*pow(temp,2))/(B+2);
			}

			//SDMOP1
			void SDMOP1(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
			{
				double gx = 1.0,Gt = cos(0.5*PI*T);
				for(unsigned int i = 1;i<X.size();i++)
					gx += pow(X[i]-Gt,2);
				F[0] = X[0];
				F[1] = gx*(1-pow(F[0]/gx,0.5));
			}

			//SDMOP2
			void SDMOP2(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
			{
				double gx = 1.0,Ht = 1.2+0.8*cos(0.5*PI*T);
				for(unsigned int i=1;i<X.size();i++)
					gx += 9*X[i]*X[i];
				F[0] = pow(X[0],Pt());
				//F[0]=X[0];
				F[1] = gx*(1-pow(F[0]/gx,Ht));
			}

			//SDMOP3
			void SDMOP3(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
			{
				double gx = 1.0,Gt = cos(0.5*PI*T),Ht = 1.2+0.8*cos(0.5*PI*T);
				for(unsigned int i=1;i<X.size();i++)
					gx += (X[i]-Gt)*(X[i]-Gt);
				F[0] = S_decept(X[0]);
				F[1] = gx*(1-pow(F[0]/gx,Ht));
			}

			//SDMOP4
			void SDMOP4(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
			{
				double Gt = cos(0.5*PI*T),gx = 1.0+fabs(Gt),sum = 0,mulit = 0;
				unsigned int i;
				for(i=0;i<5;i++)
				{sum += i;mulit += X[i]*i;}
				for(i=5;i<X.size();i++)
					gx += (X[i]-Gt)*(X[i]-Gt);
				F[0] = mulit/sum;
				//F[0] = X[0];
				F[1] = gx*(1-F[0]-cos(10*PI*F[0]+0.5*PI)/(10*PI));
			}

			//SDMOP5
			void SDMOP5(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
			{
				assert(1);
				unsigned int M = F.size(),N = X.size(),i,j;
				double gx = 0.0,temp,Gt = cos(0.5*PI*T);
				for(i=M-1;i<N;i++)
					//gx += (X[i]-Gt)*(X[i]-Gt);
						gx += S_multi(X[i])*S_multi(X[i]);
				for(i=0;i<M;i++)
				{
					temp = 1.0;
					if(i==0)
						for(j=0;j<M-1;j++)
							temp *= cos(0.5*PI*X[j]);
					else
						if(i==M-1)
							temp *= sin(0.5*PI*X[0]);
						else
						{
							temp = sin(0.5*PI*X[M-i-1]);
							for(j=0;j<M-i-1;j++)    temp *= cos(0.5*PI*X[j]); 
						}
						F[i] = (1+gx)*temp;
				}
				/*F[0] = (1.0+gx)*cos(0.5*PI*X[0])*cos(0.5*PI*X[1]);
				F[1] = (1.0+gx)*cos(0.5*PI*X[0])*sin(0.5*PI*X[1]);
				F[2] = (1.0+gx)*sin(0.5*PI*X[0]);*/
			}

			//SDMOP6
			void SDMOP6(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
			{
				assert(1);
				unsigned int M = F.size(),N = X.size(),i,j;
				double temp,Gt = cos(0.5*PI*T),gx = 0.0;
				std::vector< double > Z(M);
				for(i=0;i<M-1;i++)
					Z[i] = pow(X[i],Pt());
				for(i=M-1;i<N;i++)
					gx += (X[i]-Gt)*(X[i]-Gt);
				for(i=0;i<M;i++)
				{
					temp = 1.0;
					if(i==0)
						for(j=0;j<M-1;j++)
							temp *= cos(0.5*PI*Z[j]);
					else
						if(i==M-1)
							temp *= sin(0.5*PI*Z[0]);
						else
						{
							temp = sin(0.5*PI*Z[M-i-1]);
							for(j=0;j<M-i-1;j++)    temp *= cos(0.5*PI*Z[j]); 
						}
						F[i] = (1+gx)*temp;
				}
			}
			//SDMOP7
			void SDMOP7(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
			{
				assert(1);
				unsigned int M = F.size(),N = X.size(),i,j;
				double temp,Gt = cos(0.5*PI*T),gx = fabs(Gt);
				std::vector< double > Z(M);
				for(i=0;i<M-1;i++)
					Z[i] = S_decept(X[i]);
				for(i=M-1;i<N;i++)
					gx += (X[i]-Gt)*(X[i]-Gt);
				for(i=0;i<M;i++)
				{
					temp = 1.0;
					if(i==0)
						for(j=0;j<M-1;j++)
							temp *= cos(0.5*PI*Z[j]);
					else
						if(i==M-1)
							temp *= sin(0.5*PI*Z[0]);
						else
						{
							temp = sin(0.5*PI*Z[M-i-1]);
							for(j=0;j<M-i-1;j++)    temp *= cos(0.5*PI*Z[j]); 
						}
						F[i] = (1+gx)*temp;
				}
			}

			//SDMOP8
			void SDMOP8(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
			{
				assert(1);
				unsigned int M = F.size(),N = X.size(),i,j;
				double temp = 0,Ht = fabs(cos(0.5*PI*T)),gx = 0;
				for(i=2;i<N;i++)
					gx += 9*X[i]*X[i];
				F[0] = X[0];F[1] = X[1];

				for(i=2;i<M;i++)
				{
					temp = (X[0]/(1+gx))*(1+Ht*cos(PI*X[0]*X[0]*X[0])*cos(PI*X[0]*X[0]*X[0])) + (X[1]/(1+gx))*(1+Ht*cos(PI*X[1]*X[1]*X[1])*cos(PI*X[1]*X[1]*X[1]));
					F[i] = (1+gx)*(i-1)*(2-0.5*temp);
				}
			}

			//"JY1","JY2","JY3","JY4","JY5","JY6","JY7","JY8","JY9","JY10"
            void JY1(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
                {
                	unsigned int j, count1, count2;
                	double a, w, g, sum1;
                
                	a = 0.05, w = 6;
                	g = sin(0.5*PI*T);
                	sum1 = 0;
                	for (j = 1; j < X.size(); j++)
                	{
                		double x = X[j];
                		double yj = x - g;
                		yj = yj*yj;
                		sum1 += yj;
                	}
                	F[0] = (1 + sum1)*(X[0] + a*sin(w*PI*X[0]));
                	F[1] = (1 + sum1)*(1 - X[0] + a*sin(w*PI*X[0]));
                }
            void JY2(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
                {
                	unsigned int j, count1, count2;
                	double a, w, g, sum1;
                
                	a = 0.05, w = floor(6 * sin(0.5*PI*(T - 1)));
                	g = sin(0.5*PI*T);
                	sum1 = 0;
                
                	for (j = 1; j < X.size(); j++)
                	{
                		double x = X[j];
                		double yj = x - g;
                		yj = yj*yj;
                		sum1 += yj;
                	}
                	F[0] = (1 + sum1)*(X[0] + a*sin(w*PI*X[0]));
                	F[1] = (1 + sum1)*(1 - X[0] + a*sin(w*PI*X[0]));
                }
            void JY3(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
                {
                	unsigned int j, count1, count2;
                	double a, aa, y1, w, g, sum1;
                
                	aa = floor(100 * sin(0.5*PI*T)*sin(0.5*PI*T));
                	y1 = fabs(X[0] * sin((2 * aa + 0.5)*PI*X[0]));
                
                	a = 0.05, w = floor(6 * sin(0.5*PI*(T - 1)));
                	
                
                	sum1 = 0;
                
                	for (j = 1; j < X.size(); j++)
                	{
                		double x = 2*X[j]-1;
                		double x0 = 2*X[j - 1]-1;
                		double yj;
                		if (j == 1){ yj = x*x - y1; }
                		else{ yj = x*x - x0; }
                
                		yj = yj*yj;
                		sum1 += yj;
                	}
                	F[0] = (1 + sum1)*(y1 + a*sin(w*PI*y1));
                	F[1] = (1 + sum1)*(1 - y1 + a*sin(w*PI*y1));
                }
            void JY4(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
                {
                	unsigned int j, count1, count2;
                	double a, w, g, sum1;
                
                	g = sin(0.5*PI*T);
                	a = 0.05, w = pow(10.0, 1.0 + fabs(g));
                
                	sum1 = 0;
                
                	for (j = 1; j < X.size(); j++)
                	{
                		double x = X[j];
                		double yj = x - g;
                		yj = yj*yj;
                		sum1 += yj;
                	}
                	F[0] = (1 + sum1)*(X[0] + a*sin(w*PI*X[0]));
                	F[1] = (1 + sum1)*(1 - X[0] + a*sin(w*PI*X[0]));
					if (F[0]<0)  F[0]=0;
					if (F[1]<0)  F[1]=0;
                }
            void JY5(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
                {
                	unsigned int j, count1, count2;
                	double a, w, g, sum1;
                
                	g = sin(0.5*PI*T);
                	a = 0.3*sin(0.5*PI*(T - 1)), w = 1.0;
                
                	sum1 = 0;
                
                	for (j = 1; j < X.size(); j++)
                	{
                		double x =  X[j];
                		double yj = x;
                		yj = yj*yj;
                		sum1 += yj;
                	}
                	F[0] = (1 + sum1)*(X[0] + a*sin(w*PI*X[0]));
                	F[1] = (1 + sum1)*(1 - X[0] + a*sin(w*PI*X[0]));
                }
            void JY6(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
                {
                	unsigned int j, count1, count2;
                	double a, w, g, k, sum1;
                
                	g = sin(0.5*PI*T);
                	a = 0.1, w = 3;
                	k = 2 * floor(10 * fabs(g));
                
                	sum1 = 0;
                
                	for (j = 1; j < X.size(); j++)
                	{
                		double x =  X[j];
                		double yj = x - g;
                		yj = 4 * yj*yj - cos(k*PI*yj) + 1;
                		sum1 += yj;
                	}
                	F[0] = (1 + sum1)*(X[0] + a*sin(w*PI*X[0]));
                	F[1] = (1 + sum1)*(1 - X[0] + a*sin(w*PI*X[0]));
                
                }
            void JY7(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
                {
                	unsigned int j, count1, count2;
                	double a, w, g, s, t, sum1;
                
                	g = sin(0.5*PI*T);
                	a = 0.1, w = 3;
                	s = t = 0.2 + 2.8*fabs(g);
                
                	sum1 = 0;
                
                	for (j = 1; j < X.size(); j++)
                	{
                		double x =  X[j];
                		double yj = x - g;
                		yj = yj*yj - 10 * cos(2 * PI*yj) + 10;
                		sum1 += yj;
                	}
                	F[0] = (1 + sum1)*pow(X[0] + a*sin(w*PI*X[0]), s);
                	F[1] = (1 + sum1)*pow(1 - X[0] + a*sin(w*PI*X[0]), t);
                }
            void JY8(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
                 {
                 	unsigned int j, count1, count2;
                 	double a, w, g, s, t, sum1;
                 
                 	g = sin(0.5*PI*T);
                 	a = 0.05, w = 6;
                 	t = 10.0 - 9.8*fabs(g);
                 	s = 2.0 / t;
                 
                 	sum1 = 0;
                 
                 	for (j = 1; j < X.size(); j++)
                 	{
                 		double x = X[j];
                 		double yj = x;
                 		yj = yj*yj;
                 		sum1 += yj;
                 	}
                 	F[0] = (1 + sum1)*pow(X[0] + a*sin(w*PI*X[0]), s);
                 	F[1] = (1 + sum1)*pow(1 - X[0] + a*sin(w*PI*X[0]), t);
                 }
            void JY9(std::vector< double >& F, std::vector< double >& e, std::vector< double >& I, std::vector< double >& X)
            {
            	unsigned int j, d, count1, count2, nvar = X.size();
            		double a, w, g, sum1;
            
            		d = (int)floor(T/5) % 3;
            		g = fabs(sin(0.5*PI*T));
            		a = 0.05, w = floor(6 * pow(sin(0.5*PI*(T - 1)), d));
            
            		sum1 = 0;
            
            		for (j = 1; j < nvar; j++)
            		{
            			double x = 2*X[j]-1 ;
            			double yj = x + d - g;
            			yj = yj*yj;
            			sum1 += yj;
            		}
            		F[0] = (1 + sum1)*(X[0] + a*sin(w*PI*X[0]));
            		F[1] = (1 + sum1)*(1 - X[0] + a*sin(w*PI*X[0]));
            
            		return;
            
            }
			
			void DFRange(std::vector<double>& low, std::vector<double>& upp, std::string& name ,unsigned int mObj)
			{
				if(	name == std::string("FDA1")   || name == std::string("FDA2")   || name == std::string("FDA3")||
					name == std::string("DMOP1")  || name == std::string("DMOP2")  || name == std::string("DMOP3")||
					name == std::string("SDMOP1") || name == std::string("SDMOP2") || name == std::string("SDMOP3")
					||name == std::string("JY1")   || name == std::string("JY2")   || name == std::string("JY3")||
					name == std::string("JY4")  || name == std::string("JY5")  || name == std::string("JY6")||
					name == std::string("JY7") || name == std::string("JY8") || name == std::string("JY9") 
					)
				{
					low[0] = 0;upp[0] =  1;
					for(unsigned int i=1; i<(unsigned int)(low.size()); i++){low[i] = -1; upp[i] = 1;}
				}
				else if(name == std::string("FDA4") || name == std::string("SDMOP5") || name == std::string("SDMOP8"))
				{
					low[0] = 0;upp[0] =  1;
					for(unsigned int i=1; i<(unsigned int)(low.size()); i++){low[i] =  0; upp[i] = 1;}
				}
				else if(name == std::string("DMOPF"))
				{
					low[0] = 0; upp[0] =  1;
					low[1] = 0; upp[1] =  1;
					for(unsigned int i=2; i<(unsigned int)(low.size()); i++){low[i] = -1.0; upp[i] = 2.0;}
				}
				else if(name == std::string("SDMOP4"))
				{
					for(unsigned int i=0;i<5;i++){low[i] = 0;upp[i] = 1.0;}
					for(unsigned int i=5;i<(unsigned int)(low.size());i++){low[i] = -1.0;upp[i] = 1.0;}
				}
				else if(name == std::string("SDMOP6") || name == std::string("SDMOP7"))
				{
					for(unsigned int i=0;i<mObj-1;i++){low[i] = 0;upp[i] = 1.0;}
					for(unsigned int i=mObj-1;i<(unsigned int)(low.size());i++){low[i] = -1.0;upp[i] = 1.0;}
				}
				else
				{	//DMOPA - DMOPE
					for(unsigned int i=0; i<(unsigned int)(low.size()); i++){low[i] = 0.0; upp[i] = 5.0;}
				}
			}
			
			DMOO::DMOO(
				unsigned int	strategy,
				std::string&	optimizer,
				unsigned int	popsize	,
				unsigned int	stepmax	,
				unsigned int	taot	,
				unsigned int	nt		,
				unsigned int	torder	,
				double			t0		,
				double			alpha   ,
				CParameter&		par		,
			std::vector<int> &t_order)
				:mPop(par),mBest0(par),mBest1(par),observe(par),memory(par)
			{
				switch(strategy)
				{
				case 1:
					mStrategy = INIRIS;
					break;
				case 2:
					mStrategy = INIFPS;
					break;
				case 3:
					mStrategy = INIPPS;
					break;
				case 4:
					mStrategy = INIPZ;
					break;
				case 5:
					mStrategy = INIEGS;
					break;
				case 6:
					mStrategy = INIRG;
					break;
				case 7:
					mStrategy = INITCS;
					break;
				case 8:
					mStrategy = INITCS_2;
					break;
				case 9:
				    mStrategy = hierarchy;
				    break;
				case 10:
					mStrategy = INITEMPTY;
					break;
				default:
					mStrategy = INIRIS;
					break;
				}
				mOptimizer = optimizer;
				numMemory = 0;
				mPopSize = popsize;
				mMaxStep = stepmax;
				mTaoT	 = taot;
				mDelT	 = 1.0/nt;
				T0		 = mT0 = t0;
				mMaxOrder= torder;
				mAlpha   = alpha;
				pPar	 = &par;
				
				for(unsigned int i=0;i<t_order.size();i++)
					T_order.push_back(t_order[i]);

				if(pPar->Problem() == std::string("FDA1"))
				{
					P().Evaluator( FDA1 );
					P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
				}
				else if(pPar->Problem() == std::string("FDA2"))
				{
					P().Evaluator( FDA2 );
					P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
				}
				else if(pPar->Problem() == std::string("FDA3"))
				{
					P().Evaluator( FDA3 );
					P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
				}
				else if(pPar->Problem() == std::string("FDA4"))
				{
					P().Evaluator( FDA4 );
					P().FSize( 3 ); P().ESize( 0 ); P().ISize( 0 );
				}
				else if(pPar->Problem() == std::string("DMOP1"))
				{
					P().Evaluator( DMOP1 );
					P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
				}
				else if(pPar->Problem() == std::string("DMOP2"))
				{
					P().Evaluator( DMOP2 );
					P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
				}
				else if(pPar->Problem() == std::string("DMOP3"))
				{
					P().Evaluator( DMOP3 );
					P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
				}
				else if(pPar->Problem() == std::string("DMOPA"))
				{
					P().Evaluator( DMOPA );
					P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
				}
				else if(pPar->Problem() == std::string("DMOPB"))
				{
					P().Evaluator( DMOPB );
					P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
				}
				else if(pPar->Problem() == std::string("DMOPC"))
				{
					P().Evaluator( DMOPC );
					P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
				}
				else if(pPar->Problem() == std::string("DMOPD"))
				{
					P().Evaluator( DMOPD );
					P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
				}
				else if(pPar->Problem() == std::string("DMOPE"))
				{
					P().Evaluator( DMOPE );
					P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
				}
				else if(pPar->Problem() == std::string("DMOPF"))
				{
					P().Evaluator( DMOPF );
					P().FSize( 3 ); P().ESize( 0 ); P().ISize( 0 );
				}
				else if(pPar->Problem() == std::string("SDMOP1"))
				{
					P().Evaluator( SDMOP1 );
					P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
				}
				else if(pPar->Problem() == std::string("SDMOP2"))
				{
					P().Evaluator( SDMOP2 );
					P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
				}
				else if(pPar->Problem() == std::string("SDMOP3"))
				{
					P().Evaluator( SDMOP3 );
					P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
				}
				else if(pPar->Problem() == std::string("SDMOP4"))
				{
					P().Evaluator( SDMOP4 );
					P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
				}
				else if(pPar->Problem() == std::string("SDMOP5"))
				{
					P().Evaluator( SDMOP5 );
					P().FSize( 3 ); P().ESize( 0 ); P().ISize( 0 );
				}
				else if(pPar->Problem() == std::string("SDMOP6"))
				{
					P().Evaluator( SDMOP6 );
					P().FSize( 3 ); P().ESize( 0 ); P().ISize( 0 );
				}
				else if(pPar->Problem() == std::string("SDMOP7"))
				{
					P().Evaluator( SDMOP7 );
					P().FSize( 3 ); P().ESize( 0 ); P().ISize( 0 );
				}
				else if(pPar->Problem() == std::string("SDMOP8"))
				{
					P().Evaluator( SDMOP8 );
					P().FSize( 3 ); P().ESize( 0 ); P().ISize( 0 );
				}
				
				else if(pPar->Problem() == std::string("JY1")){
					P().Evaluator(JY1);
					P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
				}
				else if(pPar->Problem() == std::string("JY2")){
					P().Evaluator(JY2);
					P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
				}
				else if(pPar->Problem() == std::string("JY3")){
					P().Evaluator(JY3);
					P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
				}
				else if(pPar->Problem() == std::string("JY4")){
					P().Evaluator(JY4);
					P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
				}
				else if(pPar->Problem() == std::string("JY5")){
					P().Evaluator(JY5);
					P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
				}
				else if(pPar->Problem() == std::string("JY6")){
					P().Evaluator(JY6);
					P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
				}
				else if(pPar->Problem() == std::string("JY7")){
					P().Evaluator(JY7);
					P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
				}
				else if(pPar->Problem() == std::string("JY8")){
					P().Evaluator(JY8);
					P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
				}
				else if(pPar->Problem() == std::string("JY9")){
					P().Evaluator(JY9);
					P().FSize( 2 ); P().ESize( 0 ); P().ISize( 0 );
				}
				std::vector<double>point(pPar->FSize());
				az::wv::NBI_pop(pPar->XSize(),point,lambda,mPopSize);
				Reset();
			}

			void DMOO::Reset()
			{
				mStep		= 0;
				mEvas		= 0;
				T			= mT0;
				mbToChange	= false;
				hC.clear();
				DFRange(P().XLow(), P().XUpp(), P().Problem(),P().FSize());	

				g.resize(P().XSize(),0.0);
				az::rnd::seed((long) time(NULL));

			}

			
			unsigned int DMOO::Step()
			{
				if(mStep%mTaoT == 0)
				{
					switch(mStrategy)
					{
					case INIRIS:
						InitRIS(true);
						break;
					case INIFPS:
						InitFPS();
						break;
					case INIPPS:
						InitPPS();
						break;
					case INITCS:
						InitTCS();
						break;
					case INITCS_2:
					     InitTCS_2();
					     break;
					case INIPZ:
						{
							InitGIPS();
							if (mStep!=0)
							{
								mStep = mStep + 1;
							}
							break;
						}
					case INIRG:
						{
							InitRG();
							break;
						}
                    case hierarchy:
						{
							Inithierarchy();
							break;
						}

					case INIEGS:
						InitEGS();
						mStep = mStep+2;
						break;
					case INITEMPTY:
						InitEmpty();
						break;
					}
					
					Population().Evaluate(); 
					mEvas += Population().Size();
				}
				else
				{
					if(mOptimizer == std::string("MOEAD") ){
						
						az::mea::sel::MOEAD moead;
						if (CurStep() == 1){
							moead.initializeNeighborhood(neighborhood,lambda);
						}
						moead.updatePopulation(Population(),neighborhood,lambda,2);

					}else{
					
					az::mea::sel::SCrowd4 sel;
					sel.Set(0.2,1);
					sel.InitDBEA(Population(), mPopSize, get_points(), get_subPop_index(),get_subPop_index_number());
					points_fit(get_points_fitness(),get_points());	
					get_P_archive(Population(),get_archive());
					CPopulationMO popnew(P());
                    Generate(popnew, mPopSize);

					Check(popnew);
					Check(popnew, Population());

					popnew.Evaluate(); mEvas += popnew.Size();

				
					Population().Combine(popnew);					
					az::mea::sel::SCrowd2 sel1;
					sel1.Select(Population(), mPopSize);
					}
				}
				mStep++;
				if(mStep > 1 && (mStep % mTaoT == 0))
				{
					if(P().Problem().substr(0,5) == "SDMOP")
						T=mDelT*T_order[(mStep/mTaoT)%21];						
					else
						T += mDelT;													

					mbToChange = true;
				}
				else
				{
					mbToChange = false;
				}

				return mStep;
			}
			
		    void DMOO::get_P_archive(CPopulationMO pop,std::vector<int>& archive){
			archive.clear();
			unsigned int start, end;
		    start = end = 0;		       
		    while(end<pop.Size() && pop[end].Rank() == pop[start].Rank()) 
			{				
				archive.push_back(end);
				end++;
			}
		}
			void DMOO::Check(CPopulationMO& pop)
			{
				for(unsigned int i=0; i<pop.Size(); i++) pop[i].Check();
			}
			void DMOO::Check(CPopulationMO& popnew, CPopulationMO& pop)
			{
				CPopulationMO tmp(P());
				for(unsigned int i=0; i<popnew.Size(); i++) 
					if(!pop.IsContain(popnew[i])) tmp.Combine(popnew.At(i));
				popnew.Clear();
				popnew.Combine(tmp);
			}

			bool DMOO::EnvironmentChange()
			{
				if(mStep == 0) return true;

				CIndividualMO ind(P());
				unsigned int i, j, in, cnum = (unsigned int)(Population().Size()*0.05); if(cnum<1) cnum = 1;
				for(i=0; i<cnum; i++)
				{
					in  = rnd::rand((unsigned int)0, Population().Size());
					ind = Population()[in];
					ind.Evaluate(); mEvas++;
					for(j=0; j<Population().P().FSize(); j++)
						if(fabs(ind.F(j)-Population()[in].F(j))>1.0E-10) 
							return true;
				}
				return false;
			}

			void DMOO::EnvironmentChangeDegree(CPopulationMO& pop)
			{
				unsigned int i,j,k;
				double tempVar = 0.0, tempInd = 0.0;
				unsigned int N = pop.Size(),M = pop.P().FSize();
				CPopulationMO temp_pop(pop);
				pop.degreeAssign(0.0);
				std::vector<std::vector <double>> mObj;
				std::vector<double> temp;

				std::vector<double> tempVec(N);

				for(i=0;i<N;i++)
				{
					for(j=0;j<M;j++)
					{
						temp.push_back(temp_pop[i].F(j));
					}
					mObj.push_back(temp);
					temp.clear();
				}

				temp_pop.Evaluate();

				for(i=0;i<N;i++)
				{
					for(j=0;j<M;j++)
					{
						tempInd += fabs(temp_pop[i].F(j)-mObj[i][j]);
					}	
					pop[i].degreeAssign(tempInd);
					tempVar += pop[i].indDegree;
				}
				tempVar /= N;
				pop.degreeAssign(tempVar);
			}

			
			void DMOO::Predict()
			{
				unsigned int i, k, d, order, dim=(unsigned int)(mC.size());
				std::list< std::vector<double> >::iterator it, it0;

				while(hC.size()>20+mMaxOrder) hC.pop_back();

				order = std::min((unsigned int)hC.size()-1, (unsigned int)mMaxOrder);
				double **px, **pa, *pv;
				px = new double*[dim]; for(i=0; i<dim; i++) px[i] = new double[hC.size()];
				pa = new double*[dim]; for(i=0; i<dim; i++) pa[i] = new double[order];
				pv = new double[dim];
				k  = (unsigned int)hC.size()-1;
				it = hC.begin();
				d  = 0;
				while(it!=hC.end())
				{
					for(i=0; i<dim; i++) px[i][k] = (*it)[i];
					k--;
					it++;
					d++;
				}
				alg::aruv(px, dim, (unsigned int)hC.size(), order, pa, pv);
				mStdC.resize(dim); for(i=0; i<dim; i++) mStdC[i] = pv[i];
				it = hC.begin();
				for(k=0; k<order; k++)
				{
					for(i=0; i<dim; i++) pC[i] += (*it)[i]*pa[i][k];
					it++;
				}
				it = hC.begin();
				for(i=0; i<dim; i++)
				{
					if(     pC[i]>P().XUpp(i%P().XSize())) pC[i] = (*it)[i];
					else if(pC[i]<P().XLow(i%P().XSize())) pC[i] = (*it)[i];
					if(     pC[i]>P().XUpp(i%P().XSize())) pC[i] = P().XUpp(i%P().XSize());
					else if(pC[i]<P().XLow(i%P().XSize())) pC[i] = P().XLow(i%P().XSize());
				}

				for(i=0; i<dim; i++) delete []px[i]; delete []px;
				for(i=0; i<dim; i++) delete []pa[i]; delete []pa;
				delete []pv;
			}

			void DMOO::InitEGS()
			{
				//他人预测策略
			}

			void DMOO::InitPZ()
			{	
				//他人预测策略
			}

			void DMOO::InitGIPS()
			{
				//他人预测策略
			}

			void DMOO::InitPMS()
			{
				//他人预测策略
			}

			void DMOO::InitRG()
			{
				//他人预测策略
			}

			void DMOO::InitPPS()
			{	
				//他人预测策略
			}

			void DMOO::InitFPS()
			{
				//他人预测策略
			}

			void DMOO::InitRIS(bool all)
			{
				//他人预测策略
			}

			void DMOO::InitEmpty(){
				unsigned int i,j;
				if(CurStep()==0)
				{
					Population().Resize(mPopSize);
					for(i=0; i<mPopSize; i++) for(j=0; j<P().XSize(); j++) Population()[i][j] = rnd::rand(P().XLow(j),P().XUpp(j));
				}
				else
				{
					Population().Evaluate();
				}
			}

			// generate new trial solutions by any offspring generator
			CPopulationMO& DMOO::Generate(CPopulationMO& popnew, unsigned int size)
			{
				if(mOptimizer == std::string("GTM"))
				{
					az::mea::gen::mod::ModelGTM2 gen;
					gen.Set(25,2,popnew.P().FSize()-1,30,0.25);
					gen.Generate(size, popnew, Population());
				}
				else if(mOptimizer == std::string("PCA"))
				{
					az::mea::gen::mod::RM gen;
					gen.Set(popnew.P().FSize()-1, 5, 30, 0.25);
					gen.Generate(size, popnew, Population());
				}
				else if(mOptimizer == std::string("NSDE"))
				{
					popnew.P().ETA_PM() = 20.0;
					popnew.P().ETA_SBX()= 20.0;
					popnew.P().Pc()		= 0.8;
					popnew.P().Pm()		= 0.05;
					az::mea::gen::XNSDE	gen;
					gen.Set(0.5,1.0);
					gen.Generate(size, popnew, Population());
				}
				else if(mOptimizer == std::string("NSGA"))
				{
					popnew.P().ETA_PM() = 20.0;
					popnew.P().ETA_SBX()= 20.0;
					popnew.P().Pc()		= 0.8;
					popnew.P().Pm()		= 0.05;
					az::mea::gen::XSBX gen;
					gen.Generate(size, popnew, Population());
				}
				else if(mOptimizer == std::string("RM2"))
				{
					az::mea::gen::mod::RM2 gen;
					gen.Set(1.0, 1.0, 0.8, 1,30);
					gen.Generate(size, popnew, Population());
				}
				else if (mOptimizer == std::string("DEE"))
				{
					gd gen;
					gen.generater(size, popnew, Population());
					mStep = mStep + (mTaoT - 2);
				}
				else if (mOptimizer == std::string("DBEA"))
				{
					az::mea::gen::XNSDE	gen;
					gen.Set(0.1,1.0);		
					gen.Generate(size, popnew, Population(), get_points(), get_points_fitness(), get_subPop_index(),get_archive());
				}
				return popnew;
			}

			int DMOO::AddToAppSet(CPopulationMO &pop,std::vector<CIndividualMO> &set,int &grad,CIndividualMO ind)
			{
				int sizes=set.size();
				if (sizes==0)
				{
					set.push_back(ind);
					grad++;
					return grad;
				}
				else
				{
					int count,um;  count=0;
					for (int j = 0; j < sizes; j++)
					{
						if(ind.CDominate(set[j])!=(-1))
						{
							um=0;
							for (int i = 0; i < pop.P().XSize(); i++)
							{
								if(ind[i]==set[j][i])
									um++;
							}
							if(um!=pop.P().XSize() || um==0)
								count++;
						}
						if (ind.CDominate(set[j])==1)
						{
							set.erase(set.begin()+j);
							j--;
							sizes--;
						}
					}  
					if(count==sizes)
					{
						set.push_back(ind);

						grad++;
					}
					return grad;
				}
			}

			int DMOO::Grad(CPopulationMO &pop,CIndividualMO ind,int h,int grad,bool &improve ,std::vector<CIndividualMO> &set)
			{grad=AddToAppSet(pop,set,grad,pop[h]);
				if(pop[h].CDominate(ind)==1)
				{
					grad++;
					improve=true;
				}
				else
				{
					if(pop[h].CDominate(ind)==0)
					{
						if(rnd::rand()<0.95)
						{
							if(pop[h].TDominate1(ind)==1)
							{
								grad++;
								improve=true;
							}
						}
						else
						{
							unsigned int i, j, k, kk, df = pop.P().FSize();
							std::vector< std::vector< unsigned int > >	index(df);
							// sort in each dimension
							for(i=0; i<df; i++)
							{
								index[i].resize(pop.Size()); 
								for(j=0; j<pop.Size(); j++) index[i][j] = j;
								for(j=0; j<pop.Size(); j++) 
									for(k=j+1; k<pop.Size(); k++)
										if(pop[index[i][j]].F(i) > pop[index[i][k]].F(i))
											std::swap(index[i][j], index[i][k]);
							}

							float distantp,distanti;distantp=0,distanti=0;
							for (int i = 0; i < df; i++)
							{
								for (int j = 0; j < pop.Size(); j++){
									if(h==index[i][j]){
										if((j==0)||(j==pop.Size()-1))
										distantp=2.0E100;
									else 
										distantp=distantp+(pop[j-1][i]-pop[j+1][i])*(pop[j-1][i]-pop[j+1][i]);
									}
								}
							}
							CIndividualMO temp(P());
							temp=pop[h];
							pop[h]=ind;							
							for (int i = 0; i < df; i++)
							{
								for (int j = 0; j < pop.Size(); j++)
							{
								if(h==index[i][j]){if((j==0)||(j==pop.Size()-1))distantp=2.0E100;else distanti=distanti+(pop[j-1][i]-pop[j+1][i])*(pop[j-1][i]-pop[j+1][i]);}
								}
							}
							pop[h]=temp;
							if(distantp>distanti)
							{
								grad++;
								improve=true;
							}
						}						
					}				
				}
				return  grad;
			}

			int DMOO::RegionSearch1(CPopulationMO &pop,int k,bool &improve ,std::vector<CIndividualMO> &set ,std::vector<double> SR)
			{		
				int seed_init;
				seed_init = 123456789 + rnd::rand(0,1000);
				std::vector<int> orientation;
				orientation.resize(pop.P().XSize());
				for (int i = 0; i < orientation.size(); i++)
					orientation[i]=rnd::randoneortwo();
				int grad=0;
				for (int i = 0; i < pop.P().XSize(); i++)
				{
					double hg;
					hg=normal::r4_normal(0,fabs(SR[i]),&seed_init);
					
					CIndividualMO ind(pop.P());
					ind=pop[k];		
					
					pop[k][i]=pop[k][i]+hg*orientation[i];
					if(pop[k][i]>P().XUpp(i)) pop[k][i] = P().XUpp(i) - (pop[k][i]-P().XUpp(i))/10.0;
					else if(pop[k][i]<P().XLow(i)) pop[k][i] = P().XLow(i) + (P().XLow(i)-pop[k][i])/10.0;					
					pop[k].Evaluate();
					ind.Evaluate();		
					improve=false;
					grad=Grad(pop,ind,k,grad,improve,set);				
					if (improve==false)
					{
						pop[k]=ind;
						
						pop[k][i]=pop[k][i]-hg*orientation[i];
						
						if(pop[k][i]>P().XUpp(i)) 
							pop[k][i] = P().XUpp(i) - (pop[k][i]-P().XUpp(i))/10.0;
						else if(pop[k][i]<P().XLow(i)) 
							pop[k][i] = P().XLow(i) + (P().XLow(i)-pop[k][i])/10.0;							
						pop[k].Evaluate();

						//}

						grad=Grad(pop,ind,k,grad,improve,set);
						if (improve==false)
							pop[k]=ind;
					}
				}
				return grad;
			}

		


				
			void DMOO::OuInitDBEA(){
				if(mStep == 0) {InitRIS(true); return;}
				int i, j, k,index, xdim=P().XSize(),fsize=P().FSize(), psize = Population().Size();	
				assert(fsize>1);
				std::vector<std::vector<double>> points;
				std::vector<double>point(fsize);
				_normal(0,xdim,point,points);
	
				std::vector<std::vector<int>> subPop_index;
				subPop_index.resize(points.size());
	
				Population().RankSort();

				getSubPop(subPop_index,points,psize,fsize);

				std::vector<std::vector<double>> subMc;
				subMc.resize(points.size());
				for(i=0;i<subMc.size();i++)
					subMc[i].resize(xdim);
				std::vector<double> sMc;
				sMc.resize(xdim);
				int nodominateCount=0;
				std::vector<int> istrue;
				istrue.resize(points.size());

				for(i=0;i<subPop_index.size();i++){
					for(k=0;k<xdim;k++){
						sMc[k]=0.0;
						nodominateCount=0;
						for(j=0;j<subPop_index[i].size();j++){
							if(Population()[subPop_index[i][j]].Rank()==1){
								nodominateCount++;
								sMc[k] +=Population()[subPop_index[i][j]][k];
							}
						}
						sMc[k] /= double(nodominateCount);
						subMc[i][k]=sMc[k];
						if(nodominateCount!=0)
							istrue[i]=1;
						else
						{
							istrue[i] =0;
						}
					}
				}
	
			
				for(i=0;i<points.size();i++) { //		
					if(istrue[i]==0){ 
						int index1=i;  int index2=i;
						while (1)
						{
							if(istrue[index1]==0){
								index1=rnd::rand(0,i);
								index1=index1%(i);
							}
							if(istrue[index2]==0){
								index2=rnd::rand(0,i);
								index2=index2%(points.size());
							}
							if(index1==index2 ){
								index1=rnd::rand(0,i);
								index1=index1%(points.size());
							}				
							if(istrue[index1]!=0 && istrue[index2]!=0){
								break;
							}
						}
						for(k=0;k<xdim;k++){	
							subMc[i][k]=subMc[index1][k]-(index1-i)*(subMc[index1][k]-subMc[index2][k])/(double)(index1-index2);
						}
					}
				}
				for(i=0;i<subPop_index.size();i++){
					for(j=0;j<subPop_index[i].size();j++){
						if(Population()[subPop_index[i][j]].Rank()!=1){
							int direction=0;
							for(k=0;k<xdim;k++){
								if(subMc[i][k]-Population()[subPop_index[i][j]][k]>0){
									direction=1;
								}else
									direction=-1;
                                int seed_init;
							Population()[subPop_index[i][j]][k]=Population()[subPop_index[i][j]][k]+direction*normal::r4_normal(0,fabs(subMc[i][k]-Population()[subPop_index[i][j]][k]),&seed_init);				

							}
							for(k=0;k<xdim;k++)
							{
								if (Population()[subPop_index[i][j]][k]>P().XUpp(k))
								{
									Population()[subPop_index[i][j]][k] = P().XUpp(k) - (Population()[subPop_index[i][j]][k]-P().XUpp(k))/10.0;
								}
								if (Population()[subPop_index[i][j]][k]<P().XLow(k))
								{
									Population()[subPop_index[i][j]][k] = P().XLow(k) + (P().XLow(k)-Population()[subPop_index[i][j]][k])/10.0;
								}		
							}
						}
					}
				}
				Population().Evaluate();
				Population().RankSort();
				getSubPop(subPop_index,points,psize,fsize);

				
				int pops=0;
				int ranks=1;
				int sub_rank;
				while (pops<=mPopSize)
				{
					for(i=0;i<subPop_index.size();i++){
						sub_rank=mPopSize;
						for(j=0;j<subPop_index[i].size();j++){
							if(Population()[subPop_index[i][j]].Rank()==ranks){
								sub_rank=ranks;
								break;
							}else if(Population()[subPop_index[i][j]].Rank()>(ranks-1)
								&&sub_rank<Population()[subPop_index[i][j]].Rank()){
									sub_rank=Population()[subPop_index[i][j]].Rank();
							}
						}
						for(j=0;j<subPop_index[i].size();j++){
							if(Population()[subPop_index[i][j]].Rank()==sub_rank){
								pops++;
								Population()[subPop_index[i][j]].Rank(1);
							}
						}
					}
					ranks++;
				}
	
				for(i=0; i<int(Population().Size()-1); i++) 
					for(j=i+1; j<int(Population().Size()); j++)
						if(Population().In(j)<Population().In(i)) Population().Swap(i,j);

				Population().IsSort(true);
				az::mea::sel::SCrowd2 sel;
				sel.SelectSort(Population(),mPopSize).Erase(mPopSize);
				subPop_index.clear();
				points.clear();
				point.clear();
			}
			
			void DMOO::InitDBEA(){
				if(mStep == 0) {InitRIS(true); return;}
				int i, j, k,index, xdim=P().XSize(),fsize=P().FSize(), psize = Population().Size();	
				assert(fsize>1);
				std::vector<std::vector<double>> points;
				std::vector<double>point(fsize);
				_normal(0,xdim,point,points);
	
				std::vector<std::vector<int>> subPop_index;
				subPop_index.resize(points.size());
	
				Population().RankSort();

				getSubPop(subPop_index,points,psize,fsize);

				std::vector<std::vector<double>> subMc;
				subMc.resize(points.size());
				for(i=0;i<subMc.size();i++)
					subMc[i].resize(xdim);
				std::vector<double> sMc;
				sMc.resize(xdim);
				int nodominateCount=0;
				std::vector<int> istrue;
				istrue.resize(points.size());

				for(i=0;i<subPop_index.size();i++){
					for(k=0;k<xdim;k++){
						sMc[k]=0.0;
						nodominateCount=0;
						for(j=0;j<subPop_index[i].size();j++){
							if(Population()[subPop_index[i][j]].Rank()==1){
								nodominateCount++;
								sMc[k] +=Population()[subPop_index[i][j]][k];
							}
						}
						sMc[k] /= double(nodominateCount);
						subMc[i][k]=sMc[k];
						if(nodominateCount!=0)
							istrue[i]=1;
						else
						{
							istrue[i] =0;
						}
					}
				}
	

				for(i=0;i<points.size();i++) { //		
					if(istrue[i]==0){ 
						int sizes=points.size()-1;  int index1=i;  int index2=i;
						while (1)
						{
							if(istrue[index1]==0){
								index1=rnd::rand(0,sizes);
								index1=index1%(sizes);
							}
							if(istrue[index2]==0){
								index2=rnd::rand(0,sizes);
								index2=index2%(sizes);
							}
							if(index1==index2 ){
								index1=rnd::rand(0,sizes);
								index1=index1%(sizes);
							}				
							if(istrue[index1]!=0 && istrue[index2]!=0){
								break;
							}
						}
						for(k=0;k<xdim;k++){
							subMc[i][k]=subMc[index1][k]-(index1-i)*(subMc[index1][k]-subMc[index2][k])/(double)(index1-index2);
						}
					}
				}
				for(i=0;i<subPop_index.size();i++){
					for(j=0;j<subPop_index[i].size();j++){
						if(Population()[subPop_index[i][j]].Rank()!=1){
							int direction=0;
							for(k=0;k<xdim;k++){
								if(subMc[i][k]-Population()[subPop_index[i][j]][k]>0){
									direction=1;
								}else
									direction=-1;
								Population()[subPop_index[i][j]][k]=Population()[subPop_index[i][j]][k]+direction*0.001;	
                        
							}
							for(k=0;k<xdim;k++)
							{
								if (Population()[subPop_index[i][j]][k]>P().XUpp(k))
								{
									Population()[subPop_index[i][j]][k] = P().XUpp(k) - (Population()[subPop_index[i][j]][k]-P().XUpp(k))/10.0;
								}
								if (Population()[subPop_index[i][j]][k]<P().XLow(k))
								{
									Population()[subPop_index[i][j]][k] = P().XLow(k) + (P().XLow(k)-Population()[subPop_index[i][j]][k])/10.0;
								}		
							}
						}
					}
				}
				Population().Evaluate();
				Population().RankSort();
				getSubPop(subPop_index,points,psize,fsize);
				int pops=0;
				int ranks=1;
				int sub_rank;
				while (pops<=mPopSize)
				{
					for(i=0;i<subPop_index.size();i++){
						sub_rank=mPopSize;
						for(j=0;j<subPop_index[i].size();j++) {
							if(Population()[subPop_index[i][j]].Rank()==ranks){
								sub_rank=ranks;
								break;
							}else if(Population()[subPop_index[i][j]].Rank()>(ranks-1)
								&&sub_rank>Population()[subPop_index[i][j]].Rank()){
									sub_rank=Population()[subPop_index[i][j]].Rank();
							}
						}
						for(j=0;j<subPop_index[i].size();j++){
							if(Population()[subPop_index[i][j]].Rank()==sub_rank){
								pops++;
								Population()[subPop_index[i][j]].Rank(ranks);
							}
						}
					}
					ranks++;
				}
	
				for(i=0; i<int(Population().Size()-1); i++) 
					for(j=i+1; j<int(Population().Size()); j++)
						if(Population().In(j)<Population().In(i)) Population().Swap(i,j);

				Population().IsSort(true);
				az::mea::sel::SCrowd sel;
				sel.SelectSort(Population(),mPopSize).Erase(mPopSize);
				subPop_index.clear();
				points.clear();
				point.clear();
			}
			
			void DMOO::getSubPop(std::vector<std::vector<int>> &subPop_index,std::vector<std::vector<double>> points,int psize,int fsize){ 
				std::vector<double> object;
				object.resize(fsize);
				int index=0,i,j,k;
				for(i=0;i<psize;i++){
					double min_dis=1e100;
		
					for(k=0;k<fsize;k++) 
						object[k] = Population()[i].F(k);
					for(j=0;j<points.size();j++){
						if(min_dis>distance(object,points[j])){
							min_dis=distance(object,points[j]);
							index=j;
						}
					}
					subPop_index[index].push_back(i);			

				}
			}

			double DMOO::distance(std::vector<double> &referencepoint,std::vector<double> & objective){	
				double factor = std::inner_product(objective.begin(), objective.end(), referencepoint.begin(),
					 (double)0) / std::inner_product(objective.begin(), objective.end(), objective.begin(), (double)0);
				 double sum = 0;
				for (size_t i = 0; i < objective.size(); ++i)
				{				
					double temp = factor * objective[i] - referencepoint[i];
					sum += temp * temp;
				}
				assert(sum >= 0);
				return sqrt(sum);
			}

			void DMOO::_normal(const size_t component,const size_t division,std::vector<double>&point,std::vector<std::vector<double>>&points)
			{
				assert(0<=component&&component<point.size());
				if(component == point.size()-1)
				{
					point[component] = 1 - std::accumulate(point.begin(),--point.end(),(double)0);//accumulate 求数组point和
					points.push_back(point);
				}
				else
				{
					assert(division > 1);
					for (size_t i=0;i <= division; ++i)
					{
						point[component]=(double)i/division;
						if(std::accumulate(point.begin(),point.begin()+component+1,(double)0)>1)
							break;
						_normal(component+1,division,point,points);
					}
				}
			}
			 
			void DMOO::points_fit(std::vector<std::vector<double>>& points_fitness,std::vector<std::vector<double>> point)
			{
				points_fitness.clear();
				int i,j,k;
				
				double breadth=0;
				points_fitness.resize(point.size());
				
				for (i = 0; i < point.size(); i++)
				{
					points_fitness[i].resize(point.size(),0.0);
					for (j = 0; j < point.size(); j++)
					{
						if(i==j)  points_fitness[i][j]=0.0;
						else
						{
							for (k = 0; k < Population().P().FSize(); k++)
							{
							  points_fitness[i][j]=points_fitness[i][j]+std::fabs(pow(point[i][k]-point[j][k],2));
							}
							points_fitness[i][j]=points_fitness[i][j];
						}
					}
				}
				std::vector<double> row_sum,sum; 
				row_sum.resize(point.size(),0.0);
				sum.resize(point.size(),0.0);
				for (i = 0; i < point.size(); i++)
				{
					for (j = 0; j < point.size(); j++)
					{
						row_sum[i]=row_sum[i]+points_fitness[i][j];
					}
				}
				for (i = 0; i < point.size(); i++)
				{
					for (j = 0; j < point.size(); j++)
					{
						points_fitness[i][j]=points_fitness[i][j]/row_sum[i];
					}
				}
				for (i = 0; i < point.size(); i++)
				{
					for (j = 1; j < point.size(); j++)
					{
						points_fitness[i][j]=points_fitness[i][j]+points_fitness[i][j-1];
					}
				}

			}	



			void DMOO::Inithierarchy()
			{	
				
				if(mStep == 0) {
					InitRIS(true);
					oldC.resize(Population().P().XSize());
				    for(int i=0;i<Population().P().XSize();i++)oldC[i]=0.0; 
				} 
				CPopulationMO old_P(P()),old_memory(P());
				old_P=Population();
				old_memory=memory;
				oldC=mC;
				unsigned int i, j, k, xdim=P().XSize(), psize = Population().Size();
				mC.resize(xdim); 
				unsigned int num_sub1;			
				CPopulationMO population1(P());
	            CPopulationMO population2(P());
	            CPopulationMO population3(P());
                for(k=0; k<xdim; k++) 
				{
					num_sub1=0;
					mC[k]  = 0.0;
					for(i=0; i<psize; i++)
						if(old_P[i].Rank()==1)
						{ 
							mC[k] += old_P[i][k]; 
							num_sub1++;
						}
					mC[k] /= double(num_sub1);
				}
					  
				if (num_sub1<psize/20)
				{
					for (i=0;i<psize/20;i++)
						memory.Copy(old_P[i]);
				}	
				else
				{
					std::vector<size_t> vnIndex(num_sub1);
					for (size_t j = 0; j < vnIndex.size(); ++j)
						vnIndex[j] = j;
					std::random_shuffle(vnIndex.begin(), vnIndex.end());
					for (size_t j = 0; j < psize/20; ++j)
						memory.Copy(old_P[vnIndex[j]]);
				}
				if(memory.Size()>2*psize) 
				{
					memory.Pop(memory.Size()-2*psize);
				}				
				Population().Evaluate();
				Population().RankSort();
				num_sub1=0;					
				for(i=0; i<psize; i++)
					if(Population()[i].Rank()==1)
					{ 	
						population1.Combine(Population()[i]);  
						num_sub1++;
					}				
			   unsigned int num_sub2,num_sub3;
			   num_sub2=(psize-num_sub1)/2;
			   num_sub3=psize-num_sub1-num_sub2;
			   for (i = 0; i < population1.Size(); i++)
			   {
				   for (j = 0; j < P().XSize(); j++)
				   {
					   population1[i][j]=population1[i][j]+(mC[j]-oldC[j]);
					   if(population1[i][j]>P().XUpp(j)) 
							population1[i][j]  = P().XUpp(j) - (population1[i][j]-P().XUpp(j))/10.0;
						else
							if(population1[i][j]<P().XLow(j)) 
								population1[i][j]  = P().XLow(j) + (P().XLow(j)-population1[i][j])/10.0;
				   }				   
			   }
               for (i = 0; i < num_sub2; i++)
			   {
				   population2.Combine(Population()[num_sub1+i]);
			   }		   
			   
			   CIndividualMO mC2(P()),Old_mC2(P());
			  for(k=0; k<xdim; k++) 
				{
					mC2[k]  = 0.0;
					for(i=0; i<population2.Size(); i++)		mC2[k] += population2[i][k]; 							
					mC2[k] /= double(population2.Size());
				}
			  mC2.Evaluate();
			  Old_mC2=mC2;
			  bool improve;
			  std::vector<double> SR; SR.resize(P().XSize());
			  for (i = 0; i < P().XSize(); i++)
			  {
				  SR[i]=mC[i]-oldC[i];
			  }
			
			  std::vector<CIndividualMO> set;
			  for (i = 0; i < population2.Size(); i++)
			  {
				  RegionSearch1(population2,i,improve,set,SR);
			  }
			 for (int ty = 0; ty < set.size(); ty++)		population2.Combine(set[ty]);

			   for (i = 0; i < num_sub3; i++)
			   {
				   population3.Combine(Population()[num_sub1+num_sub2+i]);
			   }
			  set.clear();
			   std::vector<double> mC3,mC1; mC3.resize(P().XSize()); mC1.resize(P().XSize());
			    Check(population1, Population());
			    Population().Combine(population1);
				Check(population2, Population());
			    Population().Combine(population2);
				Check(old_memory, Population());				
				Population().Combine(old_memory);
				Population().Evaluate();
				Population().RankSort();
				for(k=0; k<xdim; k++) 
				{
					num_sub1=0;
					mC1[k]  = 0.0;
					for(i=0; i<psize; i++)
						if(Population()[i].Rank()==1)
						{ 
							mC1[k] += Population()[i][k]; 
							num_sub1++;
						}
					mC1[k] /= double(num_sub1);
				}               
				for(k=0; k<xdim; k++) 
				{
					mC3[k]  = 0.0;
					for(i=0; i<population3.Size(); i++)		mC3[k] += population3[i][k]; 							
					mC3[k] /= double(population3.Size());
				}
			   for (i = 0; i < population3.Size(); i++)
			   {
				   for (j = 0; j < population3[i].P().XSize(); j++)
				   {
					   population3[i][j]=(rnd::rand(P().XLow(j),P().XUpp(j));
					   if(population3[i][j]>P().XUpp(j)) 
							population3[i][j] = P().XUpp(j) - (population3[i][j]-P().XUpp(j))/10.0;
						else
							if(population3[i][j]<P().XLow(j)) 
								population3[i][j] = P().XLow(j) + (P().XLow(j)-population3[i][j])/10.0;
				   }
				   
			   }

				Check(population3, Population());
				Population().Combine(population3);		
				Population().Evaluate();
				population1.Clear();
				population2.Clear();
				population3.Clear();               			
			}

        } //namespace dea
	
     } //namespace mea

} //namespace az
