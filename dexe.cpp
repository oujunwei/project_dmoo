#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "AlgD.h"
#include "emo/Parameter.h"
#include "emo/Config.h"
#include <algorithm>
#include "stdlib.h"


#include <yvals.h>
#if (_MSC_VER >= 1600)
#define __STDC_UTF_16__
#endif
#include"engine.h"
#include"memory.h"
#include"../evaluate.h"
#pragma comment(lib,"libmat.lib") 
#pragma comment(lib,"libmx.lib") 
#pragma comment(lib,"libeng.lib") 
void average(unsigned int start,unsigned int ccc,unsigned int nt,unsigned int run,unsigned int predict,unsigned int metrics,
			 unsigned int testNumber,std::string instances[]);
int main(int argc, char* argv[])
{
	int ccc;
	int gh;
	//unsigned testNumber;
	unsigned int predict;	
	unsigned int strategy, runs, generation, popsize, dimension, gen, ir, i, j;
	unsigned int torder, taot, nt; 
	double		 t0, alpha;
	unsigned int run,metrics,testNumber;
	unsigned int start;
	run=20;//算法独立运行的次数
	metrics=5;//评价指标的个数

	
	//DMOPA-DMOPF : F5 F6 F7 F9 F10 F8  (7 8 9 10 11 12 )
	std::string instances[]  = {"FDA1","FDA2","FDA3","FDA4","DMOP1","DMOP2","DMOP3","DMOPA","DMOPB","DMOPC","DMOPD","DMOPE","DMOPF",
	                             "JY1","JY2","JY3","JY4","JY5","JY6","JY7","JY8","JY9"}; 
	std::string instances2[]  = {"SDMOP1","SDMOP2","SDMOP3","SDMOP4","SDMOP5","SDMOP6","SDMOP7","SDMOP8"};
	std::string instances3[]  = {"JY1","JY2","JY3","JY4","JY5","JY6","JY7","JY8","JY9"};
	std::string arraymethod[]  = {"PCA","GTM","NSDE","NSGA","RM2","MOEAD","DEE","DBEA","SBX_DE"}; 
	int arraystrategy[] = {1 , 2 , 3, 4, 5 , 6 ,7, 8,9,10};
	

	start=0;  //从第几个测试问题开始
	testNumber=1;;//测试问题的个数
	for(gh = start;gh<testNumber;gh++){
		ccc = 0;
		/*unsigned run,metrics,testNumber;*/
		
	while(ccc<run)
	{
	std::string method, problem, fname, path;
	az::mea::CParameter mPar;
	char savefilename[1024];
	
	problem = instances[gh];
	method = arraymethod[7];
	strategy = arraystrategy[8];
	runs = 1; //Number of runs
	popsize = 100;
	alpha=1;
	taot = 20;
	t0 = 0;
	nt =10;
	torder = 3;
	dimension =10;

	int a[21] = {10,14,0,19,11,16,13,5,2,8,17,12,15,4,20,9,6,1,3,7,18};
	std::vector<int> t_order;
	for(i=0;i<2*nt+1;i++) 
		t_order.push_back(a[i]);

	sprintf(savefilename,"data/%s_%s_%d_%d_%d",method.c_str(),problem.c_str(),strategy,taot,nt);
	mPar.TolF() = 1.0E-5;//TOLERANCEF
	mPar.TolX() = 1.0E-5;//TOLERANCEX
	mPar.TolC() = 1.0E-5;//TOLERANCEC

	mPar.XSize(dimension);
	mPar.Problem(problem);
	mPar.XCoding() = false;

	mPar.Pc() = 0.9;
	mPar.Pm() = 0.1;
	mPar.ETA_PM() = 20.0;
	mPar.ETA_SBX() = 20.0;

	std::vector< std::vector<double> > PF0, PF1;
	std::vector< std::vector<double> > POF;
	unsigned int ic,ipf0,ipf1;

	bool justinit = true;

	std::cout<<savefilename<<std::endl;

	std::ofstream f0,f1,pf;

	generation = 1000;
	
	/*unsigned int predict;*/
	predict=generation/taot;
	srand(123456);
	std::random_shuffle(t_order.begin(),t_order.end());
    az::mea::dea::DMOO* pEA = new az::mea::dea::DMOO(strategy, method, popsize, generation, taot, nt, torder, t0, alpha, mPar,t_order);
	unsigned int xdim = (dimension < 3) ? dimension:3;
	xdim = dimension ;
	bool flag =false ;
	for(ir=0; ir<runs; ir++)
	{
		ic = ipf1 = ipf0 = 0;

		PF0.resize(popsize*20*nt);
		PF1.resize(popsize*20*nt);
		pEA->Reset();

		unsigned int cou = 0;
		while(!pEA->IsTerminate())
		{
			Engine *ep;
			gen = pEA->Step();
			if(justinit)
			{
				justinit = false;
				for(i=0; i<pEA->Population().Size(); i++) 
				{
					PF1[ipf1+i].resize(pEA->P().FSize()+xdim);
					for(j=0; j<pEA->P().FSize(); j++) PF1[ipf1+i][j]		= pEA->Population()[i].F(j);
					for(j=0; j<xdim; j++) PF1[ipf1+i][j+pEA->P().FSize()]	= pEA->Population()[i][j];
				}
				ipf1 += pEA->Population().Size();
			}

			if(pEA->IsToChange())
			{
				char pfname[1024];
				POF.resize(pEA->Population().Size());
				for(i=0; i<pEA->Population().Size(); i++) 
				{
					PF0[ipf0+i].resize(pEA->P().FSize()+xdim);
					POF[i].resize(pEA->P().FSize());
					for(j=0; j<pEA->P().FSize(); j++)
					{
						PF0[ipf0+i][j]		= pEA->Population()[i].F(j);
						POF[i][j]      = pEA->Population()[i].F(j);
					}
					for(j=0; j<xdim; j++) PF0[ipf0+i][j+pEA->P().FSize()]	= pEA->Population()[i][j];					
				}
				if (pEA->Population().P().FSize()==2)
				{
					double f1[100];
					double f2[100];
					double mNt = nt;
					mxArray *T = mxCreateDoubleMatrix( 1 ,  pEA->Population().Size() , mxREAL  ) ; 
					mxArray *T2 = mxCreateDoubleMatrix( 1 ,  pEA->Population().Size() , mxREAL  ) ;
					mxArray *M = mxCreateDoubleMatrix( 1 ,  1 , mxREAL  ) ;
					for(int n=0; n<pEA->Population().Size(); n++)
					{ 
						f1[n]=pEA->Population()[n].F(0);
						f2[n]=pEA->Population()[n].F(1);
					}
					if(  (ep=engOpen(NULL)) )
					{
						memcpy( (char*)mxGetPr(T),(char*)f1 , pEA->Population().Size()*sizeof(double));
						memcpy( (char*)mxGetPr(T2),(char*)f2 , pEA->Population().Size()*sizeof(double));
						memcpy( (char*)mxGetPr(M),(char*)&mNt ,1*sizeof(double));
						engPutVariable(ep , "TX", T) ;
						engPutVariable(ep , "TY", T2) ;
						engPutVariable(ep , "NT", M) ;
						if(flag==false)
						{
							flag =true ;
							engEvalString( ep ,"h=plot(TX,TY,'ro');grid on");
							if(problem=="FDA1"||problem=="DMOP3"||problem=="SDMOP1")
							{
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^0.5 ;"); 
								engEvalString(ep,"hold on");
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");						   
							}

							if(problem=="JY1")
							{
								engEvalString(ep,"x=0:0.01:1 ;frontx=x+0.05*sin(6*pi*x) ;fronty=1-x+0.05*sin(6*pi*x) ;"); 
								engEvalString(ep,"hold on");
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");						   
							}

							if(problem=="FDA2")
							{
								engEvalString(ep,"hold off");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1.75*1.75)*10+0.75)) ;"); 
								engEvalString(ep,"hold on");
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");	
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*1/NT)))*(1+(0.75+0.7*sin(0.5*pi*1/NT)))*10+(0.75+0.7*sin(0.5*pi*1/NT)) )) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");	
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*2/NT)))*(1+(0.75+0.7*sin(0.5*pi*2/NT)))*10+(0.75+0.7*sin(0.5*pi*2/NT)) )) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");	
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*3/NT)))*(1+(0.75+0.7*sin(0.5*pi*3/NT)))*10+(0.75+0.7*sin(0.5*pi*3/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");	
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*4/NT)))*(1+(0.75+0.7*sin(0.5*pi*4/NT)))*10+(0.75+0.7*sin(0.5*pi*4/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");	
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*5/NT)))*(1+(0.75+0.7*sin(0.5*pi*5/NT)))*10+(0.75+0.7*sin(0.5*pi*5/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");	
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*11/NT)))*(1+(0.75+0.7*sin(0.5*pi*11/NT)))*10+(0.75+0.7*sin(0.5*pi*11/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*12/NT)))*(1+(0.75+0.7*sin(0.5*pi*12/NT)))*10+(0.75+0.7*sin(0.5*pi*12/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*13/NT)))*(1+(0.75+0.7*sin(0.5*pi*13/NT)))*10+(0.75+0.7*sin(0.5*pi*13/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*14/NT)))*(1+(0.75+0.7*sin(0.5*pi*14/NT)))*10+(0.75+0.7*sin(0.5*pi*14/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*15/NT)))*(1+(0.75+0.7*sin(0.5*pi*15/NT)))*10+(0.75+0.7*sin(0.5*pi*15/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*21/NT)))*(1+(0.75+0.7*sin(0.5*pi*21/NT)))*10+(0.75+0.7*sin(0.5*pi*21/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*22/NT)))*(1+(0.75+0.7*sin(0.5*pi*22/NT)))*10+(0.75+0.7*sin(0.5*pi*22/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*23/NT)))*(1+(0.75+0.7*sin(0.5*pi*23/NT)))*10+(0.75+0.7*sin(0.5*pi*23/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*24/NT)))*(1+(0.75+0.7*sin(0.5*pi*24/NT)))*10+(0.75+0.7*sin(0.5*pi*24/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*25/NT)))*(1+(0.75+0.7*sin(0.5*pi*25/NT)))*10+(0.75+0.7*sin(0.5*pi*25/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*31/NT)))*(1+(0.75+0.7*sin(0.5*pi*31/NT)))*10+(0.75+0.7*sin(0.5*pi*31/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*32/NT)))*(1+(0.75+0.7*sin(0.5*pi*32/NT)))*10+(0.75+0.7*sin(0.5*pi*32/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*33/NT)))*(1+(0.75+0.7*sin(0.5*pi*33/NT)))*10+(0.75+0.7*sin(0.5*pi*33/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*34/NT)))*(1+(0.75+0.7*sin(0.5*pi*34/NT)))*10+(0.75+0.7*sin(0.5*pi*34/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1/((1+(0.75+0.7*sin(0.5*pi*35/NT)))*(1+(0.75+0.7*sin(0.5*pi*35/NT)))*10+(0.75+0.7*sin(0.5*pi*35/NT)) ))  ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
							}
							if(problem=="DMOP1"||problem=="DMOP2")
							{
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^1.25 ;"); 
								engEvalString(ep,"hold on");
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");		
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*1/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*2/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*3/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*4/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*5/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*6/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*7/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*8/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*9/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*10/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");

								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*21/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*22/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*23/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*24/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*25/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*26/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*27/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*28/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*29/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.25+0.75*sin(0.5*pi*30/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
							}
							if(problem=="DMOPA"||problem=="DMOPB"||problem=="DMOPC"||problem=="DMOPD"||problem=="DMOPE")
							{
								engEvalString(ep,"m=0:0.01:1,frontx=m.^1.25 ;fronty=(1-m).^1.25 ;"); 
								engEvalString(ep,"hold on");
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");		
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(1.25+0.75*sin(pi/NT)) ;fronty=(1-m).^(1.25+0.75*sin(pi/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");		
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(1.25+0.75*sin(pi*2/NT)) ;fronty=(1-m).^(1.25+0.75*sin(pi*2/NT)) ;");  
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(1.25+0.75*sin(pi*3/NT)) ;fronty=(1-m).^(1.25+0.75*sin(pi*3/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(1.25+0.75*sin(pi*4/NT)) ;fronty=(1-m).^(1.25+0.75*sin(pi*4/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(1.25+0.75*sin(pi*5/NT)) ;fronty=(1-m).^(1.25+0.75*sin(pi*5/NT)) ;");
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(1.25+0.75*sin(pi*11/NT)) ;fronty=(1-m).^(1.25+0.75*sin(pi*11/NT)) ;");
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(1.25+0.75*sin(pi*12/NT)) ;fronty=(1-m).^(1.25+0.75*sin(pi*12/NT)) ;");
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(1.25+0.75*sin(pi*13/NT)) ;fronty=(1-m).^(1.25+0.75*sin(pi*13/NT)) ;");
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(1.25+0.75*sin(pi*14/NT)) ;fronty=(1-m).^(1.25+0.75*sin(pi*14/NT)) ;");
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(1.25+0.75*sin(pi*15/NT)) ;fronty=(1-m).^(1.25+0.75*sin(pi*15/NT)) ;");
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
							}

							if(problem=="SDMOP2"||problem=="SDMOP3")
							{
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^2 ;"); 
								engEvalString(ep,"hold on");
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");		
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*1/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*2/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*3/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*4/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*5/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*6/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*7/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*8/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*9/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*10/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");

								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*11/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*12/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*13/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*14/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*15/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*16/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*17/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*18/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*19/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=1-frontx.^(1.2+0.8*cos(0.5*pi*20/NT)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
							}
							if(problem == "SDMOP4")
							{
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=(1+1)*(1-frontx-cos(10*pi*frontx+0.5*pi)/(10*pi)) ;"); 
								engEvalString(ep,"hold on");
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=(1+cos(0.5*pi*1/NT))*(1-frontx-cos(10*pi*frontx+0.5*pi)/(10*pi)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=(1+cos(0.5*pi*2/NT))*(1-frontx-cos(10*pi*frontx+0.5*pi)/(10*pi)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=(1+cos(0.5*pi*3/NT))*(1-frontx-cos(10*pi*frontx+0.5*pi)/(10*pi)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=(1+cos(0.5*pi*4/NT))*(1-frontx-cos(10*pi*frontx+0.5*pi)/(10*pi)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=(1+cos(0.5*pi*5/NT))*(1-frontx-cos(10*pi*frontx+0.5*pi)/(10*pi)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=(1+cos(0.5*pi*6/NT))*(1-frontx-cos(10*pi*frontx+0.5*pi)/(10*pi)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=(1+cos(0.5*pi*7/NT))*(1-frontx-cos(10*pi*frontx+0.5*pi)/(10*pi)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=(1+cos(0.5*pi*8/NT))*(1-frontx-cos(10*pi*frontx+0.5*pi)/(10*pi)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=(1+cos(0.5*pi*9/NT))*(1-frontx-cos(10*pi*frontx+0.5*pi)/(10*pi)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
								engEvalString(ep,"frontx=0:0.01:1 ;fronty=(1+cos(0.5*pi*10/NT))*(1-frontx-cos(10*pi*frontx+0.5*pi)/(10*pi)) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'Linewidth',2)");
							}

							if(problem == "FDA3")
							{
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*0/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*0/NT)))))*(1+abs(sin(0.5*pi*0/NT))) ;"); 
								engEvalString(ep,"hold on");
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*1/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*1/NT)))))*(1+abs(sin(0.5*pi*1/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*2/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*2/NT)))))*(1+abs(sin(0.5*pi*2/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*3/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*3/NT)))))*(1+abs(sin(0.5*pi*3/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*4/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*4/NT)))))*(1+abs(sin(0.5*pi*4/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*5/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*5/NT)))))*(1+abs(sin(0.5*pi*5/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*11/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*11/NT)))))*(1+abs(sin(0.5*pi*11/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*12/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*12/NT)))))*(1+abs(sin(0.5*pi*12/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*13/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*13/NT)))))*(1+abs(sin(0.5*pi*13/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*14/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*14/NT)))))*(1+abs(sin(0.5*pi*14/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*15/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*15/NT)))))*(1+abs(sin(0.5*pi*15/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");

								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*21/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*21/NT)))))*(1+abs(sin(0.5*pi*21/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*22/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*22/NT)))))*(1+abs(sin(0.5*pi*22/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*23/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*23/NT)))))*(1+abs(sin(0.5*pi*23/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*24/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*24/NT)))))*(1+abs(sin(0.5*pi*24/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*25/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*25/NT)))))*(1+abs(sin(0.5*pi*25/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*31/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*31/NT)))))*(1+abs(sin(0.5*pi*31/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*32/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*32/NT)))))*(1+abs(sin(0.5*pi*32/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*33/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*33/NT)))))*(1+abs(sin(0.5*pi*33/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*34/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*34/NT)))))*(1+abs(sin(0.5*pi*34/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
								engEvalString(ep,"m=0:0.01:1,frontx=m.^(10^(2*sin(0.5*pi*35/NT))) ;fronty=(1-sqrt(frontx/(1+abs(sin(0.5*pi*35/NT)))))*(1+abs(sin(0.5*pi*35/NT))) ;"); 
								engEvalString(ep,"plot(frontx,fronty,'b.')");
							}
						}
						else
						{
							engEvalString( ep ,"set(h, 'XData',TX     ,'YData',TY   )");
						}
					}
				}

				if (pEA->Population().P().FSize()==3)
				{
					double f1[200];
					double f2[200];
					double f3[200];
					double mNt = nt;
					mxArray *T = mxCreateDoubleMatrix( 1 ,  pEA->Population().Size() , mxREAL  ) ; 
					mxArray *T2 = mxCreateDoubleMatrix( 1 ,  pEA->Population().Size() , mxREAL  ) ;
					mxArray *T3 = mxCreateDoubleMatrix( 1 ,  pEA->Population().Size() , mxREAL  ) ;
					mxArray *M = mxCreateDoubleMatrix( 1 ,  1 , mxREAL  ) ;
					for(int n=0; n<pEA->Population().Size(); n++)
					{ 
						f1[n]=pEA->Population()[n].F(0);
						f2[n]=pEA->Population()[n].F(1);
						f3[n]=pEA->Population()[n].F(2);
					}
					if(  (ep=engOpen(NULL)) )
					{
						memcpy( (char*)mxGetPr(T),(char*)f1 , pEA->Population().Size()*sizeof(double));
						memcpy( (char*)mxGetPr(T2),(char*)f2 , pEA->Population().Size()*sizeof(double));
						memcpy( (char*)mxGetPr(T3),(char*)f3 , pEA->Population().Size()*sizeof(double));
						memcpy( (char*)mxGetPr(M),(char*)&mNt ,1*sizeof(double));
						engPutVariable(ep , "TX", T) ;
						engPutVariable(ep , "TY", T2) ;
						engPutVariable(ep , "TZ", T3) ;
						engPutVariable(ep , "NT", M) ;
						if(flag==false)
						{
							flag =true ;
							engEvalString( ep ,"h=plot3(TX,TY,TZ,'r.')");
							if (problem == "FDA4"|| problem == "DMOPF")
							{
								engEvalString(ep,"t=linspace(0,pi/2,25);p=linspace(0,pi/2,25);[theta,phi]=meshgrid(t,p);x=cos(theta).*cos(phi);y=cos(theta).*sin(phi);z=sin(theta);"); 
								engEvalString(ep,"hold on");
								engEvalString(ep,"plot3(x,y,z,'bo');grid on");		
							}

							if(problem == "SDMOP5" || problem == "SDMOP6" || problem == "SDMOP7")
							{
								engEvalString(ep,"[t1 t2]=meshgrid(linspace(0,1,20),linspace(0,1,20)) ;"); 
								engEvalString(ep,"hold on");
							    engEvalString(ep,"x=cos(0.5*pi*t1).*cos(0.5*pi*t2);  y=cos(0.5*pi*t1).*sin(0.5*pi*t2);  z=sin(0.5*pi*t1);"); 
								engEvalString(ep,"alpha(surf(x,y,z),0.8);   re=[0,1,0];  colormap(re); view(69,18);");
							}
							if(problem == "SDMOP8")
							{
								/*engEvalString(ep,"[t1 t2]=meshgrid(linspace(0,1,20),linspace(0,1,20)) ;"); 
								engEvalString(ep,"hold on");
							    engEvalString(ep,"x=t1;  y=t2;  z=3-0.5*x.*x.*(1+sin(3*pi*x))-0.5*y.*y.*(1+sin(3*pi*y));));"); 
								engEvalString(ep,"alpha(surf(x,y,z),0.8);   re=[0,1,0];  colormap(re); view(69,18);"); */
							}

						}
						else
						{
							engEvalString( ep ,"set(h, 'XData',TX     ,'YData',TY  ,'ZData',TZ  )");

						}
					}
				}
				engEvalString(ep,"hold off");
				ipf0 += pEA->Population().Size();
				justinit = true;
				std::cout<<"T"<<gen<<"\n";
				sprintf(pfname,"PF/pf_%s_%d_%d",problem.c_str(),ccc,gen/taot);  
				std::stringstream spf;
				spf<<pfname<<".dat";
				pf.open(spf.str().c_str());
				pf<<std::scientific<<std::setprecision(5);
				for(i=0; i<pEA->Population().Size(); i++)
				{
					for(j=0; j<pEA->P().FSize(); j++) pf<<POF[i][j]<<"\t";
					pf<<std::endl;
				}
				POF.clear();
				pf.close();
			}
			cou++;
		}
		std::stringstream ss0,ss1;
		ss1<<savefilename<<"_"<<ir<<".pin";
		f1.open(ss1.str().c_str());
		f1<<std::scientific<<std::setprecision(5);
		for(i=0; i<ipf1; i++)
		{
			for(j=0; j<pEA->P().FSize()+xdim; j++) f1<<PF1[i][j]<<"\t";
			f1<<std::endl;
		}
		f1.close();

		ss0<<savefilename<<"_"<<ir<<".pop";
		f0.open(ss0.str().c_str());
		f0<<std::scientific<<std::setprecision(5);
		for(i=0; i<ipf0; i++)
		{
			for(j=0; j<pEA->P().FSize()+xdim; j++) f0<<PF0[i][j]<<"\t";
			f0<<std::endl;
		}
		f0.close();

		PF0.clear();
		PF1.clear();
	}
	ccc++;
	}
  }
  std::cout<<std::endl;
    average(start,ccc,nt,run,predict,metrics,testNumber,instances);
  	return 1;
}


void average(unsigned int start,unsigned int ccc,unsigned int nt,unsigned int run,unsigned int predict,unsigned int metrics,
			 unsigned int testNumber,std::string instances[]){

	std::ofstream pf;
	std::ofstream sta;
	ccc=0;
	std::cout<<predict;
	while(ccc<run)
	{
		evaluate eva;
		
		vector<double> gd;
		char filename1[1024];
		char filename2[1024];
		char    strTestInstance[256];
		for (int i=start;i<testNumber;i++)
		{			
			if(instances[i]=="JY1"||instances[i]=="JY2"||instances[i]=="JY3"||instances[i]=="JY4"||instances[i]=="JY5"
				||instances[i]=="JY6"||instances[i]=="JY7"||instances[i]=="JY8"||instances[i]=="JY9"){
				char pfname[1024];
				sprintf(pfname,"evaluate/data/%s_%d.dat",instances[i].c_str(),ccc);

				pf.open(pfname,ios::trunc);

				double ave=0.0,fangcha = 0.0;
				int gen=0;
				for(gen=1; gen<=predict;gen++)
				{
					sprintf(strTestInstance,"%s",instances[i].c_str());
					eva.getPOF(instances[i],nt,gen,eva.pf);  //取真实的点

					sprintf(filename2,"PF/pf_%s_%d_%d.dat",strTestInstance,ccc,gen);
					eva.loadpfront(filename2,eva.p,2);
                    pf<<eva.indicator_GD()<<"               ";
					pf<<eva.indicator_IGD()<<"               ";
					pf<<eva.indicator_SP(eva.p)<<"               ";
					pf<<eva.indicator_MS(eva.pf,eva.p,2)<<"               ";
					pf<<eva.indicator_hvd(eva.p,gen,instances[i].c_str(),nt)<<"\n";
					cout<<"ccc = "<<ccc<<"    "<<instances[i]<<"  gen = "<<gen<<"  "<<endl;
					eva.p.clear();
					eva.pf.clear();
				}
				pf.close();
			}
		
			if (instances[i]=="FDA1"||instances[i]=="FDA2"||instances[i]=="FDA3"||instances[i]=="DMOP1"||instances[i]=="DMOP2"||instances[i]=="DMOP3"||
				instances[i]=="DMOPA"||instances[i]=="DMOPB"||instances[i]=="DMOPC"||instances[i]=="DMOPD"||instances[i]=="DMOPE")
			{
				char pfname[1024];
				int gen=0;
				double ave=0.0;
				sprintf(pfname,"evaluate/data/%s_%d.dat",instances[i].c_str(),ccc);

				//pf.open(pfname);
				pf.open(pfname,ios::trunc);
				for(gen=1;gen<=predict;gen++)
				{
					sprintf(strTestInstance,"%s",instances[i].c_str());
					eva.getPOF(instances[i],nt,gen,eva.pf);  //取真实的点

					sprintf(filename2,"PF/pf_%s_%d_%d.dat",strTestInstance,ccc,gen);									
					eva.loadpfront(filename2,eva.p,2);
					
					pf<<eva.indicator_GD()<<"               ";
					pf<<eva.indicator_IGD()<<"               ";
					pf<<eva.indicator_SP(eva.p)<<"               ";
					pf<<eva.indicator_MS(eva.pf,eva.p,2)<<"               ";
					pf<<eva.indicator_hvd(eva.p,gen,instances[i],nt)<<"\n";
					cout<<"ccc = "<<ccc<<"    "<<instances[i]<<"  gen = "<<gen<<"  "<<endl;;
				
					eva.p.clear();
					eva.pf.clear();
				}
				pf.close();
			}

			if(instances[i]=="DMOPF" ||instances[i]=="FDA4" ||instances[i]=="FDA5")
			{
				char pfname[1024];
				double ave=0.0;
				int gen=0;
				sprintf(pfname,"evaluate/data/%s_%d.dat",instances[i].c_str(),ccc);

				pf.open(pfname,ios::trunc);

				for(gen=1; gen<=predict;gen++)
				{
					sprintf(strTestInstance,"%s",instances[i].c_str());
					
					
					eva.getPOF(instances[i],nt,gen,eva.pf);  //取真实的点

					sprintf(filename2,"PF/pf_%s_%d_%d.dat",strTestInstance,ccc,gen);										
					eva.loadpfront(filename2,eva.p,3);

					pf<<eva.indicator_GD()<<"               ";
					pf<<eva.indicator_IGD()<<"               ";
					pf<<eva.indicator_SP(eva.p)<<"               ";
					pf<<eva.indicator_MS(eva.pf,eva.p,2)<<"               ";
					pf<<eva.indicator_hvd(eva.p,gen,instances[i],nt)<<"\n";
					cout<<"ccc = "<<ccc<<"    "<<instances[i]<<"  gen = "<<gen<<"  "<<endl;			
					eva.p.clear();
					eva.pf.clear();
				}			
				pf.close();
			}
			
		}
		ccc++;
	}
	

	//求样本平均值和样本方差
	int kkk=0,iii;
	evaluate eva;
	for (int kkk = start ; kkk < testNumber; kkk++)
	{
		char pfname[1024];
		
		double ave=0.0;
		int t=0;
		sprintf(pfname,"evaluate/avg/%s.dat",instances[kkk].c_str());
		
		pf.open(pfname,ios::trunc);

		char    strTestInstance[256];
		char    strTestInstanceAv[256];

		cout<<instances[kkk]<<endl;

		pf<<setprecision(2)<<setiosflags(ios::scientific);
		for(ccc=0;ccc<run;ccc++){     //每次独立运行的求平均
			sprintf(strTestInstance,"evaluate/data/%s_%d.dat",instances[kkk].c_str(),ccc);
			eva.loadpfront(strTestInstance,eva.pf,metrics);
			eva.indicator_AVG(eva.pf,eva.in_avg);
			for(iii=0;iii<eva.in_avg.size();iii++){
				pf<<eva.in_avg[iii]<<"      ";
			}
			pf<<"\n";	
		}
		pf.close();

		//独立run次的平均和方差
		pf.open(pfname,ios::app);
		pf<<setprecision(5)<<setiosflags(ios::scientific);
		sprintf(strTestInstanceAv,"evaluate/avg/%s.dat",instances[kkk].c_str());
		eva.loadpfront(strTestInstanceAv,eva.pf,metrics);
		eva.indicator_AVG(eva.pf,eva.in_avg);
		for(iii=0;iii<eva.in_avg.size();iii++){
				pf<<eva.in_avg[iii]<<"      ";
		}
		pf<<"\n";
		int size1=eva.pf.size();
		int size2=eva.pf[0].size();
		vector<double> values;
		for(iii=0;iii<size2;iii++){
			for(ccc=0;ccc<size1;ccc++){
				values.push_back(eva.pf[ccc][iii]);
			}
			pf<<eva.STD(values)<<"      ";
			values.clear();
		}
		pf<<"\n";
		pf.close();

		//求平均IGD
		eva.GetAvgIGD(instances[kkk],run,predict);
		eva.Statistics(instances[kkk],metrics,run);//统计数据
	}
}
