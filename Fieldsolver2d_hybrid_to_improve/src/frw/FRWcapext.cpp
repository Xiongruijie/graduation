#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <fstream>
#include "FRWcapext.h"
#include <unistd.h>
using namespace std;



bool FRWControl::updateProgress(int tempWalkNum, int masterCondIndex,  double temp_VGSArea ){
	currentNumOfWalkOfCond[masterCondIndex]+=tempWalkNum;

	for(int index=0; index< FRWControl::capMatrix.size(); index++ ){
		double average_cap=FRWControl::capMatrix[index].sumOfCap/currentNumOfWalkOfCond[masterCondIndex];
		double variance_cap=FRWControl::capMatrix[index].sumOfCapSquared/currentNumOfWalkOfCond[masterCondIndex] - pow(average_cap,2.0);//?
		double std_var_cap=pow( variance_cap/currentNumOfWalkOfCond[masterCondIndex], 0.5  );
		double cap_err=100.0*3*std_var_cap/average_cap;
		std_var_cap*=temp_VGSArea;
		FRWControl::capMatrix[index].estimateErr=cap_err;
		FRWControl::capMatrix[index].std_var=std_var_cap;
	}


	int temp=currentNumOfWalkOfCond[masterCondIndex]*100 / totalNumOfWalkOfCond[masterCondIndex] ;
	if( temp> currentProgressOfCond[masterCondIndex]  ){
		currentProgressOfCond[masterCondIndex] =temp;
		// cout<<"\rConductor No "<<masterCondIndex<<"/"<<FRWControl::conductorList.size()-1<<":\tFRW Path number"<<temp<<"%"<<flush;  //<<"%\tstd_err: "<<std_var_cap<<"\t percent: "<<cap_err<<"%"<<endl;
	}

	if( currentProgressOfCond[masterCondIndex]>=100 ){
		return true;
	}else{

		for(int index=0; index< FRWControl::capMatrix.size(); index++ ){
			if( FRWControl::capMatrix[index].estimateErr>2 ){
				return false;
			}
		}
		return true;
	}

}

Face::Face(){
	x1=x2=z1=z2=0;
}

Face::Face(double x_1, double z_1, double size ){
	x1=x_1;
	z1=z_1;
	x2=x_1+size;
	z2=z_1+size;
}

Face::Face( FPoint2d &minPoint, double size ){
	x1=minPoint.x1;
	z1=minPoint.z1;
	x2=x1+size;
	z2=z1+size;
}

Face::Face( double x_1, double x_2,  double z_1, double z_2 ){
	x1=x_1;
	x2=x_2;
	z1=z_1;
	z2=z_2;
}

double Face::size(){
	return x2-x1;
}

GaussianSurface2d::GaussianSurface2d(double x_1, double x_2, double z_1, double z_2 ,std::vector<double> extensionDistance ){
	x1=x_1;
	x2=x_2;
	z1=z_1;
	z2=z_2;
	extensionDis=extensionDistance;
	setSurfaceArea();
}

void GaussianSurface2d::setSurfaceArea(){
	surfaceArea=std::vector<double> (4);
	surfaceArea[pX]= (z2-z1);
	surfaceArea[pZ]= (x2-x1);
	surfaceArea[nX]= (z2-z1);
	surfaceArea[nZ]= (x2-x1);
	// surfaceArea=std::vector<double> (1);
	// surfaceArea[pX]=surfaceArea[pX] = surfaceArea[pX] = surfaceArea[pX] = ( x2-x1 ) *(z2-z1);

}

double GaussianSurface2d::area(){
	double len=x2-x1;
	double width=z2-z1;

	// return len*width;
	return 2*(len + width);
}

bool GaussianSurface2d::isPointEnclosedBySurface( FPoint2d &point ){
	return  point.x1>x1 && point.x1<x2 && point.z1 >z1 && point.z1 <z2;
}

bool GaussianSurface2d::isPointOnSurface(FPoint2d &point){
	if(  point.x1 >=x1 && point.x1 <=x2  ){
		if( fabs( point.z1-z1)<1E-6  || fabs( point.z1-z2 )<1E-6   ){
			return true;
		}
	}
	if(  point.z1 >=z1 && point.z1 <=z2   ){
		if(  fabs(point.x1-x1 )< 1E-6 || fabs( point.x1-x2 )< 1E-6  ){
			return true;
		}
	}
	return false;

}

int GaussianSurface2d::getNormalDirOfPoint(FPoint2d  &point){
	if(  point.z1 >=z1 && point.z1 <=z2   ){
		if(  fabs(point.x1-x1 )< 1E-7 ){
			return -1;
		}else if(fabs( point.x1-x2 )<1E-7 ){
			return 1;
		}
	}
	if(  point.x1 >=x1 && point.x1 <=x2 ){
		if( fabs( point.z1-z1)<1E-7 ){
			return -3;
		}else if(fabs( point.z1-z2 )<1E-7   ){
			return 3;
		}
	}
	//point is not on the surface
	return 0;
}

bool GaussianSurface2d::operator ==( GaussianSurface2d &rhs  ){
	if(  x1==rhs.x1&&z1==rhs.z1 &&x2==rhs.x2 &&z2==rhs.z2  ){
		return true;
	}
	return false;
}


bool GaussianSurface2d::isNeighborWith( GaussianSurface2d &gaussian ){
	return x1<=gaussian.x2 && x2>=gaussian.x1 && z1<=gaussian.z2 && z2 >= gaussian.z1;
}


void *FRW( void *condIn ){

	long randSeed_uniform_01=0;
	long randSeed_int_normalDir=0;
	long randSeed_int_GVT=0;  //green function value table

	pthread_mutex_lock(&getRandSeed_mutex);

	int threadIndex=FRWControl::currentThreadIndex++;
	if(FRWControl::randSeedVec.size()==0 ){
		cout<<"Err: no randSeed in randSeedVec!"<<endl;
		exit(1);
	}
	randSeed_uniform_01=FRWControl::randSeedVec.back();
	FRWControl::randSeedVec.pop_back();

	if(FRWControl::randSeedVec.size()==0 ){
		cout<<"Err: no randSeed in randSeedVec!"<<endl;
		exit(1);
	}
	randSeed_int_normalDir = FRWControl::randSeedVec.back();
	FRWControl::randSeedVec.pop_back();	

	if(FRWControl::randSeedVec.size()==0 ){
		cout<<"Err: no randSeed in randSeedVec!"<<endl;
		exit(1);
	}
	randSeed_int_GVT =FRWControl::randSeedVec.back();
	FRWControl::randSeedVec.pop_back();
	pthread_mutex_unlock(&getRandSeed_mutex);

	//@author:Dragon
	//initialize boost random generator, use randSeed!
	boost::mt19937 gen_uni_01(randSeed_uniform_01);
	boost::mt19937 gen_int_normalDir(randSeed_int_normalDir );
	boost::mt19937 gen_int_GVT(randSeed_int_GVT );
	//here we use a group of rand seed generator!
	boost::uniform_01<> dist_uni_01;
	boost::uniform_int<> dist_int_normalDir(0,3);
	boost::uniform_int<> dist_int_GVT(0, FRWControl::GreenVT.size()-1);

	boost::variate_generator<boost::mt19937&,boost::uniform_01<> > rand_01(gen_uni_01, dist_uni_01 );
	boost::variate_generator<boost::mt19937&,boost::uniform_int<> >rand_normalDir(gen_int_normalDir, dist_int_normalDir );
	boost::variate_generator<boost::mt19937&,boost::uniform_int<> > rand_GVT(gen_int_GVT, dist_int_GVT );

	FPoint2d pointOnVGS2d, nextPoint2d, point1d;
	Face face;
	Cell2d *rootCell2d=NULL, *leafCell2d=NULL;
	bool isFinish=false;
	int loopTimer=0;
	int Nt=0, Ng=0;
	//生成本地电容矩阵
	CapacitanceMatrix capMatrix=generateCapMatrix( FRWControl::teatlist );
	list<OneConduct> ::iterator condIt;
	condIt = FRWControl::teatlist.testlist.begin();
    
	while( true  )
	{ 			
		FloatingRandomWalk( pointOnVGS2d, face, rootCell2d, leafCell2d, nextPoint2d,  point1d,  *condIt , FRWControl::teatlist ,  
						FRWControl::gridOctree2d, capMatrix ,FRWControl::GreenVT,FRWControl::GExVT,FRWControl::GEzVT ,rand_01, rand_normalDir, Nt, Ng);

		if(++loopTimer>=10000 ){
			isFinish=false;
			pthread_mutex_lock(&updateFRWData_mutex);
			condIt->conductor2d.gaussianSurfaceList.Nt+=Nt;
			condIt->conductor2d.gaussianSurfaceList.Ng+=Ng;
			condIt->conductor2d.gaussianSurfaceList.VGSArea= double(condIt->conductor2d.gaussianSurfaceList.Ng)/double(condIt->conductor2d.gaussianSurfaceList.Nt)* condIt->conductor2d.gaussianSurfaceList.SumOfBGSArea;  

			for(int i=0; i < capMatrix.size(); i++ ){
				FRWControl::capMatrix[i].update(capMatrix[i]);
				capMatrix[i].reset();
			}
			isFinish=FRWControl::updateProgress(loopTimer,0, condIt->conductor2d.gaussianSurfaceList.VGSArea);
			pthread_mutex_unlock(&updateFRWData_mutex);
			//note: betwwen lock and unlock there should be on jump out!!!
			loopTimer=0;
			if(isFinish){
				break;
			}
		}
	}
	pthread_exit(NULL);
}



void FloatingRandomWalk(  FPoint2d &pointOnVGS, Face &cube, Cell2d* &rootCell, Cell2d* &leafCell, FPoint2d &nextPoint,  FPoint2d &point2d, OneConduct &cond,  TestList &condList,  GridOctree2d &gridOctree,  CapacitanceMatrix &capMatrix 
	,std::vector<double> &GreenVT, vector<double> &GExVT, vector<double> &GEzVT ,
	boost::variate_generator<boost::mt19937&,boost::uniform_01<> > &rand_01, 
    boost::variate_generator<boost::mt19937&,boost::uniform_int<> > &rand_normalDir , int &Nt, int &Ng ){

	// cout<<"next FRW started"<<endl;
	int numOfStep=1;
	double PDFOnVGS;
	int normalDirOfVGS=0;
	int normalDirOfCube=0;
	// FPoint 
	// int Nt, Ng;
	double PDFIntegralElement;
	// Nt=Ng=0;
	PDFIntegralElement=0;
	sampleOnVGSS( pointOnVGS, cond.conductor2d.gaussianSurfaceList , PDFOnVGS, normalDirOfVGS ,  Nt, Ng, PDFIntegralElement, rand_01 ); ///--------------------------------------
	// Cube 
	generateMaxTransitionCube(cube,  rootCell, leafCell ,pointOnVGS,  gridOctree, normalDirOfCube );        ///-------------------
	
	generatePointOnCubeSurface(nextPoint, point2d ,cube, normalDirOfCube, GreenVT ,rand_01, rand_normalDir );   ///----------------

	// // double green on this 
	double greenFuncValue=surfaceGreenFunction( cube, nextPoint, normalDirOfCube, GreenVT  );
	// // double green on this 
	double greenFuncElectricField=greenFunctionForElectricField(  cube, nextPoint,  normalDirOfVGS, normalDirOfCube, GExVT, GEzVT );
	
	int unitNormalOfVGS;
	switch(normalDirOfVGS){
		case nX:
 		case nZ: unitNormalOfVGS=-1;break;
		case pX:
		case pZ: unitNormalOfVGS=1;break;
		default:cout<<"wrong normalDirOfVGS!"<<endl; exit(1);
	}
	// this 
	double capacitance= FRWControl::dielectricConstant  * greenFuncElectricField*unitNormalOfVGS / greenFuncValue;//?

	

	while(true){
		//if point is out of GDS::zone or on the border of GDS::zone, return.
		//Note that point on the GDS::zone is not permitted, or generateGridOctree will break down
		if( nextPoint.x1 <= BoundStatic::zone.x0 || nextPoint.x1 >=BoundStatic::zone.x1 || nextPoint.z1<=BoundStatic::zone.z0  || nextPoint.z1 >=BoundStatic::zone.z1    ){

			double dis2GDSZone=myDistance2d( nextPoint, BoundStatic::zone );

			if( dis2GDSZone > 1000*BoundStatic::zone.width ){
				cout<<"point is too far away from the GDS::zone! END"<<endl;
				return;
			}
			dis2GDSZone+=0.5*GridOctree2d::gridCellSize;
			cube.setCoordinate( nextPoint.x1- dis2GDSZone , nextPoint.x1+dis2GDSZone, nextPoint.z1- dis2GDSZone, nextPoint.z1+dis2GDSZone   );

			generatePointOnCubeSurface( nextPoint, point2d, cube , normalDirOfCube, GreenVT, rand_01, rand_normalDir  );

			continue;
		}

		for( SubConductorList2d::iterator neighborCondIt= cube.neighborCondList.begin(); neighborCondIt != cube.neighborCondList.end(); neighborCondIt++   ){
			//next Point is on the surface of a subCond
			
			if(  fabs(  myDistance2d( nextPoint,  *neighborCondIt ) )< 1E-7  ){
				int slaveCondID = neighborCondIt->net_name;
				int slaveCondIndex = getCondIndexByID( condList, slaveCondID );
				capMatrix[slaveCondIndex].update( capacitance );
				return;
			}
		}
		generateMaxTransitionCube(  cube,  rootCell, leafCell ,nextPoint, gridOctree, normalDirOfCube );
		generatePointOnCubeSurface( nextPoint, point2d, cube, normalDirOfCube, GreenVT ,rand_01, rand_normalDir );	
		numOfStep++;
	}

}


//@author:Dragon
//before selecting point from cube, we should use srand( unsigned( time(NULL))) to initiate the seed (just once!)
void generatePointOnCubeSurface(FPoint2d &pointOnCube, FPoint2d &point2d,  Face &cube  ,int &normalDirOfCube , vector<double> &GreenVT ,//this green
	boost::variate_generator<boost::mt19937&,boost::uniform_01<> > &rand_01, 
    boost::variate_generator<boost::mt19937&,boost::uniform_int<> > &rand_normalDir ){

	// gettimeofday(&t_p_start, NULL );

	int rand_surface=rand_normalDir();
	normalDirOfCube =rand_surface;
	switch(rand_surface){
		case nX: {

			selectPointFromPlane(point2d, cube.z1, cube.z2,  cube.size(), greenFuncProgression_NX,greenFuncProgression_NY ,GreenVT ,rand_01 );
			pointOnCube.setCoordinate(  cube.x1,  point2d.x1  );
			break ;
		}
		case pX: {

			selectPointFromPlane(point2d, cube.z1, cube.z2,  cube.size(), greenFuncProgression_NX,greenFuncProgression_NY ,GreenVT , rand_01 );
			pointOnCube.setCoordinate( cube.x2,  point2d.x1 );
			break ;
		}
		case pZ: {

			selectPointFromPlane( point2d, cube.x1, cube.x2, cube.size(), greenFuncProgression_NX, greenFuncProgression_NY ,GreenVT, rand_01 );
			pointOnCube.setCoordinate(  point2d.x1,  cube.z2  );
			break ;
		}
		case nZ:{

			selectPointFromPlane( point2d, cube.x1, cube.x2, cube.size(), greenFuncProgression_NX, greenFuncProgression_NY ,GreenVT, rand_01 );
			pointOnCube.setCoordinate(  point2d.x1,  cube.z1  );
			break ;
		}
		default:{
			cout<<"Err: wrong hop direction!"<<endl;
			exit(1);
		}
	}
}

void generatePointOnCubeSurface( FPoint2d &fp,  FPoint2d &point2d,  Face &cube , vector<double> &GreenVT ,
	boost::variate_generator<boost::mt19937&,boost::uniform_01<> > &rand_01, 
    boost::variate_generator<boost::mt19937&,boost::uniform_int<> > &rand_normalDir  ){

	int temp;
	generatePointOnCubeSurface( fp, point2d ,cube, temp, GreenVT, rand_01, rand_normalDir );

	return ;
}


//use acception-rejection sampling method to create the desired distribution.
void selectPointFromPlane( FPoint2d &fp,   double x1, double x2, double size, int nX, int nY , vector<double>  &GreenVT, boost::variate_generator<boost::mt19937&,boost::uniform_01<> > &rand_01  ){
	double x;
	double vertical_value=0;
	double greenFunc_value=0;
	int nx;
	int GreenVTSize=GreenVT.size();

	while( true ){
		//vertical_value ranges from 0 to 1. That's enough because the max value of greenFunc_value is less than 0.5
		x=rand_01()*size;
		vertical_value=rand_01()*0.5;  
		// vertical_value=rand_01() * (sqrt(3)/10);  
		nx= int( x/(size)/(1.0/GreenVTSize));
		if(nx<0||nx> GreenVTSize ){
			cout<<"Wrong green surface locate!!"<<endl;
			exit(1);
		}
  	 	nx=	(nx==GreenVTSize? nx-1: nx);
		greenFunc_value=GreenVT[nx];
		if(vertical_value <= greenFunc_value){
			fp.setCoordinate( x1+ x ,double(0.0)); 
			return;
		}

	}
}

GaussianSurface2d generateBGS( Subnet subCond, TestList condList , double scale_factor )
{
	std::vector<double> tempExtensionDis(4);
	int dirOfRect=0;
	int disMax;
	vector<Subnet>::iterator viter;


	for(  int dirOfRect  =0; dirOfRect <4; dirOfRect++  ){
		disMax=distanceToGdsBorder2d( subCond, dirOfRect  );
		for(condList.test_iter = condList.testlist.begin() ;condList.test_iter != condList.testlist.end() ;condList.test_iter++)
		{

			for(viter = condList.test_iter->net->subnets.begin(); viter != condList.test_iter->net->subnets.end();viter++)
			{

				if( subCond.x0 == viter->x0 && subCond.x1 == viter->x1 && subCond.z0 == viter->z0 && subCond.z1 == viter->z1){
					continue;
				}
				int distanceAB= myDistance2d( subCond, *viter);

				if( condList.test_iter->net_name == subCond.net_name && distanceAB==0  ){
					continue;
				}

				if(  myDistance2d( subCond, *viter, dirOfRect) == distanceAB  ){
					disMax=min(  disMax,  distanceAB  );
				}
			}

		}
		//disMax can be 0 when two subConds are adjecent. This is permitted.
		tempExtensionDis[dirOfRect]=disMax/2.0;
		if( fabs(tempExtensionDis[dirOfRect])< 1E-5  ){
			cout<<"The extensionSize of BGS cannot be 0!!"<<endl;
			exit(1);
		}
	}

// 	// cout<<"------------------"<<extensionDis[NX]<<","<<extensionDis[PX]<<","<<extensionDis[NY]<<","<<extensionDis[PY]<<","<<extensionDis[NZ]<<","<<extensionDis[PZ]<<endl;
// //@author:Dragon
// //if two subConductors are adjacent, then on this side the extension distance of Gaussian Surface is 0, this side does not need to be considered because it will not be sampled
	double extensionDistanceLimit=0;
	for( dirOfRect=0; dirOfRect <4; dirOfRect++  ){
		if( extensionDistanceLimit==0  ){
			extensionDistanceLimit=tempExtensionDis[dirOfRect];
		}else{
			if( tempExtensionDis[dirOfRect]  < extensionDistanceLimit ){
				extensionDistanceLimit=tempExtensionDis[dirOfRect];
			}
		}
	}

	extensionDistanceLimit= extensionDistanceLimit*scale_factor;

	for(dirOfRect=0; dirOfRect<4; dirOfRect++ ){
		tempExtensionDis[dirOfRect]=min( extensionDistanceLimit,  tempExtensionDis[dirOfRect]  );
	}

	return  GaussianSurface2d( subCond.x0-tempExtensionDis[nX], subCond.x1+tempExtensionDis[pX],
					subCond.z0- tempExtensionDis[nZ], subCond.z1+ tempExtensionDis[pZ]  , tempExtensionDis  );

}

void generateGaussianSurfaceOfConductor( TestList &condList, double scale_factor  )
{
	vector<Subnet>::iterator viter;

    for(condList.test_iter = condList.testlist.begin() ;condList.test_iter != condList.testlist.end() ;condList.test_iter++)
    {
		for(viter = condList.test_iter->net->subnets.begin(); viter != condList.test_iter->net->subnets.end();viter++)
		{

			GaussianSurface2d blockGaussianSurface=generateBGS( *viter, condList,scale_factor );
			condList.test_iter->conductor2d.gaussianSurfaceList.push_back( blockGaussianSurface  );
			double minExtensionDisOfBGS=min(  blockGaussianSurface.extensionDis[pX]  , 
											blockGaussianSurface.extensionDis[pZ]  ,  
											blockGaussianSurface.extensionDis[nX]  , 
											blockGaussianSurface.extensionDis[nZ]  );
			
			if(condList.test_iter->conductor2d.gaussianSurfaceList.minExtensionDis==0){
				condList.test_iter->conductor2d.gaussianSurfaceList.minExtensionDis= minExtensionDisOfBGS;
			}else{
				if(   condList.test_iter->conductor2d.gaussianSurfaceList.minExtensionDis > minExtensionDisOfBGS   ){
					condList.test_iter->conductor2d.gaussianSurfaceList.minExtensionDis = minExtensionDisOfBGS;
				}
			}
			condList.test_iter->conductor2d.gaussianSurfaceList.VGSArea+=2.0/3.0*blockGaussianSurface.area();
			condList.test_iter->conductor2d.gaussianSurfaceList.SumOfBGSArea+=blockGaussianSurface.area();
		}
	}

	//检查与一个子高斯面相邻的所有子高斯面
	for(condList.test_iter = condList.testlist.begin() ;condList.test_iter != condList.testlist.end() ;condList.test_iter++){	
		neighborGaussianSurfaceCheck(  condList.test_iter->conductor2d.gaussianSurfaceList);
	}
}

int distanceToGdsBorder2d( Subnet rectA ,int normalDirect ){
	switch(normalDirect){
		case pX:{
			return  BoundStatic::zone.x1-rectA.x1;
		}
		case nX:{
			return	rectA.x0-BoundStatic::zone.x0;
		}
		case pZ:{
			return BoundStatic::zone.z1-rectA.z1;
		}
		case nZ:{
			return rectA.z0-BoundStatic::zone.z0;
		}
		default:{
			cout<<"Error: wrong normalDirect!"<<endl;
			exit(1);
		}
	}
}

int myDistance2d( Subnet rectA,Subnet rectB, int dirOfRectA ){
	switch(dirOfRectA){
		case pX: return rectB.x0-rectA.x1;
		case pZ: return rectB.z0-rectA.z1;
		case nX: return rectA.x0-rectB.x1;
		case nZ: return rectA.z0-rectB.z1;
		default:{
			cout<<"Error: wrong direction of Rect!"<<endl;
			exit(1);
		}
	}
}

void neighborGaussianSurfaceCheck( GaussianSurfaceList2d &BGSList ){

	GaussianSurfaceList2d tempList=BGSList;
	for( GaussianSurfaceList2d::iterator BGSIt = BGSList.begin(); BGSIt != BGSList.end(); BGSIt ++  ){
		for( GaussianSurfaceList2d::iterator otherBGSIt=tempList.begin(); otherBGSIt != tempList.end() ; otherBGSIt++ ){
			if(  *BGSIt == *otherBGSIt   ){
				continue;
			}
			if( BGSIt->isNeighborWith(  *otherBGSIt )  ){
				 BGSIt->neighborList.push_back( *otherBGSIt );
			}
		}
	}
}

void sampleOnVGSS( FPoint2d &pointOnVGS, GaussianSurfaceList2d &BGSList ,double &PDFOnVGS  ,int &normalDirOfVGS , int &Nt, int &Ng , double &PDFIntegralElement, boost::variate_generator<boost::mt19937&,boost::uniform_01<> > &rand_01 )
{
	double randArea=rand_01()*double(BGSList.SumOfBGSArea);//rand is a problem
	GaussianSurfaceList2d::iterator BGSIt;
	while(true){
		for(BGSIt=BGSList.begin(); BGSIt != BGSList.end(); BGSIt ++ ){
			randArea-=BGSIt->area();
			if( randArea<=0 ){
				break;
			}
		}
		if( BGSIt != BGSList.end() ){
			break;
		}
	}
	//randomly select a point from the BGS uniformly.
	//enum NormalDirection{  PX, PY, PZ, NX, NY, NZ };
	randArea =rand_01()*BGSIt->area();
	int normalDir;
	while(true){
		for(normalDir=0; normalDir <4 ; normalDir++  ){
			randArea-=BGSIt->surfaceArea[normalDir];
			if(randArea<=0 ){
				break;
			}

		}
		if(normalDir !=4 ){
			break;
		}
	}
	double xLen=BGSIt->x2-BGSIt->x1;
	double zLen=BGSIt->z2-BGSIt->z1;
	double dX,dY,dZ;
	dX=rand_01()*xLen;
	dZ=rand_01()*zLen;
	switch(normalDir){
		case pX: pointOnVGS.setCoordinate( BGSIt->x2 , BGSIt->z1+dZ ); break;
		case pZ: pointOnVGS.setCoordinate( BGSIt->x1+dX ,  BGSIt->z2  ); break;
		case nX: pointOnVGS.setCoordinate( BGSIt->x1 ,  BGSIt->z1+dZ ); break;
		case nZ: pointOnVGS.setCoordinate( BGSIt->x1+dX ,  BGSIt->z1  ); break;
		default:{cout<<"Err: wrong normal direction"<<endl; exit(1);}
	}

	Nt++;

	int nc=1;
	int normalDirOfPointOnBGS=BGSIt->getNormalDirOfPoint(  pointOnVGS );
	int normalDirOfPointOnOtherBGS;

	for(  std::list<GaussianSurface2d>::iterator neighborBGSIt= BGSIt->neighborList.begin(); neighborBGSIt !=BGSIt->neighborList.end();  neighborBGSIt++  ){
	
		if(   neighborBGSIt->isPointEnclosedBySurface( pointOnVGS )  ){
			sampleOnVGSS( pointOnVGS, BGSList, PDFOnVGS  ,  normalDirOfVGS , Nt, Ng, PDFIntegralElement, rand_01 );
			return;
		}

		normalDirOfPointOnOtherBGS=neighborBGSIt->getNormalDirOfPoint(pointOnVGS  );
		if( normalDirOfPointOnOtherBGS ){
			if( ( normalDirOfPointOnBGS+normalDirOfPointOnOtherBGS )==0 ){
				sampleOnVGSS(pointOnVGS, BGSList ,PDFOnVGS, normalDirOfVGS  , Nt, Ng, PDFIntegralElement, rand_01 );
				return;
			}
			nc++;
		}
	}
	double r1= rand_01();
	if( r1 > 1.0/nc ){
		sampleOnVGSS(pointOnVGS, BGSList ,PDFOnVGS,normalDirOfVGS   , Nt, Ng, PDFIntegralElement, rand_01  );
		return;
	}
	//@author:Dragon
	// up to here we sample a point from the virtual gaussian surface uniformly. so at this point Ng should add 1
	Ng++;
	normalDirOfVGS=normalDir;

}

double surfaceGreenFunction(  Face &cube, FPoint2d &p , int normalDirOfCube, vector<double> &GreenVT ){
	double x;
	switch( normalDirOfCube ){
		case nX:
		case pX: x= p.z1-cube.z1; break;
		case pZ:
		case nZ: x= p.x1-cube.x1; break;
		default:{
			cout<<"Err: wrong cube dir!"<<endl;
			exit(1);
		}
	}
	double cubeSize=cube.size();
    int GreenVTSize=GreenVT.size();
	x=x/cubeSize;
	int nx= int( x /(1.0/GreenVTSize));
	if(nx<0||nx> GreenVTSize ){
		cout<<"Wrong green surface locate!!"<<endl;
		exit(1);
	}
    nx=	(nx==GreenVTSize? nx-1: nx);

	return GreenVT[nx]/double(cubeSize);////? in the calculate we divide L not L2
}

double GEx( double x ,double L, vector<double> &GExVT )
{
	x=x/L;
	int GExVTSize=GExVT.size();
	int nx= int( x /(1.0/GExVTSize));
	if(nx<0||nx> GExVTSize  ){
		cout<<"Wrong GEx surface locate!!"<<endl;
		exit(1);
	}
    nx=	(nx==GExVTSize? nx-1: nx);
    return GExVT[nx]/pow(L, 2.0);//?
}

double GEz(double x,double L ,vector<double>  &GEzVT )
{

	x=x/L;
	int GEzVTSize=GEzVT.size();
	int nx= int( x /(1.0/GEzVTSize));
	if(nx<0||nx> GEzVTSize){
		cout<<"Wrong GEz surface locate!!"<<endl;
		exit(1);
	}
    nx=	(nx==GEzVTSize? nx-1: nx);
    return GEzVT[nx]/pow(L, 2.0);//?

}

double greenFunctionForElectricField(  Face &cube,  FPoint2d &p,  int normalDirOfVGS, int normalDirOfCube ,vector<double> &GExVT, vector<double>  &GEzVT  ){	
	double GE=0;
	double L=cube.size();
	int sign;
	double x;
	double y;
	switch(normalDirOfCube ){
		case nX:
		case pX: x = p.z1-cube.z1;break;
		case nZ:
		case pZ: x = p.x1-cube.x1;break;
		default:{
			cout<<"wrong normalDirOfCube!"<<endl;
			exit(1);
		}
	}

	if( normalDirOfVGS==pZ || normalDirOfVGS==nZ  ){
		if( normalDirOfCube==pZ ||normalDirOfCube==nZ ){
			sign=normalDirOfCube==pZ?1:-1;
			GE=sign* GEz( x, L , GEzVT );
		}
		else
		{
			GE=GEx(x,L ,GExVT );
		}

	}else if( normalDirOfVGS==nX || normalDirOfVGS==pX ){
		if(normalDirOfCube==nX||normalDirOfCube==pX ){
			sign= normalDirOfCube==pX?1:-1;
			GE=sign* GEz( x, L, GEzVT );
		}
		else 
		{
			GE=GEx( x, L, GExVT );
		}

	}else{
		cout<<"wrong normal direction!!"<<endl;
		exit(1);
	}

	return GE;
}


Cell2d* locateOctreeRootCell(  FPoint2d &p, GridOctree2d &gridOctree  ){
	int gridCellSize=GridOctree2d::gridCellSize;
	int xn= int (p.x1 - BoundStatic::zone.x0 )/gridCellSize;
	int zn= int (p.z1 - BoundStatic::zone.z0 )/gridCellSize;

	Cell2d *rootCell;
	rootCell = &(gridOctree[xn][zn]);
	if(p.x1-rootCell->x1<-1E-7||p.x1-rootCell->x2>1E-7||p.z1-rootCell->z1<-1E-7||p.z1-rootCell->z2>1E-7 ){
		cout<<"p is not in root Cell2d!!!"<<endl;
		cout<<"P:"<<p.x1<<","<<p.z1<<endl;
		cout<<"CELL:"<<rootCell->x1<<","<<rootCell->x2<<","<<rootCell->z1<<","<<rootCell->z2<<endl;
		exit(1);
	}
	return rootCell;
}


Cell2d* locateCellFromRoot( FPoint2d &p, Cell2d *rootCell , Face &cube ){
	
	if(rootCell->isFilledWithCond ){
		cout<<"FRW error in  locateCellFromRoot, FRWcapext.cpp: FRW walks into the inside of a conductor! "<<endl;
		exit(1);
	}

	if( !rootCell->hasChildCell ){
		return rootCell;
	}

	for( SubCellList2d::iterator subCellIt = rootCell->subCellList.begin(); subCellIt != rootCell->subCellList.end(); subCellIt++     ){
		if(   p.x1 - subCellIt->x1>-1E-7 && p.x1 -subCellIt->x2 <1E-7 
			&& p.z1 - subCellIt->z1 > -1E-7 && p.z1 - subCellIt->z2 <1E-7 )
			{
			//@author:Dragon
			//这里放宽了定位的条件，目的是确保一定可以定位到cell里面

			return locateCellFromRoot(  p, &(*subCellIt), cube );
			
		}

	}
	cout<<"locate Cell accoring to Point failed !!!"<<endl;
	exit(1);
}





void generateMaxTransitionCube(Face &maxCube, FPoint2d &p, Cell2d* &cell ){
//@author:Dragon
// find the minimum distance from point to the candidate SubConductors, then compare the value with distance limit the lesser is the size of transition cube
//the initial value of cubeRadius is the distance from point to the border of inflated Cell.	

	double x_Dis=max( cell->x1-cell->distanceLimit-p.x1, p.x1- ( cell->x2 +cell->distanceLimit  )    );
	double z_Dis=max( cell->z1-cell->distanceLimit-p.z1, p.z1- ( cell->z2 +cell->distanceLimit )    );

	double cubeRadius= fabs( max( x_Dis,z_Dis));

	// cout<<"IIIIOOOOOOIIIIII"<<endl;
// the traditional way:
	

	for( SubConductorList2d::iterator candidateIt= cell->candidateList2d.begin(); candidateIt != cell->candidateList2d.end(); candidateIt++  ){
		//the improved way!!
		if(  myDistance2d ( *candidateIt, *cell )  - cubeRadius >1E-7   ){    //should not be >=
			break;
		}
		double tempDis=myDistance2d( p, *candidateIt) ;
		if(  tempDis -cubeRadius <= -1E-7   ){
			cubeRadius=tempDis ;
			maxCube.neighborCondList.clear();
			maxCube.neighborCondList.push_back( *candidateIt  );
		}else if( fabs(tempDis - cubeRadius)<1E-7  ){
			maxCube.neighborCondList.push_back( *candidateIt );
		}

	}

	maxCube.setCoordinate( p.x1-cubeRadius , p.x1 +cubeRadius, p.z1-cubeRadius, p.z1+cubeRadius );
}


void generateMaxTransitionCube( Face &maxCube, Cell2d* &rootCell, Cell2d* &leafCell,  FPoint2d &p, GridOctree2d &gridOctree, int &normalDirOfCube ){

	totalStepNum++;

	if( leafCell!=NULL && p.x1 - leafCell->x1> -1E-7 && p.x1 - leafCell->x2 < 1E-7 &&
		p.z1 - leafCell->z1> -1E-7 && p.z1 - leafCell->z2 < 1E-7   ){
		hitTime++;
	}	
	else{ 
		rootCell= locateOctreeRootCell( p, gridOctree  );

		int nx=(int) ((p.x1-rootCell->x1)/Cell2d::subGridSize);
		int nz=(int) ((p.z1-rootCell->z1)/Cell2d::subGridSize);
		int subGridNum=rootCell->subGridCellVec.size();///?

		nx= (nx==subGridNum? nx--: nx);
		nz= (nz==subGridNum? nz--: nz);
		leafCell= rootCell->subGridCellVec[nx][nz];

		if(leafCell->hasChildCell){
			leafCell= locateCellFromRoot( p, leafCell, maxCube );
		}

	}
	generateMaxTransitionCube(  maxCube ,p,  leafCell );	
	
}

/*
Cij  the capacitance from conductor i to conductor j

C11  C12 ... C1n
C21  C22 ... C2n
...
Cn1  Cn2 ... Cnn
*/
////////////////////////////////////////////////////////////////////////
CapacitanceMatrix generateCapMatrix(  TestList &condList ){
	int numCond= condList.testlist.size();
	if(numCond==0){
		cout<<"Err: the number of conductor is 0!"<<endl;
		exit(1);
	}
	std::vector<Capacitance> capMatrix(numCond);

	for(int i=0; i<capMatrix.size(); i++ ){
			capMatrix[i].masterCondID=i;
	}
	return capMatrix;
}

int getCondIndexByID( TestList &condList , int net_name){
	int index=-1;
	for(condList.test_iter = condList.testlist.begin() ;condList.test_iter != condList.testlist.end() ;condList.test_iter++)
		{
		index++;
		if(  condList.test_iter->net_name==net_name  ){
			return index;
		}
	}
	return -1;
}

vector<string> split( string str, string pattern){
	vector<string> ret;
	if(pattern.empty()){
		return ret;
	}
	size_t start=0,index=str.find_first_of(pattern,0);
	while(index!=string::npos)
	{
		if(start!=index){
			ret.push_back(str.substr(start,index-start));
		}
		start=index+1;
		index=str.find_first_of(pattern,start);
	}
	if(!str.substr(start).empty()){
		ret.push_back(str.substr(start));
	}
	return ret;
}

void loadDataTable(const char *dataFile, vector<vector<double> > &dataTable , int &tableSize ){
	ifstream fin(dataFile);
	string line;
	while(getline(fin, line)){
		std::vector<string> strVec=split(line,"\t");
		std::vector<double> numVec;
		for(int i=0; i<strVec.size(); i++){
			double temp;
			if(sscanf( strVec[i].c_str(),"%lf", &temp ) !=1){
				cout<<"Data Table format error!"<<endl;
				exit(1);
			}
			numVec.push_back(temp);
		}
		dataTable.push_back(numVec);
	}
	fin.close();
	tableSize=dataTable.size();
}

void loadDataTable(const char *dataFile, vector<double>  &dataTable , int &tableSize ){
	ifstream fin(dataFile);
	string line;
	while(getline(fin, line)){
		std::vector<string> strVec=split(line,"\t");
		for(int i=0; i<strVec.size(); i++){
			double temp;
			if(sscanf( strVec[i].c_str(),"%lf", &temp ) !=1){
				cout<<"Data Table format error!"<<endl;
				exit(1);
			}
			dataTable.push_back(temp);
		}
	}
	fin.close();
	tableSize=dataTable.size();
}


void copy(ofstream &fout, string &filename ){
	ifstream fin(filename.c_str());
	string str;
	while(getline(fin,str)){
		fout<<str<<"\n"<<flush;
	}
	fin.close();
}

void convertcondlist(map<int, Net*> condList)
{
	map<int, Net*>::iterator subCondIt;
	Conductor2d gaussianSurfaceList2d;
	OneConduct oneconduct;
	TestList testlist;//Testlist class will replace the conductlist
	for(subCondIt = condList.begin() ;subCondIt != condList.end() ;subCondIt++)
    {
		oneconduct.set(subCondIt->second,gaussianSurfaceList2d,subCondIt->first);
		testlist.set(oneconduct);
	}
	FRWControl::teatlist = testlist;

};