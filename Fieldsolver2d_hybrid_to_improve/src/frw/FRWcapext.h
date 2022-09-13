#include <iostream>
#include <string>
#include "spacemanagement.cpp"
#include <sys/time.h>
#include <unistd.h>
#include <boost/random.hpp>
#include <pthread.h>
#include <cstring>   // string
#include <fstream>   // ifstream
#include <regex>     // for split()
#include <vector>    // vector
#include <algorithm> // sort
#include <cmath>     // fabs
#include <map>       // map
#include <iomanip>   // setprecision
#include <sys/resource.h> // getrusage
#include <list>      // list
#include <sstream>   // istringstream, for split()

using namespace std;

const int greenFuncProgression_NX=3;////??5
const int greenFuncProgression_NY=3;////??5


// double testTime=0;
double generateMaxCubeTime=0;
double generatePointTime=0;
double locateCellTime=0;
struct timeval t_cube_start,t_cube_end, t_p_start, t_p_end, t_locatecell_start, t_locatecell_end  ;
struct timeval tv_begin,tv_end;


long hitTime=0;
long totalStepNum=0;

pthread_mutex_t getRandSeed_mutex;
pthread_mutex_t updateFRWData_mutex;
pthread_mutex_t processorIndex_mutex;  //在设定线程亲和力的时候对线程加锁

class Capacitance{
public:
	int masterCondID;
	double sumOfCap;
	double sumOfCapSquared;
	double capacitance;
	double estimateErr;
	double std_var;

	Capacitance(){
		masterCondID=-1;
		sumOfCap=0;
		sumOfCapSquared=0;
		capacitance=0;
		estimateErr=10000;
		std_var=10000;
	}
	void update( double cap ){
		sumOfCap+=cap;
		sumOfCapSquared+=pow(cap,2.0);
	}

	void update(Capacitance &capEntity){
		sumOfCap+= capEntity.sumOfCap;
		sumOfCapSquared+= capEntity.sumOfCapSquared;
	}

	void reset(){
		sumOfCap=0;
		sumOfCapSquared=0;
		capacitance=0;
		estimateErr=10000;
		std_var=10000;
	}

};

typedef  std::vector<Capacitance>  CapacitanceMatrix;
// typedef std::vector<Capacitance>  CapacitanceMatrix2d;


CapacitanceMatrix generateCapMatrix( TestList &condList );

// int getCondIndexByID( ConductorList &condList, int  ID );


class FRWControl{
public:
	static bool isRandSeedSet;
	static double dielectricConstant;
	// static ConductorList conductorList;
	static ConductorList2d  conductorList2dMove;
	static TestList teatlist;
	static map<int, Net*>  conductorList2d;
	// static GridOctree gridOctree;
	static GridOctree2d gridOctree2d;
	static CapacitanceMatrix capMatrix;
	// static CapacitanceMatrix2d capMatrix2d;
	static vector<double>  GreenVT; 
	static vector<double>  GExVT;
	static vector<double>  GEzVT;
	static int GreenVTSize;
	static int GExVTSize;
	static int GEyVTSize;
	static int GEzVTSize;

	static int totalProcessorNum;
	static int currentProcessorIndex;

	static int currentThreadIndex;

	static vector<long> totalNumOfWalkOfCond;
	static vector<long> currentNumOfWalkOfCond;
	static vector<long> currentProgressOfCond;
	static bool updateProgress(int tempWalkNum, int masterCondIndex ,double temp_VGSArea );

	static vector<long> randSeedVec;
	static double maxPermittedCapErr;


};
bool FRWControl::isRandSeedSet=false;
double FRWControl::dielectricConstant=3.9*8.85/1000;        //unit: pF/mz
double FRWControl::maxPermittedCapErr=2;
vector<long> FRWControl::totalNumOfWalkOfCond;
vector<long> FRWControl::currentNumOfWalkOfCond;
vector<long> FRWControl::currentProgressOfCond;

int FRWControl::GreenVTSize=0;
int FRWControl::GExVTSize=0;
int FRWControl::GEzVTSize=0;
int FRWControl::totalProcessorNum=0;
int FRWControl::currentProcessorIndex=0;
int FRWControl::currentThreadIndex=0;

// ConductorList FRWControl::conductorList; //= ConductorList();
ConductorList2d FRWControl::conductorList2dMove;
map<int, Net*>  FRWControl::conductorList2d;
TestList FRWControl::teatlist;
GridOctree2d FRWControl::gridOctree2d;//= GridOctree2d();
CapacitanceMatrix FRWControl::capMatrix;//= CapacitanceMatrix();
vector<double>   FRWControl::GreenVT;
vector<double>   FRWControl::GExVT;
vector<double>   FRWControl::GEzVT;
vector<long> FRWControl::randSeedVec;


vector<string> split( string str, string pattern);
void loadDataTable(const char *dataFile, vector<vector<double> > &dataTable, int &tableSize  );
void loadDataTable(const char *dataFile, vector<double>  &dataTable , int &tableSize );
void copy(ofstream &fout, string &filename );

class Face{
public:
	double x1,x2,z1,z2;

	SubConductorList2d neighborCondList;

	Face();
	Face( FPoint2d &minPoint, double size );
	Face( double x_1, double z_1, double size);
	Face(double x_1, double x_2,double z_1, double z_2  );
	Face(double x_1, double x_2 );

	double size();
	void setCoordinate(  double x_1, double x_2,  double z_1, double z_2  ){
		x1=x_1;
		x2=x_2;
		z1=z_1;
		z2=z_2;
	}

};

void generatePointOnCubeSurface( FPoint2d &pointOnCube, FPoint2d &point2d,  Face &cube , int &normalDirOfCube , vector<double>  &GreenVT,
	boost::variate_generator<boost::mt19937&,boost::uniform_01<> > &rand_01, 
    boost::variate_generator<boost::mt19937&,boost::uniform_int<> > &rand_normalDir );
void generatePointOnCubeSurface( FPoint2d &pointOnCube, FPoint2d &point2d, Face &cube  , vector<double>  &GreenVT,
	boost::variate_generator<boost::mt19937&,boost::uniform_01<> > &rand_01, 
    boost::variate_generator<boost::mt19937&,boost::uniform_int<> > &rand_normalDir );

void selectPointFromPlane( FPoint2d &fp,   double x1, double x2,  double size, int nX, int nY , vector<double> &GreenVT ,
		boost::variate_generator<boost::mt19937&,boost::uniform_01<> > &rand_01);

class GaussianSurface2d{
public:
	double x1,x2,z1,z2;
	std::list<GaussianSurface2d> neighborList;

	GaussianSurface2d();
	GaussianSurface2d( FPoint2d &minPoint, FPoint2d &maxPoint , std::vector<double> extensionDistance);
	GaussianSurface2d(double x_1, double x_2, double z_1, double z_2, std::vector<double> extensionDistance  );

	double area();
	bool isPointOnSurface(FPoint2d &point);
	bool isPointEnclosedBySurface( FPoint2d &point );
	int getNormalDirOfPoint(FPoint2d  &point);
	std::vector<double> extensionDis;
	std::vector<double> surfaceArea;
	void setSurfaceArea();
	bool operator ==( GaussianSurface2d &rhs  );
	bool isNeighborWith( GaussianSurface2d &gaussian );

};



enum NormalDirection2d{ pX, pZ, nX, nZ };
// void generateGaussianSurfaceOfConductor( ConductorList &condList, double scale_factor  );

GaussianSurface2d generateBGS(Subnet subCond, TestList condList   , double scale_factor   );
int distanceToGdsBorder2d( Subnet rectA ,int normalDirect );
// //@author:Dragon
// //myDistance take the relative direction of rectA and rectB into consideration.
int myDistance2d( Subnet rectA,Subnet rectB, int dirOfRectA);




void neighborGaussianSurfaceCheck( GaussianSurfaceList2d &BGSList );

void sampleOnVGSS( FPoint2d &pointOnVGS, GaussianSurfaceList2d &BGSList ,double &PDFOnVGS  ,int &normalDirOfVGS , int &Nt, int &Ng , double &PDFIntegralElement, boost::variate_generator<boost::mt19937&,boost::uniform_01<> > &rand_01);

double surfaceGreenFunction(  Face &cube,  FPoint2d &p, int normalDirOfCube, vector<double>  &GreenVT );

double greenFunctionForElectricField(  Face &cube,  FPoint2d &p,  int normalDirOfVGS, int normalDirOfCube ,vector<double> &GExVT, vector<double>  &GEzVT  );

Cell2d *locateOctreeRootCell( FPoint2d &p, GridOctree2d &gridOctree );

Cell2d *locateCellFromRoot(  FPoint2d &p, Cell2d *rootCell ,Face &cube );


void FloatingRandomWalk(  FPoint2d &pointOnVGS, Face &cube, Cell2d* &rootCell, Cell2d* &leafCell, FPoint2d &nextPoint,  FPoint2d &point2d, OneConduct &cond,  TestList &condList,  GridOctree2d &gridOctree,  CapacitanceMatrix &capMatrix , 
		std::vector<double> &GreenVT, vector<double> &GExVT,  vector<double>  &GEzVT,
		boost::variate_generator<boost::mt19937&,boost::uniform_01<> > &rand_01, 
    boost::variate_generator<boost::mt19937&,boost::uniform_int<> > &rand_normalDir , int &Nt, int &Ng );



void generateMaxTransitionCube(Face &maxCube, FPoint2d &p, Cell2d* &cell );
void generateMaxTransitionCube( Face &maxCube, Cell2d* &rootCell, Cell2d* &leafCell, FPoint2d &p, GridOctree2d &gridOctree,int &normalDirOfCube );


CapacitanceMatrix generateCapMatrix( map<int, Net*>  condList );

int getCondIndexByID( TestList &condList , int net_name);

