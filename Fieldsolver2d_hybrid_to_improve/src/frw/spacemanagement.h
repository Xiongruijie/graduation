#include <iostream>
#include "geo2condlist/geoloader.cpp"

using namespace std;

//这个表示从根节点（即第一级的3D数组中的一个单元）开始看，八叉树的深度
int depthOfOctreeFromRoot=4;
//这个表示从第二级3D数组的一个单元开始看，八叉树的深度
//depthOfOctreeFromSubGrid=depthOfOctreeFromRoot时表示没有第二级grid,当depthOfOctreeFromSubGrid等于0时，意味着在定位Cell的时候使用的
//是两级3D数组，并没有应用八叉树
int depthOfOctreeFromSubGrid=0;


int max( int a, int b );
double max( double a, double b  );

int min(int a, int b);
int min( int a1, int a2, int a3, int a4, int a5  , int a6  );
double min( double a1, double a2 );
double min( double a1, double a2, double a3, double a4, double a5, double a6  );


class FPoint2d{
public:
	double x1,z1;
	FPoint2d();
	FPoint2d(double x_1, double z_1 );
	void setCoordinate( double x_1, double z_1 ){
		x1=x_1;
		z1=z_1;
	}
	
};

//@author:Dragon
//为了精确计算，把Cell的坐标改成double类型
class Cell2d
{

public:
	static  int maxCandidateListLen;
	static  double minCellSize;                 //2 times of the min width (nm)
	static  double subGridSize;                 //这个变量是用来规定每一个八叉树空间的栅格的尺寸

	double x1, x2,z1, z2;

	SubConductorList2d candidateList2d;

	SubConductorList2d candidateListOfRootCell2d;

	double extensionSize;
	double distanceLimit;
	bool isFilledWithCond;  //to express if the Cell is totally filled with SubConductor. Note: subConductor, not the whole Condcutor!
	bool isFilledChecked;	//to express if the FilledWithCond is checked.
	list<Cell2d> subCellList;     //using list is more flexible than vector


	list<Cell2d*> candidateCellList;
	FPoint2d centerPoint;
	bool hasChildCell;

	vector<vector<Cell2d*> >  subGridCellVec;

	Cell2d();
	Cell2d( double x_1, double x_2, double z_1, double z_2  , double extension);
	Cell2d( double x_1, double x_2, double z_1, double z_2  , double extension, double disLimit);

	double size();
	void divideIntoSubCell( double extension );
	Cell2d getInflatedCell( double extension  );  // not used for now
	bool isInsectWith(Subnet subCond );
	bool isIncludedIn(Subnet subCond );

};

//注意在本例子中一个Cell肯定是标准的立方体

typedef std::list<Cell2d> SubCellList2d;

class GridOctree2d:public std::vector< std::vector<Cell2d> > {
public:
	static int gridCellSize;
};


int GridOctree2d::gridCellSize=1*1000 ;//320
int Cell2d::maxCandidateListLen=4;
//这个表示从根节点（即第一级的3D数组中的一个单元）开始看，八叉树的深度
int depthOfOctreeFromRoot2d=4;
//这个表示从第二级3D数组的一个单元开始看，八叉树的深度
//depthOfOctreeFromSubGrid=depthOfOctreeFromRoot时表示没有第二级grid,当depthOfOctreeFromSubGrid等于0时，意味着在定位Cell的时候使用的
//是两级3D数组，并没有应用八叉树
int depthOfOctreeFromSubGrid2d=0;

double Cell2d::minCellSize= GridOctree2d::gridCellSize /pow(2.0, double(depthOfOctreeFromRoot2d) ) ;
double Cell2d::subGridSize= pow(2.0, double(depthOfOctreeFromSubGrid2d))*Cell2d::minCellSize ;     

GridOctree2d partitionDomain(  double cellSize, int &nX, int &nY, Boundary GDS_zone);

void  insertToOctree( Subnet subCond, Cell2d &cellT);

void candidateCheck2d ( Subnet subCond,  Cell2d &cellT );

double myDistance2d(Subnet rectA, Cell2d &cellT );

int myDistance2d(Subnet rectA ,Subnet rectB );

double myDistance2d( FPoint2d  &p, Boundary &rect );

bool dominate2d(  Subnet rectA ,Subnet rectB, Cell2d &cellT );

void updateSubGridCellVec2d(Cell2d &rootCell,  Cell2d *cellP,  vector<vector<Cell2d*> > &subGridCellVec );
void updateSubGridCellVec2d(GridOctree2d &gridOctree);



