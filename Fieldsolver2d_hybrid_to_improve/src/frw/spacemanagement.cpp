#include <iostream>
#include "spacemanagement.h"

using namespace std;

FPoint2d::FPoint2d(){
	x1=z1=0;
}

FPoint2d::FPoint2d( double x_1, double z_1 ){
	x1=x_1;
	z1=z_1;
}

Cell2d::Cell2d(double x_1, double x_2, double z_1, double z_2 , double extension )
{
	x1=x_1;
	x2=x_2;
	z1=z_1;
	z2=z_2;
	distanceLimit= extensionSize = extension;
	isFilledWithCond=false;
	isFilledChecked=false;
	hasChildCell=false;
}

double Cell2d::size(){
	if( x2-x1<1E-5 ||z2-z1<1E-5 ){
		cout<<"The size of Cell is 0 !";
		exit(1);
	}
	return max(x2-x1, z2-z1) ;
}

bool Cell2d::isInsectWith(Subnet subCond){
	return x2-subCond.x0 >1E-7 && x1-subCond.x1<-1E-7&&z2-subCond.z0>1E-7 && z1-subCond.z1<-1E-7; 
}

bool Cell2d::isIncludedIn(Subnet subCond){
	return x2-subCond.x1<-1E-7&&x1-subCond.x0>1E-7&&
			 z2-subCond.z1<-1E-7&&z1-subCond.z0>1E-7 ;
}

void generateGridOctree( GridOctree2d &gridOctree2d,  int gridCellSize ,Boundary GDS_zone ,map<int, Net*> generateConductorList, TestList &condList){
	int nX=0,nZ=0;

	GDS_zone.x0 -=gridCellSize;//dext
	GDS_zone.x1 +=gridCellSize;
	GDS_zone.z0 -=gridCellSize;
	GDS_zone.z1 +=gridCellSize;
//static
	BoundStatic::zone.x0 = GDS_zone.x0;
	BoundStatic::zone.x1 = GDS_zone.x1;
	BoundStatic::zone.z0 = GDS_zone.z0;
	BoundStatic::zone.z1 = GDS_zone.z1;
	BoundStatic::zone.width = GDS_zone.x1 - GDS_zone.x0;
	BoundStatic::zone.height = GDS_zone.z0 - GDS_zone.z1;

	gridOctree2d=partitionDomain(double(gridCellSize) , nX, nZ ,GDS_zone);

	int ix1,ix2,iz1,iz2;
	int tempI;
	vector<Subnet>::iterator viter;

	for(condList.test_iter =condList.testlist.begin() ;condList.test_iter != condList.testlist.end() ;condList.test_iter++ )
    {

		for(viter = condList.test_iter->net->subnets.begin(); viter != condList.test_iter->net->subnets.end();viter++)
		{
			//扩大检查初始的candidateList的范围，尽可能产生更大的Cube
			tempI=(viter->x0 - GDS_zone.x0)/gridCellSize;
			ix1=max(((int)((viter->x0)*1e3) - (int)((GDS_zone.x0)*1e3)) % (int)(gridCellSize*1e3)==0? tempI-2:tempI-1 , 0  );

			tempI=int(ceil ((viter->x1 -  GDS_zone.x0)/(gridCellSize) ) ) +1;
			ix2=min(((int)((viter->x1)*1e3) - (int)((GDS_zone.x0)*1e3)) % (int)(gridCellSize*1e3)==0? tempI+1:tempI ,  nX -1  );

			tempI=( viter->z0 - GDS_zone.z0 )/gridCellSize;
			iz1=max(((int)((viter->z0)*1e3) - (int)((GDS_zone.z0)*1e3)) % (int)(gridCellSize*1e3)==0? tempI-2: tempI-1 , 0    );

			tempI=int(ceil ((viter->z1 - GDS_zone.z0)/(gridCellSize) ) ) +1;
			iz2=min(((int)((viter->z1)*1e3) - (int)((GDS_zone.z0)*1e3)) % (int)(gridCellSize*1e3)==0? tempI+1:tempI,  nZ -1  );

			for( int i= ix1;  i<ix2 ;i++  )
			{
				for(  int k=iz1; k<iz2; k++   )
				{
					candidateCheck2d(  *viter,  gridOctree2d[i][k]  );
				}
			}
		}
	}

	for( int i=0;  i< nX; i++   ){
		for( int k=0; k< nZ; k++ ){//this
			gridOctree2d[i][k].candidateListOfRootCell2d =gridOctree2d[i][k].candidateList2d;	
			gridOctree2d[i][k].candidateList2d.clear();				
			for(  SubConductorList2d::iterator  subCondIt2d=  gridOctree2d[i][k].candidateListOfRootCell2d.begin(); 
				subCondIt2d != gridOctree2d[i][k].candidateListOfRootCell2d.end(); subCondIt2d++   ){
				insertToOctree( *subCondIt2d ,  gridOctree2d[i][k] );//generateOctree
			}
		}
	}

// 	//改进方法：为八叉树创建索引数组（也称为2级三维数组）
	updateSubGridCellVec2d(gridOctree2d);
}

GridOctree2d partitionDomain( double cellSize,  int &nX, int &nZ ,Boundary GDS_zone){

	//inflate the original gds zone for security reason
	double x1 = GDS_zone.x0;   
	double x2 = GDS_zone.x1;
	// double y1=GDS::zone.y1;//min position (x0,y0,z0)
	// double y2=GDS::zone.y2;
	double z1=GDS_zone.z0;
	double z2=GDS_zone.z1;

	nX=  int(ceil(double(x2-x1)/double(cellSize)));
	nZ=  int(ceil(double(z2-z1)/double(cellSize)));
	
	GridOctree2d gridOctree2d;
	int i,k;
	for( i=0; i<nX; i++ ){
			vector < Cell2d > vecZ;
			for( k=0; k<nZ; k++  ){
				Cell2d cell2d=Cell2d(  double(x1+cellSize*i), double(x1+cellSize*(i+1)) , double(z1+cellSize*k), double(z1+cellSize*(k+1)) , double(cellSize) );
				vecZ.push_back(cell2d);
			}
		gridOctree2d.push_back(vecZ);
	}

	return  gridOctree2d;
}

double scale_of_extensionSize ;

void  insertToOctree( Subnet subCond, Cell2d &cellT  ){
	if(cellT.size()<1E-5){
		cout<<"Error: cellT is 2-d or 1-d. This program only handles 3-dimensional space"<<endl;
		return;
	}
	//by seting the initial value of distanceLimit, We can generate incomplete candidatelist
	if( cellT.distanceLimit !=-1 && myDistance2d( subCond,cellT )- cellT.distanceLimit>1E-7 ){
		return;
	}
	//if cellT is filled with conductor, then it's useless, omit it.
	if(  cellT.isFilledWithCond ){
		return;
	}

	//if cellT is a leaf node
	if(  cellT.subCellList.size()==0  ){
		candidateCheck2d( subCond, cellT );
		if( cellT.candidateList2d.size()>Cell2d::maxCandidateListLen && cellT.size() > Cell2d::minCellSize ){

		//In the future we may make some improvement about the extensionSize
			cellT.divideIntoSubCell(  cellT.extensionSize );
			for( SubConductorList2d::iterator candidateIt2d= cellT.candidateList2d.begin(); candidateIt2d != cellT.candidateList2d.end(); candidateIt2d++  ){
				for( SubCellList2d::iterator subCellIt2d=cellT.subCellList.begin(); subCellIt2d != cellT.subCellList.end(); subCellIt2d++ ){
					insertToOctree( *candidateIt2d, *subCellIt2d );
				}
			}
		}

	}else{
		for( SubCellList2d::iterator subCellIt2d=cellT.subCellList.begin(); subCellIt2d!= cellT.subCellList.end(); subCellIt2d++ ){
			insertToOctree( subCond, *subCellIt2d);
		}
	}
}

void candidateCheck2d(Subnet subCond,  Cell2d &cellT){
	if( cellT.isFilledWithCond ){
		cout<<"Error: try to do candidateCheck on a cellT which is filled with subConductor"<<endl;
		return;
	}
	if( cellT.isFilledChecked ==false ){
		if( cellT.isInsectWith(subCond) ){
			if(  cellT.isIncludedIn( subCond )  ){
				cellT.isFilledWithCond=true;
			}else{
				cellT.isFilledWithCond=false;
			}
			cellT.isFilledChecked=true;
		}
	}

	double dis=myDistance2d( subCond, cellT );
	double cellSize=cellT.size();
	//dis==distanceLimit, still need to recondisder rect, cannot omit it simplely
	if(  cellT.distanceLimit!=-1 && dis - cellT.distanceLimit>1E-7  ) {    
		return;
	}
	// if you conduct a delete or add operation on a list, you should break out the loop! or there will be wrong pointer!
	SubConductorList2d::iterator candidateRectIt2d;
	bool isContinue=false;
	// bool isContinue1=false;
	vector<Subnet>::iterator viter;
	// double x0,x1,z0,z1;

	while(true){
					isContinue=false;
		for (candidateRectIt2d=cellT.candidateList2d.begin();   candidateRectIt2d != cellT.candidateList2d.end(); candidateRectIt2d++ )
		{

			
				if(  dominate2d(  *candidateRectIt2d , subCond, cellT )  ){
					return;
				}
				else if(  dominate2d( subCond, *candidateRectIt2d, cellT )   ){
			
					cellT.candidateList2d.erase( candidateRectIt2d );
					isContinue=true;
					break;
				}
		}	
		if( isContinue)
		{
			continue;
		}
		if( candidateRectIt2d == cellT.candidateList2d.end() ){
			break;
		}
	}
	
//sorting the candidate conductors in the ascending order of their distance to the cell
	if(cellT.candidateList2d.size()==0){
		cellT.candidateList2d.push_back( subCond );
	}else{
		double disSubCond_CellT=myDistance2d(subCond, cellT);
		// SubConductorList::iterator candidateIt;
		SubConductorList2d::iterator candidateIt2d;
		for( candidateIt2d=cellT.candidateList2d.begin(); candidateIt2d != cellT.candidateList2d.end(); candidateIt2d++   ){
			if(   disSubCond_CellT - myDistance2d(  *candidateIt2d, cellT )<-1E-7  ){
				cellT.candidateList2d.insert( candidateIt2d, subCond  );
				break;
			}
		}
		if( candidateIt2d == cellT.candidateList2d.end() ){
			cellT.candidateList2d.push_back( subCond );
		}
	}

	
	if(    cellT.distanceLimit==-1 ||  (dis + cellSize ) - cellT.distanceLimit < -1E-7  ){
		cellT.distanceLimit=dis+cellSize;
	}

}

int max(int a, int b){
	return a>b? a:b;
}

double max( double a, double b  ){
	return a-b>1E-7 ?a:b;
}


int min( int a, int b  ){ 
	return a<b? a:b;
}

int min( int a1,int a2, int a3, int a4, int a5  , int a6  ){
	return min(min(min(min( min(a1,a2), a3), a4), a5),a6);

}

double min( double a1, double  a2){
	return a1-a2<-1E-7 ?a1:a2;
}

double min( double a1, double a2, double a3, double a4, double a5, double a6  ){
	return min(min(min(min( min(a1,a2), a3), a4), a5),a6);
}

double min( double a1, double a2, double a3, double a4 ){
	return min(min( min(a1,a2), a3), a4);
}

double myDistance2d( FPoint2d  &p, Boundary &rect ){
	double x_dis= max( rect.x0-p.x1, p.x1-rect.x1 );
	double z_dis=max(rect.z0-p.z1, p.z1-rect.z1 );
	return max(x_dis,z_dis);
}

double myDistance2d( FPoint2d  &p, Subnet rect  ){
	double x_dis= max( rect.x0-p.x1, p.x1-rect.x1 );
	double z_dis= max( rect.z0-p.z1, p.z1-rect.z1 );
	return max(x_dis,z_dis);
}
double myDistance2d(Subnet rectA, Cell2d &cellT ){
	double x_dis=max(   rectA.x0-cellT.x2, cellT.x1-rectA.x1);
	double z_dis=max(   rectA.z0-cellT.z2, cellT.z1-rectA.z1);
	// double temp=max(x_dis, y_dis );
	return max( x_dis,z_dis );

}

//这个函数用来计算两个subConductor之间的距离
//distance(rectA,rectB) represents the minimum distance of d( p in rectA, rectB  )
// int myDistance(Rectangle &rectA, Rectangle &rectB ){
// 	int x_dis=max(   rectA.x1- rectB.x2, rectB.x1-rectA.x2 );
// 	int y_dis=max(  rectA.y1-rectB.y2, rectB.y1-rectA.y2  );
// 	int z_dis=max(   rectA.z1-rectB.z2,  rectB.z1-rectA.z2 );
// 	int temp=max(x_dis, y_dis );
// 	return max( temp,z_dis );
// }

int myDistance2d(Subnet rectA ,Subnet rectB){
	int x_dis=max(   rectA.x0 - rectB.x1, rectB.x0 - rectA.x1 );
	int z_dis=max(   rectA.z0 - rectB.z1, rectB.z0 - rectA.z1 );
	return max( x_dis,z_dis );
}
bool dominate2d(Subnet rectA ,Subnet rectB, Cell2d &cellT)
{

	int nVertex=8;
	int i=0;
	double disTA=0;
	double disTB=0;
	double maxDistanceTA=0;
	double minDistanceTB=0;
	vector<FPoint2d> vertex(nVertex);
	//counter-clock

	vertex[0]=FPoint2d(cellT.x1, cellT.z1 );   
	vertex[1]=FPoint2d(cellT.x2, cellT.z1 );
	vertex[2]=FPoint2d(cellT.x1, cellT.z2 );   
	vertex[3]=FPoint2d(cellT.x2, cellT.z2 );
	if(cellT.isInsectWith(rectB) ){

		return false;                           
	}

	for(i=0;i<nVertex;i++){
		disTA=myDistance2d(vertex[i], rectA );
		disTA=disTA<0?0:disTA;
		if(disTA - maxDistanceTA>1E-7 ){
			maxDistanceTA=disTA;
		}
	}

	minDistanceTB=myDistance2d( rectB, cellT );
	//if maxDistanceTA < minDistanceTB , rectA will dominate rectB with respect to rectT
	return (maxDistanceTA-minDistanceTB<-1E-7);

}

void Cell2d::divideIntoSubCell( double extension){
	int nSubCell=4;

	double centerX= (x1+x2)/2.0;
	double centerZ= (z1+z2)/2.0;

	centerPoint.setCoordinate(centerX, centerZ );
	hasChildCell=true;
	subCellList.push_back(Cell2d( x1, centerX,  z1, centerZ ,extension ));
	subCellList.push_back(Cell2d( centerX, x2,  z1, centerZ ,extension ));
	subCellList.push_back(Cell2d( x1, centerX,  centerZ,z2 ,extension ));
	subCellList.push_back(Cell2d( centerX, x2,  centerZ,z2 ,extension ));
}

void updateSubGridCellVec2d(Cell2d &rootCell,  Cell2d *cellP,  vector<vector<Cell2d*> > &subGridCellVec ){
	if(!cellP->hasChildCell || cellP->size()-Cell2d::subGridSize<1E-5  ){
		int startXn=(int) nearbyint((cellP->x1-rootCell.x1)/Cell2d::subGridSize);
		int startZn=(int) nearbyint((cellP->z1-rootCell.z1)/Cell2d::subGridSize);
		int endXn=((int) nearbyint((cellP->x2-rootCell.x1)/Cell2d::subGridSize))-1 ;
		int endZn=((int) nearbyint((cellP->z2-rootCell.z1)/Cell2d::subGridSize))-1;
		if( fabs((cellP->x1-rootCell.x1)/Cell2d::subGridSize - startXn) > 1E-5 ||  
			fabs((cellP->z1-rootCell.z1)/Cell2d::subGridSize - startZn) > 1E-5 ||  
			fabs((cellP->x2-rootCell.x1)/Cell2d::subGridSize -1 - endXn) > 1E-5 ||  
			fabs((cellP->z2-rootCell.z1)/Cell2d::subGridSize -1 - endZn) > 1E-5   ){
			cout<<"Err: wrong Cell location!"<<endl;
			exit(1);
		}

		for(int i=startXn; i<= endXn; i++ ){
				for(int k=startZn; k<=endZn; k++ ){
					subGridCellVec[i][k]=cellP;
				}
			
		}
	}
	else{
		for(SubCellList2d::iterator subCellIt=cellP->subCellList.begin(); subCellIt!= cellP->subCellList.end(); subCellIt++   ){
			updateSubGridCellVec2d( rootCell, &(*subCellIt), subGridCellVec );
		}

	}
}

void updateSubGridCellVec2d(GridOctree2d &gridOctree)
{
	for(int i=0; i<gridOctree.size(); i++ )
	{
		for(int k=0; k<gridOctree[i].size(); k++ )
		{
			int subGridNum=(int) nearbyint( GridOctree2d::gridCellSize/Cell2d::subGridSize );
			Cell2d *cellP=NULL;
			std::vector<Cell2d*> vz(subGridNum, cellP);
			std::vector<std::vector<Cell2d*> > vy(subGridNum, vz );
			updateSubGridCellVec2d(gridOctree[i][k], &(gridOctree[i][k])  ,  vy );
			gridOctree[i][k].subGridCellVec=vy;
		}
	}
}
