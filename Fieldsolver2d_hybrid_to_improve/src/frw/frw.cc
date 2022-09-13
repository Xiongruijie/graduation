#ifndef COMMON_IS_INCLUDED
#endif

#include "FRWcapext.cpp"
using namespace std;


int frw(map<int, Net*> generateConductorList, Boundary GDS_zone, double dielectric)
{

	cout << endl;
    cout << "============================================" << endl;
    cout << ">                frw is used               <" << endl;
    cout << "============================================" << endl;
    cout << endl;


	int NTHREAD = 4;
	FRWControl::maxPermittedCapErr = 1;
	long loopNum = 10000;
 	FRWControl::dielectricConstant = dielectric*0.008854 ;    //e0

	GDS_zone.x0 *= 1000;
	GDS_zone.z0 *= 1000;
	GDS_zone.x1 *= 1000;
	GDS_zone.z1 *= 1000;

	map<int, Net*>::iterator nets_iter;
	vector<Subnet>::iterator subnets_iter;
	for (nets_iter = generateConductorList.begin(); nets_iter != generateConductorList.end(); ++nets_iter) {
		nets_iter->second->x0 *= 1000;
		nets_iter->second->z0 *= 1000;
		nets_iter->second->x1 *= 1000;
		nets_iter->second->z1 *= 1000;
		for (subnets_iter = nets_iter->second->subnets.begin(); subnets_iter != nets_iter->second->subnets.end(); ++subnets_iter) {
			subnets_iter->x0 *= 1000;	
			subnets_iter->z0 *= 1000;	
			subnets_iter->x1 *= 1000;	
			subnets_iter->z1 *= 1000;	
		}
    }

	//print_input(GDS_zone, dielectric, generateConductorList);


// int main(int argc, char* argv[]  ){


// 	int NTHREAD= atoi(argv[2]) ;
// 	FRWControl::maxPermittedCapErr=atof(argv[3]); //最大容许的电容误差百分比
// 	long loopNum= atol(argv[4]);
/*
	stringstream ss;//output
// 	ss<<argv[5];
	string outFileName=ss.str();
	ss.str("capOut.txt");
	// ss<<argv[1];
	string geoFileName=ss.str();
*/

gettimeofday(&tv_begin,NULL );
	// map<string, Net*> generateConductorList = GenerateConductorList(geoFileName);
	// Boundary GDS_zone = GDS(geoFileName);//why it is stop for hear
	FRWControl::conductorList2d = generateConductorList;
	convertcondlist(generateConductorList);
gettimeofday(&tv_end,NULL);
cout<<"ConductorList generated.----- Time consumption: "<<( (tv_end.tv_sec*1000000+tv_end.tv_usec)-(tv_begin.tv_sec*1000000+tv_begin.tv_usec))/1000<<" ms" <<endl;

gettimeofday(&tv_begin,NULL );
	loadDataTable("DataTables_test/Green.txt", FRWControl::GreenVT ,FRWControl::GreenVTSize );
	loadDataTable("DataTables_test/GEx.txt", FRWControl::GExVT , FRWControl::GExVTSize );
	loadDataTable("DataTables_test/GEz.txt", FRWControl::GEzVT , FRWControl::GEzVTSize );
gettimeofday(&tv_end,NULL);
cout<<"Data table loaded. ----- Time consumption: "<<( (tv_end.tv_sec*1000000+tv_end.tv_usec)-(tv_begin.tv_sec*1000000+tv_begin.tv_usec))/1000<<" ms" <<endl;

gettimeofday(&tv_begin,NULL );
	generateGridOctree(FRWControl::gridOctree2d, GridOctree2d::gridCellSize ,GDS_zone ,generateConductorList,FRWControl::teatlist);
gettimeofday(&tv_end,NULL);
cout<<"GridOctree generated.----- Time consumption: "<<( (tv_end.tv_sec*1000000+tv_end.tv_usec)-(tv_begin.tv_sec*1000000+tv_begin.tv_usec))/1000<<" ms" <<endl;

gettimeofday(&tv_begin,NULL );
	double scale_factor=2;
	generateGaussianSurfaceOfConductor(FRWControl::teatlist, scale_factor);//change area
gettimeofday(&tv_end,NULL);
cout<<"VGS generated. ----- Time consumption: "<<( (tv_end.tv_sec*1000000+tv_end.tv_usec)-(tv_begin.tv_sec*1000000+tv_begin.tv_usec))/1000<<" ms" <<endl;
	
	FRWControl::capMatrix= generateCapMatrix( FRWControl::teatlist );
	FRWControl::totalNumOfWalkOfCond = vector<long>(FRWControl::teatlist.testlist.size(), loopNum );
	FRWControl::currentNumOfWalkOfCond =vector<long>(FRWControl::teatlist.testlist.size(),0 );
	FRWControl::currentProgressOfCond = vector<long>(FRWControl::teatlist.testlist.size(),0);
	cout<<"Global variable initialized!"<<endl;

	srand( unsigned (time(NULL)) );
	FRWControl::isRandSeedSet=true;

//     //generate a random seed array
//     //利用boost Random库为每一个线程生成一组随机数生成器，randSeedVec为其提供初始化种子
	int maxNumOfThread=500;
	for(int r=0;r<maxNumOfThread; r++ ){
		long rv=rand();
		bool isSame=false;
		for(int k=0; k<FRWControl::randSeedVec.size(); k++){
			if( FRWControl::randSeedVec[k]==rv ){
				isSame=true;
			}
		}
		if(!isSame){
			FRWControl::randSeedVec.push_back(rv);
		}
	}


	cout<<"FRW walk started."<<endl;
	//start counting
gettimeofday(&tv_begin,NULL );

	pthread_t tids[NTHREAD];
	
	pthread_attr_t attr;
	pthread_attr_init(&attr );
	pthread_attr_setdetachstate( &attr, PTHREAD_CREATE_JOINABLE );

	pthread_mutex_init( &getRandSeed_mutex, NULL );
	pthread_mutex_init( &updateFRWData_mutex, NULL );
	pthread_mutex_init(&processorIndex_mutex, NULL);

	for( int i=0; i<NTHREAD; i++ ){
		if(    pthread_create(&tids[i], &attr, FRW, NULL ) !=0  ){
			cout<< " pthread_create ERROR!!"<<endl;
			exit(1);
		}
		//使线程错开，尽量避免同时访问锁
		sleep(1);
	}
	pthread_attr_destroy( &attr );
	void *status;
	for( int i=0; i<NTHREAD; i++ ){
		if(    pthread_join(  tids[i], &status )  !=0 ){
			cout<<"pthread_join  ERROR!!"<<endl;
			exit(1);
		}
	}

	pthread_mutex_destroy( &getRandSeed_mutex );
	pthread_mutex_destroy( &updateFRWData_mutex);
	pthread_mutex_destroy(&processorIndex_mutex);

	//end counting
gettimeofday(&tv_end,NULL);
cout<<"\nFRW walk ended. ----- Time consumption: "<<(tv_end.tv_sec-tv_begin.tv_sec)<<" s" <<endl;

	cout<<"locateCellTime: "<<locateCellTime/1000000<<endl;
	cout<<"generateMaxCubeTime: "<<generateMaxCubeTime/1000000<<endl;
	cout<<"generatePointTime: "<<generatePointTime/1000000<<endl;

	cout<<"hitTime: "<<hitTime<<endl;
	cout<<"totalStepNum: "<<totalStepNum<<endl;
	cout<<"hit rate: "<<double(hitTime)/totalStepNum<<endl;


	FRWControl::teatlist.test_iter = FRWControl::teatlist.testlist.begin();
    
	for(int column=0; column<FRWControl::capMatrix.size(); column++ ){
		FRWControl::capMatrix[column].capacitance=FRWControl::capMatrix[column].sumOfCap*FRWControl::teatlist.test_iter->conductor2d.gaussianSurfaceList.VGSArea ;
		FRWControl::capMatrix[column].capacitance/=FRWControl::currentNumOfWalkOfCond[0];
	}

	/*
	ofstream capFout(outFileName.c_str());

	for(int j=0; j<FRWControl::capMatrix.size(); j++ ){
		cout<<FRWControl::capMatrix[j].capacitance;
		capFout<<FRWControl::capMatrix[j].capacitance;

		if( j<FRWControl::capMatrix.size()-1 ){
			cout<<"\t"<<flush;
			capFout<<"\t"<<flush;
		}
	}
	cout<<endl;
	capFout<<"\r\n";

	capFout.close();
	*/
	int j = 0;
//	map<int, Net*>::iterator nets_iter;

	for (nets_iter = generateConductorList.begin(); nets_iter != generateConductorList.end(); ++nets_iter,j++) {

		cout << "net" << nets_iter->first << " : " << fixed << setprecision(6) << fabs(FRWControl::capMatrix[j].capacitance) << " fF" << endl;
		// restore cout precision
        cout <<resetiosflags(ios::fixed) <<setprecision(6);
		// Save it into data structure of each net
		nets_iter->second->set_cap(fabs(FRWControl::capMatrix[j].capacitance));
	}

	cout << endl;
    cout << "============================================" << endl;
    cout << "*               frw succeeded!             *" << endl;
    cout << "============================================" << endl;
    cout << endl;


	return 1;
	
}

