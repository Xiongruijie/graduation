#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <pthread.h>
#include <vector>
using namespace std;

int numOfGrid=2000;

vector<string> outNameList(3);

double unitSurfaceGreenFunction(  double ni, int N  );
// void *generateGreenFuncTable(void * ind);
// void copy(ofstream &fout, string &filename );
void GF(int numOfGrid);
void GX(int numOfGrid);
void GZ(int numOfGrid);

double unitGEx( double ni, int N  );
double unitGEz( double ni, int N  );

int main(int argc, char *argv[] ){

	numOfGrid=atoi(argv[1]);
	GF(numOfGrid);
	GX(numOfGrid);
	GZ(numOfGrid);

}

void GF(int numOfGrid){

	string outfileName = "Green.txt";
	ofstream fout(outfileName);

	if(!fout.is_open() ){
		cout<<"file open failed!"<<endl;
		exit(1);
	}
	for (int i = 0 ; i < numOfGrid; i++)
	{	

		double ni = 0;
		double green = 0;
		int k=0;

		while(k++<10){
			ni=1.0/double(numOfGrid)*i + (rand()/(double(RAND_MAX)))/numOfGrid;
			green+=unitSurfaceGreenFunction(ni,10);
		}
		green=green/10.0;
		fout<<setprecision(9)<<green<<"\t";
	}
	fout.close();
}
void GX(int numOfGrid){

	string outfileName = "GEx.txt";
	ofstream fout(outfileName);


	if(!fout.is_open() ){
		cout<<"file open failed!"<<endl;
		exit(1);
	}
	for (int i = 0 ; i < numOfGrid; i++)
	{
		double ni = 0;
		double green = 0;
		int k=0;
		
		while(k++<10)
		{
			ni=1.0/double(numOfGrid)*i + (rand()/(double(RAND_MAX)))/numOfGrid;
			green+=unitGEx(ni,10);
		}
		green=green/10.0;
		fout<<setprecision(9)<<green<<"\t";
	}
	fout.close();
}
void GZ(int numOfGrid){

	string outfileName = "GEz.txt";
	ofstream fout(outfileName);

	if(!fout.is_open() ){
		cout<<"file open failed!"<<endl;
		exit(1);
	}
	for (int i = 0 ; i < numOfGrid; i++)
	{	
		int k = 0;
		double ni = 0;
		double green = 0;

		while(k++<10)
		{
			ni=1.0/double(numOfGrid)*i + (rand()/(double(RAND_MAX)))/numOfGrid;
			green+=unitGEz(ni,10);
		}
		green=green/10.0;
		fout<<setprecision(9)<<green<<"\t";
	}
	fout.close();
}




double unitSurfaceGreenFunction(  double ni, int N )
{//ni and n
	int n;
	double greenFuncValue=0;
	for( n=1; n<=N; n++ ){//n
			greenFuncValue+=2.0 *sin(M_PI*n/2.0)*sinh(M_PI*n/2.0)* sin(M_PI*n*ni)/sinh(M_PI*n);
	}
	return greenFuncValue;
}

double unitGEx( double ni, int N  )
{
	int n;
	double GE=0;
	for( n=1; n<=N; n++ ){
		// GE+= (-1.0)*M_PI * n * cos(M_PI*n/2.0)*sinh(M_PI*n/2.0)/sinh(M_PI*n)*sin(M_PI*n*ni);
		GE+= (-1.0)*2*M_PI * n * cos(M_PI*n/2.0)*sinh(M_PI*n/2.0)/sinh(M_PI*n)*sin(M_PI*n*ni);

	}//
	return GE;
}

double unitGEz(  double ni, int N )
{
	int n;
	double GE=0;
	for(n=1; n<=N; n++ ){
		// GE+= (-1.0)*M_PI * n * sin(M_PI*n/2.0)*cosh(M_PI*n/2.0)/sinh(M_PI*n) *sin(M_PI*n*ni);
		GE+= (-1.0)*2*M_PI * n * sin(M_PI*n/2.0)*cosh(M_PI*n/2.0)/sinh(M_PI*n) *sin(M_PI*n*ni);

	}		
	return GE;
}
