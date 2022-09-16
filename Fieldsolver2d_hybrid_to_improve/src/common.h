#include <iostream>
#include <cstring>   // string
#include <fstream>   // ifstream
//#include <regex>     // for split()
#include <vector>    // vector
#include <algorithm> // sort
#include <cmath>     // fabs
#include <map>       // map
#include <iomanip>   // setprecision setw
#include <sys/resource.h> // getrusage
#include <list>      // list
#include <sstream>   // istringstream, for split_2()
#include <getopt.h>  // getopt_long

#include <stdexcept>
#include <functional>
#include <string>
#include <sys/time.h>
#include <unistd.h>
#include <boost/random.hpp>
#include <pthread.h>

#define COMMON_IS_INCLUDED

using namespace std;

// Double type comparison
const double eps = 1e-6;
bool dbl_equ(double src, double dst) {
	return fabs(src-dst) < eps;
}

// Used for find_if()
struct dbl_cmp {
    dbl_cmp(double v, double d) : val(v), delta(d) { }
    inline bool operator()(const double &x) const {
        return abs(x-val) < delta;
    }
private:
    double val, delta;
};

// Define empty space dielectric
const double e0 = 0.008854;

// Define class for boundary
class Boundary
{
    public:
		// these are the original boundary coordinates
        double x0, x1, z0, z1;
        double width;
        double height;
        void set(double a, double b, double c, double d);
		// these are the minimum rectangle coordinates containing all conductors
		double rec_x0, rec_x1, rec_z0, rec_z1;
		void shrink(double shrink_val);
};

void Boundary::set (double a, double b, double c, double d)
{
    x0 = a;
    z0 = b;
    x1 = c;
    z1 = d;
    width = x1 - x0;
    height = z1 - z0;
	/*
	cout << "x0:"<< x0 <<endl;
	cout << "x1:"<< x1 <<endl;
	cout << "z0:"<< z0 <<endl;
	cout << "z1:"<< z1 <<endl;
	cout << "width:"<< width <<endl;
	cout << "height:"<< height <<endl;
	*/
    if (width <= 0 || height <= 0) {
        cout << "Error boundary coordinates" << endl;
        exit(1);
    }
}

void Boundary::shrink (double shrink_val)
{
	if (shrink_val > eps) {
		cout << "Original boundary: " << x0 << " ";
		cout << z0 << " ";
		cout << x1 << " ";
		cout << z1 << endl;

		cout << "The minimum rect: " << rec_x0 << " ";
		cout << rec_z0 << " ";
		cout << rec_x1 << " ";
		cout << rec_z1 << endl;

		x0 = rec_x0 - shrink_val;
		z0 = rec_z0 - shrink_val;
		x1 = rec_x1 + shrink_val;
		z1 = rec_z1 + shrink_val;

		cout << "After shrinking: " << x0 << " ";
		cout << z0 << " ";
		cout << x1 << " ";
		cout << z1 << endl;
	} 
}

// Define class for net
class Subnet
{
    public:
        double x0,x1,z0,z1;
		double gauss_x0, gauss_x1, gauss_z0, gauss_z1;
        double width;
        double height;
        int net_name;
        Subnet(double a, double b, double c, double d, int Net_name);
        bool operator == (Subnet& Subnet );
};

Subnet::Subnet(double a, double b, double c, double d, int Net_name)
{
    x0 = a;
    z0 = b;
    x1 = c;
    z1 = d;
    net_name = Net_name;
    width = x1 - x0;
    height = z1 - z0;
    if (width <= 0 || height <= 0) {
        cout << "Error net coordinates" << endl;
        exit(1);
    }
}

// Define net border
class Border
{
	public: 
		int flag;     // 0 - bottom to up; 1 - left to right; 2 - up to bottom; 3 right to left
		double common; // common coordinate between the two border points; it is the x-coordinate if dir is 0
		double start; // start point coordinate; it is the y-coordinate of the bottom point if dir is 0
		double end;   // it is the y-coordinate of the up point if dir is 0
		Border(int f, double c, double s, double e);
};
Border::Border(int f, double c, double s, double e)
{
	flag = f;
	common = c;
	start = s;
	end = e;
}

class Net
{
    public:
		// these are the coordinates of the smallest rectangle containing all subnets
		// they are the corresponding indexes of these coordinates in x_coords and z_coords
		// these parameters are used to work out # of grid points on nets
        double x0, z0, x1, z1; 
        unsigned int lo_x, lo_z, hi_x, hi_z; 

		// gauss curve exactly surrounding the net borders for compute charge
		list<Border> borders; 
		unsigned long int seq_start, seq_length; // sequential number in those element vectors such as xk, zk, lk etc.
		void set_seqlen(unsigned long int start, unsigned int length);

        double cap;  // net cap, the final result

        vector<Subnet> subnets; // for debug
        void insert(double a, double b, double c, double d, int Net_name);   // insert other subnets and update x0, x1, z0, z1

        void set_index(unsigned int a, unsigned int b, unsigned int c, unsigned int d); // set lo_x in terms of x0 through searching x_coords
        void set_cap(double c); // set cap
        Net(double a, double b, double c, double d, int Net_name); // insert the 1st subnet when initializing
};

Net::Net(double a, double b, double c, double d, int Net_name)
{
    x0 = a;
    z0 = b;
    x1 = c;
    z1 = d;
    Subnet subnet(a, b, c, d ,Net_name);
    subnets.push_back(subnet);
	Border border_l(0, a, b, d); // left
	Border border_u(1, d, a, c); // up
	Border border_r(2, c, d, b); // right
	Border border_b(3, b, c, a); // bottom
	borders.push_back(border_l);
	borders.push_back(border_u);
	borders.push_back(border_r);
	borders.push_back(border_b);
}

void Net::insert(double a, double b, double c, double d, int Net_name)
{
    x0 = min(a, x0); // update x0, z0, x1, z1 only if the subnet affects the outer rectangle boundary
    z0 = min(b, z0);
    x1 = max(c, x1);
    z1 = max(d, z1);
    Subnet subnet(a,b,c,d,Net_name);
    subnets.push_back(subnet);

	int up_overlap, bottom_overlap;

	up_overlap = 0;
	bottom_overlap = 0;

	// only up and bottom borders may overlap with prior nets	
	for(list<Border>::iterator iter=borders.begin();iter!=borders.end();iter++){ 
	//	cout << "iter->flag " << iter->flag << endl;
	//	cout << "iter->common " << iter->common << endl;
	//	cout << "iter->start " << iter->start << endl;
	//	cout << "iter->end " << iter->end << endl;
		if (dbl_equ(iter->common, d)) {
			up_overlap = 1;
			Border border_new_r(3, d, iter->end, c);   // only keep the two non-overlap parts
			Border border_new_l(3, d, a, iter->start); // of the prior bottom 
			borders.insert(iter, border_new_r);
			borders.insert(iter, border_new_l);
       		iter = borders.erase(iter);
			iter--;
		} else if (dbl_equ(iter->common, b)) {  // bottom
			bottom_overlap = 1;
			Border border_new_r(3, b, c, iter->end);
			Border border_new_l(3, b, iter->start, a);
			borders.insert(iter, border_new_r);
			borders.insert(iter, border_new_l);
       		iter = borders.erase(iter);
			iter--;
		} 
	}
	Border border_l(0, a, b, d); // left
	Border border_r(2, c, d, b); // right
	borders.push_back(border_l);
	borders.push_back(border_r);
	if (up_overlap == 0) {
		Border border_u(1, d, a, c); // up
		borders.push_back(border_u);
	}
	if (bottom_overlap == 0) {
		Border border_b(3, b, c, a); // bottom
		borders.push_back( border_b);
	}
}

void Net::set_index(unsigned int a, unsigned int b, unsigned int c, unsigned int d)
{
    lo_x = a;
    lo_z = b;
    hi_x = c;
    hi_z = d;
}

void Net::set_seqlen(unsigned long int start, unsigned int length)
{
	seq_start = start;
	seq_length = length;
}

void Net::set_cap(double c)
{
    cap = c;
}

// Split string in terms of delim
/*
std::vector<std::string> split(const std::string& input,
                               const std::string& delim)
{
    std::regex re(delim);
    std::sregex_token_iterator first {input.begin(), input.end(), re, -1}, last;
    return {first, last};
}
*/
void split_2(const std::string& input, const std::string& delim, vector<string>& res)
{
	istringstream is(input);
	std::string s;
	if (delim != " ") {
		cout <<"Error: only space delim is supported"<< endl;
		exit(1);
	}
	while(is>>s)
		res.push_back(s);
}

// Sort and unique vector
void sort_unique(std::vector<double> &v)
{
    vector<double>::iterator uniq_v;
    std::sort(v.begin(),v.end());
    uniq_v = unique(v.begin(), v.end());
    v.erase(uniq_v, v.end());
}

// Print vector
//void print_vector(std::vector<double> v)
//{
//    for (std::vector<double>::iterator it=v.begin(); it!=v.end(); ++it) {
//        std::cout << *it << ' ';
//    }
//    copy(v.begin(), v.end(), ostream_iterator<double>(cout, " "));
//    std::cout<<endl;
//}

void print_vector(std::vector<double> v, int flag=0)
{
	double prev = 0;
	double min = +100;
	double max = -100;
	double diff;
	double min_coord;
    for (std::vector<double>::iterator it=v.begin(); it!=v.end(); ++it) {
        std::cout << *it << ' ';
		if (flag) {
			// the first element
			if (prev < eps) {
				prev = *it;
			} else {
				diff = *it - prev;
				if (min > diff + eps) {
					min = diff;
					min_coord = *it;
				}
				if (max + eps < diff) {
					max = diff;
				}
				prev = *it;
			}
		}
    }
    /*
    copy(v.begin(), v.end(), ostream_iterator<double>(cout, " "));
    */
    std::cout<<endl;
	if (flag) {
		cout << "The minimum distance is: "<< min << " at " << min_coord << endl;
		cout << "The maximum distance is: "<< max << endl;
	}
}

// Print usage in case of wrong cmdline parameters
void print_usage()
{
    cout << "Usage: fieldsolver2d -in input.data -out result.out ";
	cout << "[-method auto|fdm|bem|frw]" << endl << endl;;
	cout << "NOTICE: -method is optional. By default, auto is used and the program will determine one method automatically in terms of the scenario. ";
	cout << "If one method is specified by this option, only this method will be used even if it would fail. ";
	cout << "Here fdm represents finite difference method, bem stands for boundary element method and frw means floating random walking." << endl;
}

const std::string helpInfo = 
R"(
Usage: fieldsolver2d --in input.data --out result.out [options]
  --help            Display this information
  --shrink <value>  Shrink boundary in terms of the locations of all conductors and this value has to be greater than 0 
  --method <auto|fdm|bem|frw> Specify method. By default auto is used and it will determine one automatically in terms of the input data
)";

const struct option LongOptions[] =
{
	{"in", required_argument, NULL, 'i'},
	{"out", required_argument, NULL, 'o'},
	{"method", required_argument, NULL, 'm'},
	{"shrink", required_argument, NULL, 's'},
	{"help", no_argument, NULL, 'h'},
	{0,0,0,0}
};

int parse_options(int argc, char* argv[], string &input_data, string &result_out, string &use_method, double &shrink_val)
{
	int error_code = 0;
	int opt = 0;
	int optIndex = 0;
	string shrink_str;

	input_data = "";
	result_out = "";
	shrink_val = 0;
	use_method = "auto";
	opt = getopt_long(argc, argv, "i:o:m:s:h", LongOptions, &optIndex);

//	if (opt == -1) {
//		std::cout << helpInfo << std::endl;
//		error_code = 1;
//	}

	while (-1 != opt) {
		switch (opt) 
		{
			case 'i':
				input_data = optarg; 
				break;

			case 'o':
				result_out = optarg; 
				break;

			case 'h':
				std::cout << helpInfo << std::endl;
				error_code = 1;
				break;

			case 'm':
				use_method = optarg;
				break;

			case 's':
				shrink_str = optarg;
				shrink_val = atof(shrink_str.c_str());
				if (shrink_val < eps) {
					cout << "shrink value has to be greater than 0" << endl;
					error_code = 1;
				}
				break;

			default:
				std::cout << helpInfo << std::endl;
				error_code = 1;
				break;
		}
		opt = getopt_long(argc, argv, "i:o:m:s:h", LongOptions, &optIndex);
	}

	// these two arguments are must
	if (input_data == "" || result_out =="") {
		std::cout << helpInfo << std::endl;
		error_code = 1;
	}
	return error_code;
}

// Print runtime and memory usage
void print_runtime()
{
    cout << "------------- Runtime and Memory -----------" << endl;
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    long user = usage.ru_utime.tv_sec * 1000000 + usage.ru_utime.tv_usec; // user time used
    long sys  = usage.ru_stime.tv_sec * 1e6 + usage.ru_stime.tv_usec; // sys time used
    long mem  = usage.ru_maxrss;
    cout << "User time  : " << setw(10) << user << " us" << endl;
    cout << "Sys time   : " << setw(10) << sys << " us" << endl;
    cout << "Total time : " << setw(10) << user+sys << " us" << endl;
    cout << "Max memory : " << setw(10) << mem << " kB" << endl;
}

// Print input data
//void print_input(Boundary b, double d, map<string,Net*> nets)
void print_input(Boundary b, double d, map<int,Net*> nets)
{
    //map<string, Net*>::iterator nets_iter;
    map<int, Net*>::iterator nets_iter;

    cout << "------------------Dielectric----------------" << endl;
    cout << "Value: " << d << "   Epislon0: " << e0 << " ff/um" << endl;
    cout << "------------------ Boundary ----------------" << endl;
    cout << "Width: " << b.width << " Height: " << b.height;
    cout << " Coord: " << b.x0 << " " << b.z0 << " " << b.x1 << " " << b.z1 << endl;
    cout << "-------------------- Nets ------------------" << endl;
    cout << "#Nets: " << nets.size() << endl;
    for (nets_iter = nets.begin(); nets_iter != nets.end(); ++nets_iter) {
        //cout << "Name: " << nets_iter->first;
        cout << "Name: net" << nets_iter->first;
        cout << " Rect: " << nets_iter->second->x0 << " "
                          << nets_iter->second->z0 << " "
                          << nets_iter->second->x1 << " "
                          << nets_iter->second->z1 << endl;
    }
    cout << "------------------- Subnets ----------------" << endl;
    for (nets_iter = nets.begin(); nets_iter != nets.end(); ++nets_iter) {
        vector<Subnet> subnets = nets_iter->second->subnets;
        for (unsigned int j=0; j<subnets.size(); j++) {
            //cout << "Name: " << nets_iter->first;
            cout << "Name: net" << nets_iter->first;
            cout << " Width: " << subnets[j].width
                 << " Height: "<< subnets[j].height;
            cout << " Coord: " << subnets[j].x0 << " "
                               << subnets[j].z0 << " "
                               << subnets[j].x1 << " "
                               << subnets[j].z1 << endl;
        }
    }
}

// Print gauss border list for each net
//void print_gauss(map<string, Net*> nets)
void print_gauss(map<int, Net*> nets)
{
    //map<string, Net*>::iterator nets_iter;
    map<int, Net*>::iterator nets_iter;
    cout << "-------------------- Gauss -----------------" << endl;
    for (nets_iter = nets.begin(); nets_iter != nets.end(); ++nets_iter) {
        //cout << "Name: " << nets_iter->first << endl;
        cout << "Name: net" << nets_iter->first << endl;
//		cout << nets_iter->second->borders.size() << endl;
		for (list<Border>::iterator iter = nets_iter->second->borders.begin(); 
				iter != nets_iter->second->borders.end();
					++iter) {
			switch(iter->flag) {
				case 0:
					cout << "V: (" << iter->common << ", " << iter->start << "), (" << iter->common << ", " << iter->end << ")" << endl;
					break;
				case 1:
					cout << "H: (" << iter->start << ", " << iter->common << "), (" << iter->end << ", " << iter->common << ")" << endl;
					break;
				case 2:
					cout << "V: (" << iter->common << ", " << iter->start << "), (" << iter->common << ", " << iter->end << ")" << endl;
					break;
				case 3:
					cout << "H: (" << iter->start << ", " << iter->common << "), (" << iter->end << ", " << iter->common << ")" << endl;
					break;
				default:
					cout << "wrong flag: " << iter->flag << endl;
					exit(1);
			}
		}
		/*
		list<Border>::iterator bdr_iter;
		while(bdr_iter != nets_iter->second->borders.end()) {
			cout << bdr_iter->flag << endl;
		}
		*/
	}
}

// Print result to result file
// b represents bits number after dot, i.e. precision
//void print_result(map<string, Net*> nets, string result_out, int b)
/*
void print_result(map<int, Net*> nets, string result_out, int b)
{
    cout << "--------------- Output Result --------------" << endl;
    //map<string, Net*>::iterator nets_iter;
    map<int, Net*>::iterator nets_iter;
    ofstream fout(result_out);

    cout << setiosflags(ios::left) << setw(5) << setfill(' ') << " " << resetiosflags(ios::left);
    fout << setiosflags(ios::left) << setw(5) << setfill(' ') << " " << resetiosflags(ios::left);

    for (nets_iter = nets.begin(); nets_iter != nets.end(); ++nets_iter) {
        cout << setiosflags(ios::right) << setw(b+6) << setfill(' ') << setprecision(b) << "net" << nets_iter->first;
        fout << setiosflags(ios::right) << setw(b+6) << setfill(' ') << setprecision(b) << "net" << nets_iter->first;
    }    

    cout << resetiosflags(ios::right) << endl;
    fout << resetiosflags(ios::right) << endl;

    cout << setiosflags(ios::left) << setw(5) << setfill(' ') << "net0:" << resetiosflags(ios::left);
    fout << setiosflags(ios::left) << setw(5) << setfill(' ') << "net0:" << resetiosflags(ios::left);

    // b+4 being 2 less than b+6 is due to the tail 'ff'
    for (nets_iter = nets.begin(); nets_iter != nets.end(); ++nets_iter) {
        cout << setiosflags(ios::right) << setw(b+4) << setfill(' ') << setprecision(b) << nets_iter->second->cap <<"fF";
        fout << setiosflags(ios::right) << setw(b+4) << setfill(' ') << setprecision(b) << nets_iter->second->cap <<"fF";
        // delete space allocated by new
        delete nets_iter->second;
    }    
    // restore cout precision
    cout << resetiosflags(ios::right) << setprecision(6) << endl;
    fout << resetiosflags(ios::right) << setprecision(6) << endl;

    fout.close();
}
*/

void print_result_2(map<int, Net*> nets, string result_out, int b)
{
    cout << "--------------- Output Result --------------" << endl;
    //map<string, Net*>::iterator nets_iter;
    map<int, Net*>::iterator nets_iter;
    ofstream fout(result_out);

	cout << "# length = 1um" << endl;
	fout << "# length = 1um" << endl;

	cout << "# master = net0" << endl;
	fout << "# master = net0" << endl;

    for (nets_iter = nets.begin(); nets_iter != nets.end(); ++nets_iter) {
        cout << "          " << "net" << nets_iter->first;
        fout << "          " << "net" << nets_iter->first;
    }    

	cout << endl;
	fout << endl;

    cout << " net0:"; 
    fout << " net0:";

    for (nets_iter = nets.begin(); nets_iter != nets.end(); ++nets_iter) {
		if (nets_iter->first == 0) {
			cout << " " << nets_iter->second->cap  <<"fF";
			fout << " " << nets_iter->second->cap  <<"fF";
		} else {
			cout << "   " << nets_iter->second->cap  <<"fF";
			fout << "   " << nets_iter->second->cap  <<"fF";
		}
        // delete space allocated by new
        delete nets_iter->second;
    }    
	cout << endl;
	fout << endl;

    fout.close();
}

// Solver

/*
int IDAMAX(int N,vector<double> dx,int incx)
{
    int i,ix,IDAMAX;
    double dmax;

    IDAMAX = 0 ;
    if((N < 1) || (incx <= 0))return IDAMAX;
    IDAMAX = 1 ;
    if(N == 1) return IDAMAX;
    if(!(incx == 1)){
        ix = 1;
        dmax = abs(dx[1]);
        ix = ix + incx;
        for (int i = 2; i <= N; i++)
        {
            if (!(abs(dx[ix]) <= dmax))
            {
                IDAMAX=i;
                dmax=abs(dx[ix]);
            }
            ix=ix+incx;
        }
        return IDAMAX;
    }
    dmax=abs(dx[1]);
    for (int i = 2; i <= N; i++)
    {
        if(!(abs(dx[i])<=dmax))
        {
        IDAMAX=i;//this is add one
        dmax=abs(dx[i]);
        }
    }
    return IDAMAX;
}

vector<double> DSCAL(int N,double da,vector<double> dx, int incx){
    int i,m,mp1,nincx;
    if((N <= 0) || (incx <= 0)) return dx;
    if(!(incx == 1))
    {
        nincx = N*incx;
        for (int i = 0; i < nincx; i+incx)
        {
            dx.push_back(da*dx[i]);
        }
        return dx;
    }
    m = N%5;

    if(!(m == 0))
    {
        for (int i = 1; i <= m; i++)
        {
            dx.insert(dx.begin()+i,dx[i]*da);
            dx.erase(dx.begin()+1+i,dx.begin()+2+i);
        }
        if(N < 5) return dx;
    }
    mp1=m+1;

    for (int i = mp1; i <= N; i = i+5)
    {
        dx.insert(dx.begin()+i,dx[i]*da);
        dx.erase(dx.begin()+1+i,dx.begin()+2+i);

        dx.insert(dx.begin()+i+1,dx[i+1]*da);
        dx.erase(dx.begin()+1+i+1,dx.begin()+2+i+1);

        dx.insert(dx.begin()+i+2,dx[i+2]*da);
        dx.erase(dx.begin()+1+i+2,dx.begin()+2+i+2);
        
        dx.insert(dx.begin()+i+3,dx[i+3]*da);
        dx.erase(dx.begin()+1+i+3,dx.begin()+2+i+3);

        dx.insert(dx.begin()+i+4,dx[i+4]*da);
        dx.erase(dx.begin()+1+i+4,dx.begin()+2+i+4);  
    }
    return dx;
};

list<vector<double>> DAXPY(int N ,double da ,vector<double> dx,int incx,vector<double> dy,int incy)//!
{
    int i,ix,iy,m,mp1;
    list<vector<double>> doubleReplace;

    if(N <= 0) 
    {
        doubleReplace.push_front(dx);
        doubleReplace.push_back(dy); 
        return doubleReplace;
    }
    if((da < 1e-6) && (da > -1e-6))//==0 can't be use
    {
        doubleReplace.push_front(dx); 
        doubleReplace.push_back(dy);
        return doubleReplace;
    }
    if(!(incx == 1) && (incy == 1))
    {
        ix=1;
        iy=1;
        if(incx < 0) ix=(-N+1)*incx+1;
        if(incy < 0) iy=(-N+1)*incy+1;
        for (int i = 1; i <= N; i++)
        {
            dy.insert(dy.begin()+iy,dy[iy]+da*dx[ix]);//need move 1 distance
            dy.erase(dy.begin()+1+iy,dy.begin()+2+iy);
            ix=ix+incx;
            iy=iy+incy;
        }
        doubleReplace.push_front(dx); 
        doubleReplace.push_back(dy);
        return doubleReplace;
    }
    m = N%4;

    if(!(m == 0))
    {
        for (int i = 1; i <= m; i++)
        {
            dy.insert(dy.begin()+i,dy[i]+da*dx[i]);
            dy.erase(dy.begin()+1+i,dy.begin()+2+i);
        }
        if(N <= 4)
        {
        doubleReplace.push_front(dx); 
        doubleReplace.push_back(dy);
        return doubleReplace;           
        }
        
    }
    mp1=m+1;
    for (int i = mp1; i <= N; i = i+4)
    {
        dy.insert(dy.begin()+i,dy[i]+da*dx[i]);
        dy.erase(dy.begin()+1+i,dy.begin()+2+i);

        dy.insert(dy.begin()+i+1,dy[i+1]+da*dx[i+1]);
        dy.erase(dy.begin()+1+i+1,dy.begin()+2+i+1);

        dy.insert(dy.begin()+i+2,dy[i+2]+da*dx[i+2]);
        dy.erase(dy.begin()+1+i+2,dy.begin()+2+i+2);

        dy.insert(dy.begin()+i+3,dy[i+3]+da*dx[i+3]);
        dy.erase(dy.begin()+1+i+3,dy.begin()+2+i+3);
    }
    doubleReplace.push_front(dx); 
    doubleReplace.push_back(dy);
    return doubleReplace;   
    
}

vector<double> solver( map<vector<int>,double> A,vector<double> B ,int N ,int lud)
{
    int lda,info,j,k,kp1,l,nm1,kb,number;
    double t,AMD_value;
    vector<int> ADD;
    vector<int> ipvt;
    vector<double> Z,ZC,dx,dy,replace,insect;
    map<vector<int>,double> AMD;
    list<vector<double>> doubleReplace;

    nm1=N-1;//129
    Z = B;
    ipvt.insert(ipvt.begin(),0);

    if(!(lud == 0))
    {
        AMD = A;
        info=0;
        if (!(nm1 < 1)){
            for(int k = 1; k <= nm1 ; k++)// low
            {
                kp1=k+1;// low
                for(int getNumber = k; getNumber <= N; getNumber++){
                    ADD.push_back(k); 
                    ADD.push_back(getNumber);
                    dx.push_back(AMD[ADD]);
                    ADD.clear();
                }
                dx.insert(dx.begin(),0);
                l=IDAMAX(N-k+1,dx,1)+k-1 ;//this is remove +1 
                
                ipvt.push_back(l);// low id 66 is 67number
                dx.clear();

                ADD.push_back(k); 
                ADD.push_back(l);
                AMD_value = AMD[ADD];
                ADD.clear();
                bool con = true;
                if (!(round(AMD_value*1e6) == 0))
                {   
                    if (!(l == k))//need reduce one k is low too
                    {   
                        ADD.push_back(k); 
                        ADD.push_back(l);//recover the value
                        AMD_value = AMD[ADD];
                        ADD.clear();
                        t=AMD_value;

                        ADD.push_back(k); 
                        ADD.push_back(k);
                        AMD_value = AMD[ADD];
                        ADD.clear();

                        ADD.push_back(k); 
                        ADD.push_back(l);
                        AMD[ADD] = AMD_value;
                        ADD.clear();
                            
                        ADD.push_back(k); 
                        ADD.push_back(k);
                        AMD[ADD] = t;
                        ADD.clear();
                    }
                    ADD.push_back(k); 
                    ADD.push_back(k);
                    t = -1.0/AMD[ADD];
                    ADD.clear();

                    for(int getNumber = 1+k; getNumber <= N; getNumber++)
                    {
                        ADD.push_back(k); 
                        ADD.push_back(getNumber);
                        dx.push_back(AMD[ADD]);
                        ADD.clear();
                    }
                    dx.insert(dx.begin(),0);
                    replace = DSCAL(N-k,t,dx,1);
                    dx.clear();
                    for(int getNumber = 1; getNumber <= N-k; getNumber++)
                    {
                        ADD.push_back(k); 
                        ADD.push_back(getNumber+k);
                        AMD[ADD] = replace[getNumber];//need -1
                        ADD.clear();
                    }


                    for(int j = kp1 ; j <= N;j++){
                        ADD.push_back(j); 
                        ADD.push_back(l);
                        t = AMD[ADD];
                        ADD.clear();
                        if(!(l == k)) 
                        {
                            ADD.push_back(j); 
                            ADD.push_back(k);
                            AMD_value = AMD[ADD];
                            ADD.clear();

                            ADD.push_back(j); 
                            ADD.push_back(l);
                            AMD[ADD] = AMD_value;
                            ADD.clear();
                            
                            ADD.push_back(j); 
                            ADD.push_back(k);
                            AMD[ADD] = t;
                            ADD.clear();
                        }
                       
                    for(int getNumber = 1+k; getNumber <= N; getNumber++)
                    {
                        ADD.push_back(k); 
                        ADD.push_back(getNumber);
                        dx.push_back(AMD[ADD]);
                        ADD.clear();
                    }
                    dx.insert(dx.begin(),0);
                    for(int getNumber = k+1; getNumber <= N; getNumber++)
                    {   //input is same but output be delect
                        ADD.push_back(j); 
                        ADD.push_back(getNumber);
                        dy.push_back(AMD[ADD]);
                        ADD.clear();
                    }
                    dy.insert(dy.begin(),0);
                    doubleReplace = DAXPY(N-k,t,dx,1,dy,1);//size is tooooooo low
                    dx.clear();
                    dy.clear();
                    for(int getNumber = 1; getNumber <= N-k; getNumber++)
                    {
                        replace = doubleReplace.front();
                        ADD.push_back(k); 
                        ADD.push_back(getNumber+k);
                        AMD[ADD] = replace[getNumber];
                        ADD.clear();
                    }
                    for(int getNumber = 1; getNumber <= N-k; getNumber++)
                    {
                        replace = doubleReplace.back();
                        ADD.push_back(j); 
                        ADD.push_back(getNumber+k);
                        AMD[ADD] = replace[getNumber];
                        ADD.clear();
                    }

                    }
                    bool con = false;
                }
                if(con)info=k;  //low
            }
        }
        // ipvt.insert(ipvt.begin()+N,N);
        // ipvt.erase(ipvt.begin()+N+1,ipvt.begin()+N+2);
        ipvt.push_back(N);

        ADD.push_back(N); 
        ADD.push_back(N);
        AMD_value = AMD[ADD];
        ADD.clear();
        if ((AMD_value < 1e-6) && (AMD_value > -1e-6)) info=N;
        if (info == 0) cout << "Division by zero in SOLVER!" << endl;
    }
    if (!(nm1 < 1)) 
    {
        for (int k = 1; k <= nm1; k++)
        {
           l=ipvt[k];
           t=Z[l];
            if (!(l == k))
            {
                Z.insert(Z.begin()+l,Z[k]);
                Z.erase(Z.begin()+1+l,Z.begin()+2+l);

                Z.insert(Z.begin()+k,t);
                Z.erase(Z.begin()+1+k,Z.begin()+2+k);
            }
                    
            for(int getNumber = 1+k; getNumber <= N; getNumber++)
            {
                ADD.push_back(k); 
                ADD.push_back(getNumber);
                dx.push_back(AMD[ADD]);
                ADD.clear();
            }
            dx.insert(dx.begin(),0);
            for(int getNumber = 1+k; getNumber <= N ; getNumber++)
            {
                ZC.push_back(Z[getNumber]);
            }
            ZC.insert(ZC.begin(),0);
            doubleReplace = DAXPY(N-k,t,dx,1,ZC,1);
            dx.clear();
            for(int getNumber = 1; getNumber <= N-k; getNumber++)
            {
                replace = doubleReplace.front();//can be better
                ADD.push_back(k); 
                ADD.push_back(getNumber+k);
                AMD[ADD] = replace[getNumber];
                ADD.clear();
            }
            number = k;
            ZC.clear();
            
            for(int getNumber = 1; getNumber <= number; getNumber++)
            {
                insect.push_back(Z[getNumber]); 
            }  
            insect.insert(insect.begin(),0);
            replace = doubleReplace.back(); 
            replace.erase(replace.begin(),replace.begin()+1);
            insect.insert(insect.end(),replace.begin(),replace.end());
            Z.swap(insect);
            insect.clear();   
        }
    }
    for(int kb = 1;kb <= N;kb++){
        k=N+1-kb;
        ADD.push_back(k); 
        ADD.push_back(k);
        AMD_value = AMD[ADD];
        ADD.clear();

        Z.insert(Z.begin()+k,Z[k]/AMD_value);//this value is null 
        Z.erase(Z.begin()+1+k,Z.begin()+2+k);

        t=-Z[k];
        for(int getNumber = 1; getNumber <= N; getNumber++)
        {
            ADD.push_back(k); 
            ADD.push_back(getNumber);
            dx.push_back(AMD[ADD]);
            ADD.clear();
        }
        dx.insert(dx.begin(),0);
        doubleReplace = DAXPY(k-1,t,dx,1,Z,1);
        dx.clear();

        replace = doubleReplace.front();
        for(int getNumber = 1; getNumber <= N; getNumber++)
        {
            ADD.push_back(k);
            ADD.push_back(getNumber);
            AMD[ADD] = replace[getNumber];
            ADD.clear();
        }
        Z = doubleReplace.back();
    }
    Z.erase(Z.begin(),Z.begin()+1);
    return Z;
}
*/
