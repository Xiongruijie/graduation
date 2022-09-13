#ifndef COMMON_IS_INCLUDED
#include "common.h"
#endif
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

const double pi = 4*atan(1);

void CPF(double xi, double eta, double xk, double yk, double nkx, double nky, double L, double& PF1, double& PF2)
{
    double A,B,E,D,BA,EA;
    A = pow(L,2);
    B = 2*L*(-nky*(xk-xi)+nkx*(yk-eta));
    E = pow((xk-xi),2)+pow((yk-eta),2);
    D = sqrt(abs(4*A*E-pow(B,2)));
    BA=B/A;
    EA=E/A;
    if (dbl_equ(D, 0)) {
        PF1=0.5*L*(log(L)+(1+0.5*BA)*log(abs(1+0.5*BA))-0.5*BA*log(abs(0.5*BA))-1);
        PF2=0;
    } else {
        PF1=0.25*L*(2*(log(L)-1)-0.5*BA*log(abs(EA))+(1+0.5*BA)*log(abs(1+BA+EA))+(D/A)*(atan((2*A+B)/D)-atan(B/D)));
        PF2=L*(nkx*(xk-xi)+nky*(yk-eta))/D*(atan((2*A+B)/D)-atan(B/D));
    }
}

int generate_elements(int flag, double common, double start, double end, double el, vector<double>&xk, vector<double>&zk,
        vector<double>&xm, vector<double>&zm, vector<double>&nx, vector<double>&nz, vector<double>&lk)
{
    // recording number of generated elements
    int num=0;
    double x, z;
    switch(flag) {
        case 0: // left side
            x = common;
            z = start;
            while (z + eps < end) {
                num++;
                xk.push_back(x);
                zk.push_back(z);
                nx.push_back(1);
                nz.push_back(0);
                xm.push_back(x);
                if (z + el > end + eps) {
                    zm.push_back((end+z)/2);
                    lk.push_back(end-z);
                } else {
                    zm.push_back(z+0.5*el);
                    lk.push_back(el);
                }
                z += el;
            }
            break;
        case 1: // up side
            z = common;
            x = start;
            while (x + eps < end) {
                num++;
                xk.push_back(x);
                zk.push_back(z);
                nx.push_back(0);
                nz.push_back(-1);
                zm.push_back(z);
                if (x + el > end + eps) {
                    xm.push_back((end+x)/2);
                    lk.push_back(end-x);
                } else {
                    xm.push_back(x+0.5*el);
                    lk.push_back(el);
                }
                x += el;
            }
            break;
        case 2: // right side
            x = common;
            z = start;
            while (z > end + eps) {
                num++;
                xk.push_back(x);
                zk.push_back(z);
                nx.push_back(-1);
                nz.push_back(0);
                xm.push_back(x);
                if (z - el + eps < end) {
                    zm.push_back((z+end)/2);
                    lk.push_back(z-end);
                } else {
                    zm.push_back(z-0.5*el);
                    lk.push_back(el);
                }
                z -= el;
            }
            break;
        case 3: // bottom side
            x = start;
            z = common;
            while (x > end + eps) {
                num++;
                xk.push_back(x);
                zk.push_back(z);
                nx.push_back(0);
                nz.push_back(1);
                zm.push_back(z);
                if (x - el + eps < end) {
                    xm.push_back((x+end)/2);
                    lk.push_back(x-end);
                } else {
                    xm.push_back(x-0.5*el);
                    lk.push_back(el);
                }
                x -= el;
            }
            break;
        default:
            cout << "wrong flag: " << flag << endl;
            exit(1);
    }
    return num;
}

int generate_net_elements(double el, vector<double>&xk, vector<double>&zk, vector<double>&xm, vector<double>&zm,
        vector<double>&nx, vector<double>&nz, vector<double>&lk, vector<double>&phi, map<int, Net*> nets)
{
    // recording number of generated elements
    int num=0;

    map<int, Net*>::iterator nets_iter;
    for (nets_iter = nets.begin(); nets_iter != nets.end(); ++nets_iter) {
        int subnum = 0;
        for (list<Border>::iterator iter = nets_iter->second->borders.begin();
                iter != nets_iter->second->borders.end();
                    ++iter) {
            subnum += generate_elements(iter->flag, iter->common, iter->start, iter->end,
                    el, xk, zk, xm, zm, nx, nz, lk);
        }
        // record staring position and length in the element array for each net
        nets_iter->second->set_seqlen(phi.size(), subnum);
//        if (nets_iter->first == "net0") {
        if (nets_iter->first == 0) {
            phi.insert(phi.end(), subnum, 1);
        } else {
            phi.insert(phi.end(), subnum, 0);
        }
        num += subnum;
    }
    return num;
}


int generate_boundary_elements(double el, vector<double>&xk, vector<double>&zk, vector<double>&xm, vector<double>&zm,
        vector<double>&nx, vector<double>&nz, vector<double>&lk, vector<double>&phi, Boundary boundary)
{
    // recording number of generated elements
    int num=0;
    // left side
    num += generate_elements(0, boundary.x0, boundary.z0, boundary.z1, el, xk, zk, xm, zm, nx, nz, lk);
    // up side
    num += generate_elements(1, boundary.z1, boundary.x0, boundary.x1, el, xk, zk, xm, zm, nx, nz, lk);
    // right side
    num += generate_elements(2, boundary.x1, boundary.z1, boundary.z0, el, xk, zk, xm, zm, nx, nz, lk);
    // bottom side
    num += generate_elements(3, boundary.z0, boundary.x1, boundary.x0, el, xk, zk, xm, zm, nx, nz, lk);
    // phi value for out boundary is 0
    phi.insert(phi.end(), num, 0);
    return num;
}

int eigen_solver(map<int, Net*>nets, map<int, Net*>::iterator nets_iter, int num_tot,
    vector<double> xk, vector<double> zk, vector<double> xm, vector<double> zm,
	vector<double> nx, vector<double> nz, vector<double> lk, vector<double> phi, 
	double dielectric)
{
	// generate A & B
    map<vector<int>,double> A; // A is 2-dimensional array
    vector<double> B, Z;
    vector<int> ADD;
    double PF1, PF2;
    double F1, F2;
    double del;
    double bsum;

	// eigen
	Matrix<double, Dynamic, Dynamic> matrix_NNd(num_tot,num_tot);
    VectorXd matrix_1Nd(num_tot);
    VectorXd matrix_N1d;

	for(int m=0; m<num_tot; m++) {
        bsum = 0;
        for(int k=0; k<num_tot; k++)
        {
            CPF(xm[m],zm[m],xk[k],zk[k],nx[k],nz[k],lk[k], PF1, PF2);
            F1 = PF1/pi;
            F2 = PF2/pi;
            if (k==m) {
                del = 1.0;
            } else {
                del = 0.0;
            }
            bsum += phi[k]*(-F2 + 0.5*del);
            matrix_NNd(m,k) = -F1;//this change   matrix_NNd(k,m)
        }
        matrix_1Nd(m) = bsum;
    }

	matrix_N1d = matrix_NNd.fullPivHouseholderQr().solve(matrix_1Nd);

	for (nets_iter = nets.begin(); nets_iter != nets.end(); ++nets_iter) {
        double charge=0;
        int m = nets_iter->second->seq_start;
        int k = m + nets_iter->second->seq_length;
        while (m < k) {
            charge += matrix_N1d(m)*lk[m];
            m++;
        }
        charge *= (e0 * dielectric);
		charge = fabs(charge);

		cout << "net" << nets_iter->first << " : " << fixed << setprecision(6) << charge << " fF" << endl;
        // restore cout precision
        cout <<resetiosflags(ios::fixed) <<setprecision(6);
//        nets_iter->second->set_cap(fabs(charge));
        nets_iter->second->set_cap(charge);
    }
	return 1;
}

int bem(map<int, Net*>nets, map<int, Net*>::iterator nets_iter, Boundary boundary, double dielectric, double shrink_val)
{
	cout << endl;
    cout << "============================================" << endl;
    cout << ">                bem is used               <" << endl;
    cout << "============================================" << endl;
    cout << endl;

	boundary.shrink(shrink_val);

	// Show gauss curve
	print_gauss(nets);

	// ----------------- Generate elements --------------------
    // xk,zk: starting coordinates of each element
    // xm,zm: mid point coordinates of each element
    // nx,nz: normal coorindates of each element
    // lk: length of each element
    // phi: value of phi of each element
    //-----------------------------------------------------
    // nx(k) = (z(k+1) - z(k))/l(k)
    // nz(k) = (x(k) - x(k+1))/l(k)
    //-----------------------------------------------------

    vector<double> xk,zk,xm,zm,nx,nz,lk,phi;
    int num_bnd = 0;
    int num_net = 0;
    int num_tot; // number of total elements
    // Generate elements for out boundary by giving a larger length
    num_bnd = generate_boundary_elements(1,xk,zk,xm,zm,nx,nz,lk,phi,boundary);
    // Generate elements for net boundary by giving a small length
    num_net = generate_net_elements(0.01,xk,zk,xm,zm,nx,nz,lk,phi,nets);

    num_tot = num_bnd + num_net;

    cout << "-------------- Generate elements -----------" << endl;
    cout << "# of elements on out boundary: "<< num_bnd << endl;
    cout << "# of elements on net boundary: "<< num_net << endl;
    cout << "# of elements in a whole: "<< num_tot << endl;

    for (nets_iter = nets.begin(); nets_iter != nets.end(); ++nets_iter) {
        cout << "Postion and length of net" << nets_iter->first << ": "
            << nets_iter->second->seq_start << ", "
            << nets_iter->second->seq_length << endl;
    }

	int is_ok;
	is_ok = eigen_solver(nets, nets_iter, num_tot, xk, zk, xm, zm, nx, nz, lk, phi, dielectric);

	cout << endl;
	cout << "============================================" << endl;
	cout << "*               bem succeeded!             *" << endl;
	cout << "============================================" << endl;

	return is_ok;
}
