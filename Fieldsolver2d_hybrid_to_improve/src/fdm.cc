#ifndef COMMON_IS_INCLUDED
#include "common.h"
#endif
#include "cholmod.h"

using namespace std;

// Define spacing
const double dz = 0.2;             // uniform spacing
const double dz_1_10 = dz / 10;    // spacing of out boundary
const double dz_1_100 = dz / 100;  // minimum spacing of net boundary

const double dx = 0.2;             // uniform spacing
const double dx_1_10 = dx / 10;    // spacing of out boundary
const double dx_1_100 = dx / 100;  // minimum spacing of net boundary

const double inc_multiple = 1.5;   // multiples of spacing increment

/*
 * interplate between exclusive intervals.
 * start is the net boundary, if start < end, coordinate increases; otherwise, decreases.
 */

void print_dense_matrix(cholmod_dense *A)
{
    int nrow = A->nrow;
    int ncol = A->ncol;
    for (int i=0; i<nrow; i++)
        for (int j=0; j<ncol; j++) {
            printf("%6.3f  ",((double *)A->x)[i*ncol+j]);
            if (j == ncol-1)
                printf("\n");
        }
}

void interpolate_exclusive(double start, double end, double dmin, double multiple, vector<double> &v)
{
    vector<double> v_tmp;
    double z;
    double d = dmin;
    // from start to end
    if (start + eps < end) {
        z = start + d;
        //while (z < end) {
        //while (z + eps < end) {
        while ((z + eps < end) && (end - z > dmin + eps)) {
            v.push_back(z);
            d *= multiple;
            z += d;
        }
    // from end to start
    } else {
        z = start - d;
        //while (z > end) {
        //while (z > end + eps) {
        while ((z > end + eps) && ( z - end > dmin + eps)) {
            v_tmp.push_back(z);
            d *= multiple;
            z -= d;
        }
        reverse(v_tmp.begin(), v_tmp.end());
        v.insert(v.end(), v_tmp.begin(), v_tmp.end());
    }
}

// Define class for sentry
class Sentry
{
    public:
        double lo_val();
        double hi_val();
        void inc();
        Sentry(vector<double> v, double lo, double hi);
    private:
        unsigned int lo_idx;
        unsigned int hi_idx;
        vector<double> nets_coord;
        double boundary_lo_coord;
        double boundary_hi_coord;
        double val(unsigned int idx);
};

Sentry::Sentry(vector<double> v, double lo, double hi)
{
    lo_idx = 0;
    hi_idx = 1;
    nets_coord = v;
    boundary_lo_coord = lo;
    boundary_hi_coord = hi;
}

void Sentry::inc()
{
    lo_idx++;
    hi_idx++;
}

double Sentry::val(unsigned int idx)
{
    double v;
    if (idx == 0) {
        v = boundary_lo_coord;
    } else if (idx > nets_coord.size()) {
        v = boundary_hi_coord;
    } else {
        v = nets_coord[idx-1];
    }
    return v;
}

double Sentry::lo_val() {
    return val(lo_idx);
}

double Sentry::hi_val() {
    return val(hi_idx);
}

//-------------------------------------------------------------------
//   generate grids for one coordinate
//   input: nets, lo, hi
//   output: grids
//
//   nets: the coordinates of all nets along one direction
//   lo: low boundary coordinate along one direction
//   hi: high boundary coordinate along one direction
//
//   don't interpolate near the out boundary
//-------------------------------------------------------------------
void generate_grids_2(vector<double> &grids, vector<double> nets, double lo, double hi)
{
    double t;
    Sentry sentry(nets, lo, hi);
    t = lo;
//    while (t < hi) {
    while (t + eps < hi) {
        grids.push_back(t);

//        if (sentry.hi_val() - sentry.lo_val() < 2*dx) {
        if (sentry.hi_val() - sentry.lo_val() + eps < 2*dx) {
            t = 0.5*(sentry.hi_val() + sentry.lo_val());
            interpolate_exclusive (sentry.lo_val(), t, dx_1_100, inc_multiple, grids);
            grids.push_back(t);
            interpolate_exclusive (sentry.hi_val(), t, dx_1_100, inc_multiple, grids);
            t = sentry.hi_val();
            sentry.inc();
//        } else if ((t == sentry.lo_val()) && (sentry.lo_val() != lo) ){
        } else if ((abs(t - sentry.lo_val()) < eps) && (abs(sentry.lo_val() - lo) > eps) ){
            t = sentry.lo_val() + dx;
            interpolate_exclusive (sentry.lo_val(), t, dx_1_100, inc_multiple, grids);
//        } else if ((sentry.hi_val() - t < 2*dx) && (sentry.hi_val() != hi)) {
        } else if ((sentry.hi_val() - t + eps < 2*dx) && (abs(sentry.hi_val() - hi) > eps) ) {
            t = sentry.hi_val() - dx;
            grids.push_back(t);
            interpolate_exclusive (sentry.hi_val(), t, dx_1_100, inc_multiple, grids);
            t = sentry.hi_val();
            sentry.inc();
        } else {
            t += dx;
        }
    }
    grids.push_back(hi);
}

//bool on_boundary(double x, double z, Boundary b, map<string, Net*> nets, int &u) {
bool on_boundary(double x, double z, Boundary b, map<int, Net*> nets, int &u) {
    unsigned int i;
    bool on;
    //map<string, Net*>::iterator nets_iter;
    map<int, Net*>::iterator nets_iter;
    on = false;
    u = 0;
//  if (x == b.x0 || x == b.x1 || z == b.z0 || z == b.z1) {
    if ((fabs(x-b.x0) < eps) || (fabs(x - b.x1) < eps) || (fabs(z - b.z0) < eps) || (fabs(z - b.z1) < eps)) {
        on = true;
    } else {
        for (nets_iter = nets.begin(); nets_iter != nets.end(); ++nets_iter) {
            vector<Subnet> n = nets_iter->second->subnets;
            for (i = 0; i < n.size(); i++) {
//              if (x >= n[i].x0 && x <= n[i].x1 && z >= n[i].z0 && z <= n[i].z1) {
            if ((x > n[i].x0 + eps || fabs(x - n[i].x0) < eps) &&
                        (x < n[i].x1 - eps || fabs(x - n[i].x1) < eps) &&
                        (z > n[i].z0 + eps || fabs(z - n[i].z0) < eps) &&
                        (z < n[i].z1 - eps || fabs(z - n[i].z1) < eps)) {
                    on = true;
                    //u = nets_iter->first == "net0" ? 1 : 0;
                    u = nets_iter->first == 0 ? 1 : 0;
                    break;
                }
            }
        }
    }
    return on;
}

void derive_cofficient(double &aa, double &al, double &ar, double &au, double &ad,
        int i, int j, vector<double> x_coords, vector<double> z_coords) {
    double dl, dr, du, dd;
    double xl, xr, zu, zd;
    double xx, zz;
    xx = x_coords[i];
    zz = z_coords[j];
    xl = x_coords[i-1];
    xr = x_coords[i+1];
    zu = z_coords[j+1];
    zd = z_coords[j-1];
    dl = xx - xl;
    dr = xr - xx;
    du = zu - zz;
    dd = zz - zd;
    aa = (1/(dl*dr) + 1/(du*dd)) * (dl + dr) * (du + dd);
    al = -1 * (du + dd) / dl;
    ar = -1 * (du + dd) / dr;
    au = -1 * (dl + dr) / du;
    ad = -1 * (dl + dr) / dd;
}

// rectangle method for integration and centered difference for du
double compute_charge_3(int lo_x, int lo_z, int hi_x, int hi_z, int Nz,
        vector<double> x_coords, vector<double> z_coords, double *ux, double epislon)
{
    int i,j;
    double du; // u2 - u1
    double dl; // distance along the voltage direction
    double ds; // distance othorgonal to the voltage direction
    double sum = 0;

    //left boundary, E is towards left, from bottom to up
    i = lo_x;
    //cout << "lo_x:" << endl;
    for (j = lo_z; j < hi_z; j++) {
        du = fabs(ux[(i+1)*Nz+j] - ux[(i-1)*Nz+j]);
        dl = fabs(x_coords[i+1] - x_coords[i-1]);
        ds = z_coords[j+1] - z_coords[j];
        sum += (du * ds) / dl ;
    }

    //up boundary, E is towards up, from left to right
    j = hi_z;
    //cout << "hi_z:" << endl;
    for (i = lo_x; i < hi_x; i++) {
        du = fabs(ux[i*Nz+(j+1)] - ux[i*Nz+(j-1)]);
        dl = fabs(z_coords[j+1] - z_coords[j-1]);
        ds = x_coords[i+1] - x_coords[i];
        sum += (du * ds) / dl ;
    }
        //right boundary, E is towards right, from up to bottom
    i = hi_x;
    //cout << "hi_x:" << endl;
    for (j = hi_z; j > lo_z; j--) {
        du = fabs(ux[(i+1)*Nz+j] - ux[(i-1)*Nz+j]);
        dl = fabs(x_coords[i+1] - x_coords[i-1]);
        ds = z_coords[j] - z_coords[j-1];
        sum += (du * ds) / dl ;
    }

    //bottom boundary, E is towards down, from right to left
    j = lo_z;
    //cout << "lo_z:" << endl;
    for (i = hi_x; i > lo_x; i--) {
        du = fabs(ux[i*Nz+(j+1)] - ux[i*Nz+(j-1)]);
        dl = fabs(z_coords[j+1] - z_coords[j-1]);
        ds = x_coords[i] - x_coords[i-1];
        sum += (du * ds) / dl ;
    }

    sum *= (epislon * e0);
    return sum;
}

//int cholmod_solver(size_t Nx, size_t Nz, size_t NN, size_t Nb, size_t Nb_net0, 
//		Boundary boundary, vector<double> x_coords, vector<double> z_coords,
//		map<string, Net*> nets, map<string, Net*>::iterator nets_iter, double dielectric) {
int cholmod_solver(size_t Nx, size_t Nz, size_t NN, size_t Nb, size_t Nb_net0, 
		Boundary boundary, vector<double> x_coords, vector<double> z_coords,
		map<int, Net*> nets, map<int, Net*>::iterator nets_iter, double dielectric) {
	// dimension of A
    size_t A_nrows = NN;
    size_t A_ncols = NN;
    size_t A_nnz = Nb + (NN - Nb)*5;

	// dimension of b
    size_t b_nrows = NN;
    size_t b_ncols = 1;
    size_t b_nnz = Nb_net0;

	// A common struct that cholmod always needs
    cholmod_common c;

	// start cholmod
    cholmod_start (&c);
	//c.print=5 will print out the detailed matrix content
	//c.print = 5;
	//c.supernodal = CHOLMOD_SIMPLICIAL;
	
	// a triplet that is very helpful to build a sparse matrix
    cholmod_triplet *TA;
    cholmod_triplet *Tb;

    // create sparse matrix
    cholmod_sparse *A;
    cholmod_sparse *b;

	// allocate a triplet TA, that we will use to generate A matrix
    // cholmod_allocate_triplet(nrow, ncol, nzmax, stype, xtype)
    // if assign nzmax with A_nows*A_ncols, it will report the problem is too large, so A_nrows*5 is a proper value
    // but here is better since we directly give it the correct nz numbers
    TA = cholmod_allocate_triplet(A_nrows, A_ncols, A_nnz, 0, CHOLMOD_REAL, &c);

    // allocate a triplet Tb, that we will use to generate b matrix
    Tb = cholmod_allocate_triplet(b_nrows, b_ncols, b_nnz, 0, CHOLMOD_REAL, &c);

	int *TAi = (int *) TA->i;
    int *TAj = (int *) TA->j;
    double *TAx = (double *) TA->x;

    int *Tbi = (int *) Tb->i;
    int *Tbj = (int *) Tb->j;
    double *Tbx = (double *) Tb->x;

    int ka = 0;
    int kb = 0;
    int row;
    double aa, al, ar, au, ad;

	int u;
	unsigned int i, j;
	for (i = 0; i < Nx; i++) {
        for (j = 0; j < Nz; j++) {
            row = i * Nz + j;
            // grid boundary point has only one cofficient with value 1
            if (on_boundary(x_coords[i], z_coords[j], boundary, nets, u)) {
                TAi[ka] = row; TAj[ka] = row; TAx[ka] = 1; ka++; //A_nnz++;
                if (u == 1) {
                    Tbi[kb] = row; Tbj[kb] = 0; Tbx[kb] = 1; kb++; //b_nnz++;
                }
            } else { // other grid points have five cofficients
                derive_cofficient(aa, al, ar, au, ad, i, j, x_coords, z_coords);
                // al
                TAi[ka] = row; TAj[ka] = (i-1)*Nz+j; TAx[ka] = al; ka++; //A_nnz++;
                // ad
                TAi[ka] = row; TAj[ka] = i*Nz+(j-1); TAx[ka] = ad; ka++; //A_nnz++;
                // aa
                TAi[ka] = row; TAj[ka] = row;        TAx[ka] = aa; ka++; //A_nnz++;
                // au
                TAi[ka] = row; TAj[ka] = i*Nz+(j+1); TAx[ka] = au; ka++; //A_nnz++;
                // ar
                TAi[ka] = row; TAj[ka] = (i+1)*Nz+j; TAx[ka] = ar; ka++; //A_nnz++;
            }
        }
    }

	// T must be told what nnz is
    TA->nnz = A_nnz;
    Tb->nnz = b_nnz;

	// create sparse matrix from triplet
    A = cholmod_triplet_to_sparse(TA, A_nnz, &c);
    b = cholmod_triplet_to_sparse(Tb, b_nnz, &c);

    // print the matrix
    cholmod_print_sparse(A, "A", &c);
    cholmod_print_sparse(b, "b", &c);

    // create dense matrics
    cholmod_dense *x;
    // the cholesky factor
    cholmod_factor *L;

	//---- solve A'Ax=A'b -----
    // compute At_A = A'A
    cholmod_sparse *At = cholmod_transpose(A, 2, &c);
    cholmod_sparse *At_A = cholmod_ssmult(At, A, 1, 1 ,1, &c);
    cholmod_print_sparse(At_A, "At_A", &c);

    printf("cholmod_check_sparse: %d\n", cholmod_check_sparse(At_A, &c));

	// compute At_b = A'b
    double alpha[2] = {1,1};
    double beta[2] = {1,1};
    cholmod_dense *At_b = cholmod_zeros(A->nrow, 1, A->xtype, &c);
    // 0: compute A*b
    // 1: compute A'*b
    cholmod_dense *bb = cholmod_sparse_to_dense(b, &c);
    cholmod_sdmult(A, 1, alpha, beta, bb, At_b, &c);

    // compute At_A \ At_b
    L = cholmod_analyze(At_A, &c);
    cholmod_factorize(At_A, L, &c);
    x = cholmod_solve(CHOLMOD_A, L, At_b, &c);

	// print dense matrix
	//cholmod_print_dense(x,"x",&c);
	//print_dense_matrix(x);
	

	// print residual error
	// Here we can also determine if cholmod solver is successful by checking if the result is nan
    double one [2] = {1,0}, m1 [2] = {-1,0} ;
    cholmod_dense *r = cholmod_copy_dense (bb, &c) ;
    cholmod_sdmult (A, 0, m1, one, x, r, &c) ;
    double residual_error = cholmod_norm_dense (r, 0, &c);

    // free memory
	cholmod_free_dense (&r, &c);
    cholmod_free_factor (&L, &c);
    cholmod_free_sparse (&At, &c);
    cholmod_free_sparse (&At_A, &c);
    cholmod_free_dense (&At_b, &c);
    cholmod_free_triplet (&TA, &c);
    cholmod_free_sparse (&A, &c);
    cholmod_free_triplet (&Tb, &c);
    cholmod_free_sparse (&b, &c);
    cholmod_free_dense (&bb, &c);

	int is_ok = 1;
	// check nan
    if (isnan(residual_error)) {
        printf ("WARNING: nan appeared and cholmod solver failed!\n") ;
		is_ok = 0;
    } else {
        printf ("norm(b-Ax) %8.1e\n", residual_error) ;
		// -----------------------compute charge----------------------//
   		double *ux = (double *) x->x;
   		int gau_lo_x, gau_lo_z, gau_hi_x, gau_hi_z; // gauss boundary index
   		double charge;
   		cout << "------------------ net cap -----------------" << endl;
   		for (nets_iter = nets.begin(); nets_iter != nets.end(); ++nets_iter) {
   		    gau_lo_x = nets_iter->second->lo_x - 2;
   		    gau_lo_z = nets_iter->second->lo_z - 2;
   		    gau_hi_x = nets_iter->second->hi_x + 2;
   		    gau_hi_z = nets_iter->second->hi_z + 2;

   		    charge = compute_charge_3(gau_lo_x, gau_lo_z, gau_hi_x, gau_hi_z, Nz,
   		            x_coords, z_coords, ux, dielectric);
   		    cout << "net" << nets_iter->first << " : " << fixed << setprecision(6) << charge << " fF" << endl;
			// restore cout precision
		    cout <<resetiosflags(ios::fixed) <<setprecision(6); 

   		    nets_iter->second->set_cap(charge);
   		}
   		// -----------------------------------------------------------//
    }

 	// end cholmod
	cholmod_finish (&c);
	cholmod_free_dense (&x, &c);

	return is_ok;
}

//int fdm(map<string, Net*> nets, map<string, Net*>::iterator nets_iter, Boundary boundary, 
//		vector<double> nets_x, vector<double> nets_z, double dielectric) {
int fdm(map<int, Net*> nets, map<int, Net*>::iterator nets_iter, Boundary boundary, 
		vector<double> nets_x, vector<double> nets_z, double dielectric) {

	cout << endl;
    cout << "============================================" << endl;
    cout << ">                fdm is used               <" << endl;
    cout << "============================================" << endl;
    cout << endl;

	// Sort&uniq nets_x & nets_z for later grid generation
    sort_unique(nets_x);
    sort_unique(nets_z);

	// Show nets_x & nets_z
	/*
    cout << "---------------- nets_x_coord --------------" << endl;
    print_vector(nets_x);
    cout << "---------------- nets_z_coord --------------" << endl;
    print_vector(nets_z);
	*/

	// Generate z_coords for all grids
    vector<double> z_coords;
    generate_grids_2(z_coords, nets_z, boundary.z0, boundary.z1);

	/*
    // Show z_coords
    cout << "------------- generated z_coords -----------" << endl;
    cout << "z_coords size: " << z_coords.size() << endl;
    print_vector(z_coords, 1);
	*/

    // Generate x_coords for all grids
    vector<double> x_coords;
    generate_grids_2(x_coords, nets_x, boundary.x0, boundary.x1);

	/*
    // Show x_coords
    cout << "------------- generated x_coords -----------" << endl;
    cout << "x_coords size: " << x_coords.size() << endl;
    print_vector(x_coords, 1);
	*/

	// Show grids statistics
    // NN is necessary for cholmod solver to define nrow or ncol
	// Nb, Nb_net0 are used to define nzmax. Actually here they are also used to specify nnz.
	// Although nnz can be worked out when assigning non-element values for triplet matrix,
	// but here they are worked out in prior, which helps to debug for printing information.
	// The experiment shows that working them out in prior doesn't consume too much extra time
    cout << "-------------- grids statistics ------------" << endl;
    unsigned int i, j;
    // use long type to accommodate big numbers
    unsigned long int Nx, Nz, NN;
    int u; // boundary conditions
    Nx = x_coords.size();
    Nz = z_coords.size();
    NN = Nx * Nz; // # of all grid points

	unsigned long int Nb; // # of grid points on boundaries
    unsigned long int Nb_net0; // # of grid points on net0

    vector<double>::iterator iter;
    unsigned int net_lo_x, net_lo_z, net_hi_x, net_hi_z; // nets boundary index

    Nb_net0 = 0;
    // first count points on out boundary
    // each virtex is recounted once so need to substract 4
    Nb = (Nx + Nz) * 2 - 4;

	for (nets_iter = nets.begin(); nets_iter != nets.end(); ++nets_iter) {
        //net_lo_x
//      iter = find(x_coords.begin(), x_coords.end(), nets_iter->second->x0);
        iter = find_if(x_coords.begin(), x_coords.end(), dbl_cmp(nets_iter->second->x0,eps));
        net_lo_x = distance(x_coords.begin(), iter);
        //net_hi_x
//      iter = find(x_coords.begin(), x_coords.end(), nets_iter->second->x1);
        iter = find_if(x_coords.begin(), x_coords.end(), dbl_cmp(nets_iter->second->x1,eps));
        net_hi_x = distance(x_coords.begin(), iter);
        //net_lo_z
//      iter = find(z_coords.begin(), z_coords.end(), nets_iter->second->z0);
        iter = find_if(z_coords.begin(), z_coords.end(), dbl_cmp(nets_iter->second->z0,eps));
        net_lo_z = distance(z_coords.begin(), iter);
        //net_hi_z
//      iter = find(z_coords.begin(), z_coords.end(), nets_iter->second->z1);
        iter = find_if(z_coords.begin(), z_coords.end(), dbl_cmp(nets_iter->second->z1,eps));
        net_hi_z = distance(z_coords.begin(), iter);

        nets_iter->second->set_index(net_lo_x, net_lo_z, net_hi_x, net_hi_z);

		for (i = net_lo_x; i <= net_hi_x; i++) {
            for (j = net_lo_z; j <= net_hi_z; j++) {
                if (on_boundary(x_coords[i], z_coords[j], boundary, nets, u)) {
                    Nb++;
                    if (u == 1) {
                        Nb_net0++;
                    }
                }
            }
        }
	}

	cout << "# of all grid points = " << NN <<endl;
    cout << "# of boundary points = " << Nb <<endl;
    cout << "# of boundary points on net0 = " << Nb_net0 <<endl;
    cout << "# of all elements in matrix = " << NN*NN << endl;
    cout << "# of all non-zero elements in matrix = " << Nb + (NN-Nb)*5 << endl;

	int is_ok = 1;
	is_ok = cholmod_solver(Nx, Nz, NN, Nb, Nb_net0, 
			boundary, x_coords, z_coords, nets, nets_iter, dielectric); 
	if (is_ok) {
        cout << endl;
        cout << "============================================" << endl;
        cout << "*               fdm succeeded!             *" << endl;
        cout << "============================================" << endl;
        cout << endl;
    } else {
        cout << endl;
        cout << "============================================" << endl;
        cout << "#                 fdm failed!              #" << endl;
        cout << "============================================" << endl;
        cout << endl;
    }

	return is_ok;
}
