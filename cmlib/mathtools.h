#ifndef __MATHTOOLS_H__       // __MATHTOOLS_H__
#define __MATHTOOLS_H__

#include "types.h"
#include "io.h"

int ludcmp(fp_t **a,int ll,int ul,int *indx,fp_t &d,io_class &io);
void lubksb(fp_t **a,int ll,int ul,int *indx,fp_t *b);
void Shipley_Inversion(fp_t **MI,int ll,int ul);
int svdcmp(fp_t **a, int m, int n, fp_t *w, fp_t **v);
int svd_analyze(fp_t * M, int start, int N);
fp_t ** solve_svd(fp_t ** Matrix, fp_t ** rhs, int start, int N, int N_equations);
void invert(fp_t **M,fp_t **MI,int m,int n);

int multiply_4x4(fp_t ** A, fp_t ** B, fp_t ** result);
int multiply_4x4(fp_t ** A, fp_t *B, fp_t * result);
fp_t ** multiply_square(fp_t ** A, fp_t ** B, int dim);
fp_t * multiply_vector(fp_t ** A, fp_t * B, int dim);
fp_t * multiply_vector(fp_t ** A, fp_t * B, int N_rows, int N_columns); // same as above but for arbitrary matrices, i.e. not only square ones
fp_t * multiply_with_vector(fp_t * A, fp_t ** B, int dim);
fp_t ** transpose (fp_t ** A, int N_rows, int N_columns);
fp_t **** transpose(fp_t ****A, int N1, int N2, int N3, int N4);
fp_t ** multiply_with_transpose(fp_t ** A, int N_rows, int N_columns);
fp_t ** make_from_diagonal(fp_t * diag, int from, int to);

fp_t * solve(fp_t ** A, fp_t * rhs, int from, int to); // Solve linear system using LU decomposition

// Milic: 09/02/2015: From now on, let us document things a bit. 

fp_t gamma_f(fp_t arg); // No real nneed, it will be used seldom enough that built in c function will do just fine 

fp_t max_1d(fp_t * array, int begin, int end);
fp_t min_1d(fp_t * array, int begin, int end);
fp_t max_2d(fp_t ** array, int beginx, int endx, int beginy, int endy);
fp_t min_2d(fp_t ** array, int beginx, int endx, int beginy, int endy);
int max_1d_index(fp_t * array, int begin, int end);
int max_1d_index(int * array, int begin, int end);

fp_t aps(fp_t); // Why is abs always so problematic?

void Crout(int d,fp_t*S,fp_t*D);
void solveCrout(int d,fp_t*LU,fp_t*b,fp_t*x);
void Crout(int d,__float128*S,__float128*D);
void solveCrout(int d,__float128*LU,__float128*b,__float128*x);

// Routines regarding the computation of the exponential integral:

double      Exponential_Integral_Ei( double x );
long double xExponential_Integral_Ei( long double x );

static long double Continued_Fraction_Ei( long double x );
static long double Power_Series_Ei( long double x );
static long double Argument_Addition_Series_Ei( long double x);

fp_t interpol_2d(fp_t ** f, fp_t * x, fp_t * y, int Nx, int Ny, fp_t x_to_interpol, fp_t y_to_interpol);
fp_t interpol_2d(fp_t * f, fp_t * x, fp_t * y, int Nx, int Ny, fp_t x_to_interpol, fp_t y_to_interpol);

fp_t interpol_1d(fp_t * f, fp_t * x, int N, fp_t x_to_interpol);
fp_t interpol_1d_linear(fp_t * f, fp_t * x, int N, fp_t x_to_interpol);

fp_t * add_to_1d_array(fp_t * x, int &n, fp_t to_add);

int atmospheric_interpolation(fp_t *, fp_t *, int, fp_t *, fp_t *, int, int, int);

fp_t Planck_f(fp_t lambda, fp_t T);
fp_t Planck_f_derivative(fp_t lambda, fp_t T);
fp_t compute_E(fp_t y0, fp_t y1, fp_t y2, fp_t h1, fp_t h2);
fp_t compute_F(fp_t y0, fp_t y1, fp_t y2, fp_t h1, fp_t h2);

fp_t compute_w0(fp_t);
fp_t compute_w1(fp_t);

fp_t vactoair(fp_t);
fp_t airtovac(fp_t);
fp_t * airtovac(fp_t *, int);
fp_t * vactoair(fp_t *, int);

int convolve_spectra_with_gauss(fp_t **, fp_t *,int, fp_t);
int convolve_spectra_with_psf(fp_t **, fp_t *,int, int, fp_t *);

int convolve_response_with_gauss(fp_t *** response, fp_t * lambda, int N_parameters, int N_lambda, fp_t);
int convolve_response_with_psf(fp_t *** response, fp_t * lambda, int N_parameters, int N_lambda, int, fp_t*);
int convolve_response_with_gauss_tau(fp_t ***, int, int, fp_t *, int, fp_t);

int convolve_with_gauss(fp_t *,fp_t *, int,fp_t with);

int set_to_zero_except(fp_t * x, int N, int to_keep);
int stupid_sort_indices_by_abs(fp_t * A, int * indices, int N);
int index_of_max_abs(fp_t * x, int N);

fp_t * svd_treshold_square_matrix(fp_t ** A, int N, fp_t treshold, int return_w);
fp_t * svd_invert_square_matrix(fp_t ** A, int N, fp_t treshold, int return_w);
int shifted_svd(fp_t ** U, int M, int N, fp_t * w, fp_t **V);

double P0(double);
double P1(double);
double P2(double);
double Pn(unsigned int, double);
fp_t * project_on_legendre_basis(fp_t * y, fp_t * x, int N, int up_to);
fp_t * reconstruct_from_legendre_basis(fp_t * coeff, fp_t * x, int N, int up_to);

static double PYTHAG(double a, double b);
int dsvd(double **a, int m, int n, double *w, double **v);

// Various merit functions and such:
fp_t chi_sqr(fp_t ** obs, fp_t ** fit, fp_t * noise, int nlambda, int * stokes_to_fit, int n_stokes_to_fit, fp_t * ws);

#endif                        // __MATHTOOLS_H__
