/*=============================================================================

 Solver.cpp

 Copyright (C) 2015 by Yun Ling

 Clean Version

 =============================================================================*/

 // breakable implementation
 
 // changes: 
 
 // parameter save change
 // initialized parameters change
 // imputed values input

/*
 -----------------------------------------------------------------------------
 Headers
 -----------------------------------------------------------------------------
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <ctime>
#include <map>
#include <set>
#include <numeric>
#include <algorithm>
#include <limits>
#include <omp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
using namespace std;

/*
 -----------------------------------------------------------------------------
 Function
 -----------------------------------------------------------------------------
 */

// (1) Sampling Truncated Normals
// mu = 0, sigma = 1, (-\infty, u]
double randn_tu(const double u, gsl_rng * r) {
	double x;
	if (u > 0.0) {
		do {
			x = gsl_ran_ugaussian_ratio_method(r);
		} while (x > u);
	} else if (u > -0.6) {
		do {
			x = gsl_ran_ugaussian_ratio_method(r);
		} while (x > u && x < -u);
		if (x >= -u) {
			x = -x;
		}
	} else {
		do {
			x = log(gsl_rng_uniform(r)) / u;
		} while (gsl_rng_uniform(r) > exp(-x * x / 2.0));
		x = u - x;
	}
	return x;
}

// mu = 0, sigma = 1, [l, \infty)
double randn_tl(const double l, gsl_rng * r) {
	return -randn_tu(-l, r);
}

// mu = 0, sigma = 1, [l, u]
double randn_tlu(const double l, const double u, gsl_rng * r) {
	double x;
	if (u > 0.0) {
		if (l < 0.0) {
			do {
				x = gsl_rng_uniform(r) * (u - l) + l;
			} while (gsl_rng_uniform(r) > exp(-x * x / 2.0) );
		} else {
			do {
				x = gsl_rng_uniform(r) * (u - l) + l;
			} while (gsl_rng_uniform(r) > exp((l * l - x * x) / 2.0) );
		}
	} else {
		do {
			x = gsl_rng_uniform(r) * (u - l) + l;
		} while (gsl_rng_uniform(r) > exp((u * u - x * x)/ 2.0) );
	}
	return x;
}

// mu = mu, sigma = sigma, (-\infty, u]
double rn_tu(const double u, const double mu, const double sigma2, gsl_rng * r) {
    double sigma = sqrt(sigma2);
	double uu = (u - mu) / sigma;
	return mu + sigma * randn_tu(uu, r);
}

// mu = mu, sigma = sigma, [l, \infty)
double rn_tl(const double l, const double mu, const double sigma2, gsl_rng * r) {
    double sigma = sqrt(sigma2);
	double ll = (l - mu) / sigma;
	return mu + sigma * randn_tl(ll, r);
}

// mu = mu, sigma = sigma, [l, u]
double rn_tlu(const double l, const double u, const double mu, const double sigma2, gsl_rng * r) {
    double sigma = sqrt(sigma2);
	double ll = (l - mu) / sigma;
	double uu = (u - mu) / sigma;
	return mu + sigma * randn_tlu(ll, uu, r);
}

// (2) Bayesian Regressions

void Bayesupdate_sigknown(const gsl_matrix* X, const gsl_vector* y, gsl_matrix* A, gsl_vector* B) {
	// rule: A1 = A0 + X'*X, B1 = B0 + X'*y

	// A1 = A + X'*X	
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, X, X, 1.0, A);	
	// B1 = B + X'*y
	gsl_blas_dgemv(CblasTrans, 1.0, X, y, 1.0, B);
}

void Bayesdraw_sigknown(const gsl_matrix* A, const gsl_vector* B, const int k, const double sigma2, gsl_vector* phi, gsl_rng * r) {
	// phi: ~ N(inv(A)*B, inv(A)*sigma2)
	
	gsl_matrix * temp_mat = gsl_matrix_alloc(k, k);
	gsl_matrix_memcpy(temp_mat, A);
	gsl_linalg_cholesky_decomp(temp_mat);
	gsl_linalg_cholesky_solve(temp_mat, B, phi);
	// z ~ N(0,I)
	gsl_vector * z = gsl_vector_alloc(k);
	for (int i = 0; i < k; i++) {
		gsl_vector_set(z, i, gsl_ran_ugaussian_ratio_method(r));
	}
	// phi + [sigma*inv(A1^(1/2))] * z
	gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, A, z);
	gsl_vector_scale(z, sqrt(sigma2));
	gsl_vector_add(phi, z);	
	// free
	gsl_vector_free(z);
}

void Bayesupdate(const gsl_matrix* X, const gsl_vector* y, gsl_matrix* A, gsl_vector* B, double& a, double& b) {
	// rule: A1 = A0 + X'*X, B1 = B0 + X'*y, a1 = a0 + n/2, b1 = 1.0 / ((1.0 / b0) + ee + bb)

	int n = X->size1;
	int k = X->size2;
	
	gsl_matrix * XX = gsl_matrix_calloc(k, k);
	gsl_vector * beta0 = gsl_vector_calloc(k);
	gsl_vector * beta1 = gsl_vector_calloc(k);
	gsl_vector * temp  = gsl_vector_calloc(k);
	gsl_vector * e     = gsl_vector_calloc(n);
	
	// a1 = a0 + n/2
	a = a + n/2.0;	
	
	// beta0 = inv(A0)*B0
	gsl_matrix_memcpy (XX, A);
	gsl_linalg_cholesky_decomp (XX);
	gsl_linalg_cholesky_solve (XX, B, beta0);	
	
	// A1 = A0 + X'*X
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, X, X, 1.0, A);	
	// B1 = B0 + X'*y
	gsl_blas_dgemv(CblasTrans, 1.0, X, y, 1.0, B);
	
	// beta1 = inv(A1)*B1
	gsl_matrix_memcpy (XX, A);
	gsl_linalg_cholesky_decomp (XX);
	gsl_linalg_cholesky_solve (XX, B, beta1);
	
	// ee = e'*e, e = Y - X*beta1 
	double ee;
	gsl_vector_memcpy (e, y);
	gsl_blas_dgemv (CblasNoTrans, -1.0, X, beta1, 1.0, e);
	gsl_blas_ddot (e, e, &ee);
	// bb = (beta1-beta0)'*A1*(beta1-beta0)
	double bb;
	gsl_vector_sub (beta1, beta0);	
	gsl_blas_dsymv (CblasUpper, 1.0, A, beta1, 0.0, temp);
	gsl_blas_ddot (temp, beta1, &bb);	
	// b1 = 1.0 / ((1.0 / b0) + (ee + bb)/2)
	b = 1.0 / ((1.0 / b) + (ee + bb)/2.0);
	
	// free
	gsl_matrix_free (XX);
	gsl_vector_free (beta0);
	gsl_vector_free (beta1);
	gsl_vector_free (temp);
	gsl_vector_free (e);

}

void Bayesdraw(const gsl_matrix* A, const gsl_vector* B, const double a, const double b, const int k, double& sigma2, gsl_vector* phi, gsl_rng * r) {
	// sigma2 ~ IG(a, b), phi: ~ N(inv(A)*B, inv(A)*sigma2)
	
	gsl_matrix * temp_mat = gsl_matrix_alloc(k, k);
	gsl_matrix_memcpy(temp_mat, A);
	gsl_linalg_cholesky_decomp(temp_mat);
	gsl_linalg_cholesky_solve(temp_mat, B, phi);
	// sigma2
	sigma2 = 1.0 / gsl_ran_gamma(r, a, b);
	// z ~ N(0,I)
	gsl_vector * z = gsl_vector_alloc(k);
	for (int i = 0; i < k; i++) {
		gsl_vector_set(z, i, gsl_ran_ugaussian_ratio_method(r));
	}	
	// phi + [sigma*inv(A1^(1/2))] * z
	gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, A, z);
	gsl_vector_scale(z, sqrt(sigma2));
	gsl_vector_add(phi, z);	
	// free
	gsl_vector_free(z);	
}

// (3) vector -> gsl_vector or gsl_matrix

void gsl_vector_from(gsl_vector * tempv, const vector<double> x, const int start) {
	for (unsigned int k = 0; k < x.size(); k++) {
		gsl_vector_set(tempv, start + k, x[k]);
	}
}

void gsl_vector_from_to(gsl_vector * tempv, const vector<double> x, const int start, const int end) {
	for (unsigned int k = 0; k < end - start; k++) {
		gsl_vector_set(tempv, start + k, x[k]);
	}
}

void gsl_matrix_from(gsl_matrix * tempm, const vector<double> x, const int row, const int start) {
	for (unsigned int k = 0; k < x.size(); k++) {
		gsl_matrix_set(tempm, row, start + k, x[k]);
	}
}

void gsl_matrix_from_to(gsl_matrix * tempm, const vector<double> x, const int row, const int start, const int end) {
	for (unsigned int k = 0; k < end - start; k++) {
		gsl_matrix_set(tempm, row, start + k, x[k]);
	}
}

// (4) input

void gsl_matrix_input(gsl_matrix * tempm, const char * fname, const int G, const int K) {
	ifstream in;
	in.open(fname);
	double temp;
	for (int g = 0; g < G; g++) {
		for (int k = 0; k < K; k ++) {
			in >> temp;
			gsl_matrix_set(tempm, g, k, temp);
		}
	}
	in.close();
}

void vector_input(vector<double> &x, const char * fname, int G) {
	ifstream in;
	in.open(fname);
	double temp;
	for (int g = 0; g < G; g++) {
		in >> temp;
		x[g] = temp;
	}
	in.close();
}

// (5) output

void gsl_matrix_output(const gsl_matrix * tempm, const char * fname, const int G, const int K) {
	ofstream out;
	out.open(fname);
	for (int g = 0; g < G; g++) {
		for (int k = 0; k < K - 1; k++) {
			out << gsl_matrix_get(tempm, g, k) << "\t";
		}
		out << gsl_matrix_get(tempm, g, K - 1) << "\n";
	}
	out.close();
}

void vector_output(const vector<double> x, const char * fname, int G) {
	ofstream out;
	out.open(fname);
	for (int g = 0; g < G; g++) {
		out << x[g] << "\n";
	}
	out.close();
}

/*
 -----------------------------------------------------------------------------
 Main
 -----------------------------------------------------------------------------
 */

 // first is the start iteration
 // second is the end iteration
 
 int main(int argc, char** argv)

{
	int G1 = atoi(argv[1]); // starting iteration (0)
	int G2 = atoi(argv[2]); // ending iteration
	double delta = 3; // for upper and lower truncated normal
	int percent = 1; // display progress
	double SparseMax; // # table_tjs_sparse
	
	// ------------ Metadata ------------ //
	ifstream in;
	in.open("metadata.txt");
	double temp;

	// I, J, T
	int I, J, T;
	in >> temp;
	I = (int) temp;
	in >> temp;
	J = (int) temp;
	in >> temp;
	T = (int) temp;
	
	// K's (new)
	int K_r, K_r_y, K_r_n, K_s, K_s_y, K_s_n;
	in >> temp;
	K_r = (int) temp;
	in >> temp;
	K_r_y = (int) temp;
	in >> temp;
	K_r_n = (int) temp;
	in >> temp;
	K_s = (int) temp;
	in >> temp;
	K_s_y = (int) temp;
	in >> temp;
	K_s_n = (int) temp;

	// K's
	int K_t, K_j, K_tj, K_i, K_ti, K_ji, K_tji;
	in >> temp;
	K_t = (int) temp;
	in >> temp;
	K_j = (int) temp;
	in >> temp;
	K_tj = (int) temp;
	in >> temp;
	K_i = (int) temp;
	in >> temp;
	K_ti = (int) temp;
	in >> temp;
	K_ji = (int) temp;
	in >> temp;
	K_tji = (int) temp;
	in >> temp;
	SparseMax = temp;

	in.close();

	// K's
	K_r_y = K_r_y + 1;
	K_r_n = K_r_n + 1;
	K_s_y = K_s_y + 2;
	K_s_n = K_s_n + 2;
	int K_v = (K_j + K_tj) + (K_i + K_ti + K_ji + K_tji) + 1; // no constant, no t-vars
	cout << "K_r_y: " << K_r_y << endl;
	cout << "K_r_n: " << K_r_n << endl;
	cout << "K_s_y: " << K_s_y << endl;
	cout << "K_s_n: " << K_s_n << endl;
	cout << "K_v: " << K_v << endl;
	
	// ------------ Hyperparameter ------------ //

	// mu
	gsl_vector * mu_r_y = gsl_vector_calloc(K_r_y); // initialize all to zero
	gsl_vector * mu_r_n = gsl_vector_calloc(K_r_n);
	gsl_vector * mu_s_y = gsl_vector_calloc(K_s_y);
	gsl_vector * mu_s_n = gsl_vector_calloc(K_s_n);
	gsl_vector * mu_v = gsl_vector_calloc(K_v);
	// A
	gsl_matrix * A_r_y = gsl_matrix_calloc(K_r_y, K_r_y);
	for (int k = 0; k < K_r_y; k++) {
		gsl_matrix_set(A_r_y, k, k, 1.0/100.0);
	}
	gsl_matrix * A_r_n = gsl_matrix_calloc(K_r_n, K_r_n);
	for (int k = 0; k < K_r_n; k++) {
		gsl_matrix_set(A_r_n, k, k, 1.0/100.0);
	}
	gsl_matrix * A_s_y = gsl_matrix_calloc(K_s_y, K_s_y);
	for (int k = 0; k < K_s_y; k++) {
		gsl_matrix_set(A_s_y, k, k, 1.0/100.0);
	}
	gsl_matrix * A_s_n = gsl_matrix_calloc(K_s_n, K_s_n);
	for (int k = 0; k < K_s_n; k++) {
		gsl_matrix_set(A_s_n, k, k, 1.0/100.0);
	}
	gsl_matrix * A_v = gsl_matrix_calloc(K_v, K_v);
	for (int k = 0; k < K_v; k++) {
		gsl_matrix_set(A_v, k, k, 1.0/100.0);
	}
	// (a, b)
	double a_r_y = 2.0;
	double b_r_y = 1.0;
	double a_r_n = 2.0;
	double b_r_n = 1.0;
	double a_s_y = 2.0;
	double b_s_y = 1.0;
	double a_s_n = 2.0;
	double b_s_n = 1.0;

	// ------------ Parameter Initialization ------------ //

	// (1) allocate space to save parameters
	// phi_save
	gsl_matrix * phi_r_y_save = gsl_matrix_alloc(G2, K_r_y);
	gsl_matrix * phi_r_n_save = gsl_matrix_alloc(G2, K_r_n);
	gsl_matrix * phi_s_y_save = gsl_matrix_alloc(G2, K_s_y);
	gsl_matrix * phi_s_n_save = gsl_matrix_alloc(G2, K_s_n);
	gsl_matrix * phi_v_save = gsl_matrix_alloc(G2, K_v);
	// sigma2_save
	vector<double> sigma2_r_y_save(G2);
	vector<double> sigma2_r_n_save(G2);
	vector<double> sigma2_s_y_save(G2);
	vector<double> sigma2_s_n_save(G2);	
	
	// (2) input estimated parameters (in to save)
	if (G1 > 0) {
		gsl_matrix_input(phi_r_y_save, "phi_r_y.txt", G1, K_r_y);
		gsl_matrix_input(phi_r_n_save, "phi_r_n.txt", G1, K_r_n);
		gsl_matrix_input(phi_s_y_save, "phi_s_y.txt", G1, K_s_y);
		gsl_matrix_input(phi_s_n_save, "phi_s_n.txt", G1, K_s_n);
		gsl_matrix_input(phi_v_save, "phi_v.txt", G1, K_v);
		vector_input(sigma2_r_y_save, "sigma2_r_y.txt", G1);
		vector_input(sigma2_r_n_save, "sigma2_r_n.txt", G1);
		vector_input(sigma2_s_y_save, "sigma2_s_y.txt", G1);
		vector_input(sigma2_s_n_save, "sigma2_s_n.txt", G1);
	}
	
	// (3) starting values
	gsl_vector * phi_r_y = gsl_vector_alloc(K_r_y);
	gsl_vector * phi_r_n = gsl_vector_alloc(K_r_n);
	gsl_vector * phi_s_y = gsl_vector_alloc(K_s_y);
	gsl_vector * phi_s_n = gsl_vector_alloc(K_s_n);
	gsl_vector * phi_v = gsl_vector_alloc(K_v);
	double sigma2_r_y, sigma2_r_n, sigma2_s_y, sigma2_s_n;
	if (G1 > 0) {
		gsl_vector_view phi_r_y_save_g = gsl_matrix_row(phi_r_y_save, G1-1);
		gsl_vector_memcpy(phi_r_y, &phi_r_y_save_g.vector);
		gsl_vector_view phi_r_n_save_g = gsl_matrix_row(phi_r_n_save, G1-1);
		gsl_vector_memcpy(phi_r_n, &phi_r_n_save_g.vector);	
		gsl_vector_view phi_s_y_save_g = gsl_matrix_row(phi_s_y_save, G1-1);
		gsl_vector_memcpy(phi_s_y, &phi_s_y_save_g.vector);
		gsl_vector_view phi_s_n_save_g = gsl_matrix_row(phi_s_n_save, G1-1);
		gsl_vector_memcpy(phi_s_n, &phi_s_n_save_g.vector);		
		gsl_vector_view phi_v_save_g = gsl_matrix_row(phi_v_save, G1-1);
		gsl_vector_memcpy(phi_v, &phi_v_save_g.vector);	
		sigma2_r_y = sigma2_r_y_save[G1-1];
		sigma2_r_n = sigma2_r_n_save[G1-1];
		sigma2_s_y = sigma2_s_y_save[G1-1];
		sigma2_s_n = sigma2_s_n_save[G1-1];
	}
	else {
		gsl_vector_memcpy(phi_r_y, mu_r_y);
		gsl_vector_memcpy(phi_r_n, mu_r_n);
		gsl_vector_memcpy(phi_s_y, mu_s_y);
		gsl_vector_memcpy(phi_s_n, mu_s_n);
		sigma2_r_y = 1.0;
		sigma2_r_n = 1.0;
		sigma2_s_y = 1.0;
		sigma2_s_n = 1.0;
	}


	/* 2. Input Data */

	// ------------ Startup Table ------------ //
	
	cout << "Starting Input...\n";
	
	// Table0: t, j, status, match, r -> index0, status, match, r, (born, dead)
	int Nrows = 0;
	vector<int> index_0;
	map<int, int> born; // document when j is born
	map<int, int> dead; // document when j is dead
	map<int, int> status;
	map<int, int> match;
	map<int, double> r_obs;
	// J = U (unmatched) + M (matched)
	// E = J + N
	map<int, vector<int> > SetJ_U;
	map<int, vector<int> > SetJ_M;
	map<int, vector<int> > SetN;
	map<int, vector<int> > SetD;

	ifstream in_0("table_0_b.txt");
	while (in_0 >> temp) {
		// t
		int t = (int) temp;
		// j
		in_0 >> temp;
		int j = (int) temp;
		// index_0
		index_0.push_back(t * J + j);
		// status
		in_0 >> temp;
		status[t * J + j] = (int) temp;
		double temp2;
		// match
		in_0 >> temp2;
		match[t * J + j] = (int) temp2;
		// Set: status = -2 (N), 0 (J), -1/1 (D)
		if (temp == -2) {
			SetN[t].push_back(j);
			born[j] = t;
			// fake death at t = T-1 for j, but not included in SetD
			if (t == T - 1) {
				dead[j] = T - 1;
			}
		} else if (temp == 0) {
			// SetJ_U and SetJ_M
			if (temp2 == -1) {
				SetJ_U[t].push_back(j);
			} else {
				SetJ_M[t].push_back(j);
			}
			// fake death at t = T-1 for j, but not included in SetD
			if (t == T - 1) {
				dead[j] = T - 1;
			}
		} else {
			SetD[t].push_back(j);
			dead[j] = t;
		}
		// r
		in_0 >> temp;
		r_obs[t * J + j] = temp;
		Nrows++;
	}
	in_0.close();	
	cout << "finishing table_0\n";	
	
	// ------------ Independent Tables ------------ //
	
	// Table_r
	vector<int> index_r;
	map<int, vector<double> > indep_r;
	Nrows = 0;
	ifstream in_r("table_r.txt");
	while (in_r >> temp) {
		// t
		int t = (int) temp;
		// j
		in_r >> temp;
		int j = (int) temp;
		// index_r
		index_r.push_back(t * J + j);
		// indep_r
		for (int k = 0; k < K_r; k++) {
			in_r >> temp;
			indep_r[t * J + j].push_back(temp);
		}
		Nrows++;
	}
	in_r.close();	
	cout << "finishing table_r\n";
	
	// Table_s
	vector<int> index_s;
	map<int, vector<double> > indep_s;
	Nrows = 0;
	ifstream in_s("table_s.txt");
	while (in_s >> temp) {
		// t
		int t = (int) temp;
		// j
		in_s >> temp;
		int j = (int) temp;
		// index_s
		index_s.push_back(t * J + j);
		// indep_s
		for (int k = 0; k < K_s; k++) {
			in_s >> temp;
			indep_s[t * J + j].push_back(temp);
		}
		Nrows++;
	}
	in_s.close();	
	cout << "finishing table_s\n";

	// Table_t
	vector<int> index_t;
	map<int, vector<double> > indep_t;
	Nrows = 0;
	ifstream in_t("table_t.txt");
	while (in_t >> temp) {
		// t
		int t = (int) temp;
		// index_t
		index_t.push_back(t);
		// indep_t
		for (int k = 0; k < K_t; k++) {
			in_t >> temp;
			indep_t[t].push_back(temp);
		}
		Nrows++;
	}
	in_t.close();	
	cout << "finishing table_t\n";

	// Table_j
	vector<int> index_j;
	map<int, vector<double> > indep_j;
	Nrows = 0;
	ifstream in_j("table_j.txt");
	while (in_j >> temp) {
		// j
		int j = (int) temp;
		// index_j
		index_j.push_back(j);
		// indep_j
		for (int k = 0; k < K_j; k++) {
			in_j >> temp;
			indep_j[j].push_back(temp);
		}
		Nrows++;
	}
	in_j.close();	
	cout << "finishing table_j\n";

	// Table_tj
	vector<int> index_tj;
	map<int, vector<double> > indep_tj;
	Nrows = 0;
	ifstream in_tj("table_tj.txt");
	while (in_tj >> temp) {
		// t
		int t = (int) temp;
		// j
		in_tj >> temp;
		int j = (int) temp;
		// index_tj
		index_tj.push_back(t * J + j);
		// indep_tj
		for (int k = 0; k < K_tj; k++) {
			in_tj >> temp;
			indep_tj[t * J + j].push_back(temp);
		}
		Nrows++;
	}
	in_tj.close();	
	cout << "finishing table_tj\n";

	// table_i
	vector<int> index_i;
	map<int, vector<double> > indep_i;
	Nrows = 0;
	ifstream in_i("table_i.txt");
	while (in_i >> temp) {
		// i
		int i = (int) temp;
		// index_i
		index_i.push_back(i);
		// indep_i
		for (int k = 0; k < K_i; k++) {
			in_i >> temp;
			indep_i[i].push_back(temp);
		}
		Nrows++;
	}
	in_i.close();	
	cout << "finishing table_i\n";

	// table_ti
	vector<int> index_ti;
	map<int, vector<double> > indep_ti;
	Nrows = 0;
	ifstream in_ti("table_ti.txt");
	while (in_ti >> temp) {
		// t
		int t = (int) temp;
		// i
		in_ti >> temp;
		int i = (int) temp;
		// index_ti
		index_ti.push_back(t * I + i);
		// indep_ti
		for (int k = 0; k < K_ti; k++) {
			in_ti >> temp;
			indep_ti[t * I + i].push_back(temp);
		}
		Nrows++;
	}
	in_ti.close();	
	cout << "finishing table_ti\n";

	// table_ji
	vector<int> index_ji;
	map<int, vector<double> > indep_ji;
	Nrows = 0;
	ifstream in_ji("table_ji.txt");
	while (in_ji >> temp) {
		// j
		int j = (int) temp;
		// i
		in_ji >> temp;
		int i = (int) temp;
		// index_ji
		index_ji.push_back(j * I + i);
		// indep_ji
		for (int k = 0; k < K_ji; k++) {
			in_ji >> temp;
			indep_ji[j * I + i].push_back(temp);
		}
		Nrows++;
	}
	in_ji.close();	
	cout << "finishing table_ji\n";
	cout << "SparseMax: " << SparseMax << "\n";

	// table_tji: sparse
	vector<int> index_tji;
	map<int, vector<double> > indep_tji;
	Nrows = 0;
	ifstream in_tji("table_tji_sparse.txt");
	while (in_tji >> temp) {
		// t
		int t = (int) temp;
		// j
		in_tji >> temp;
		int j = (int) temp;
		// i
		in_tji >> temp;
		int i = (int) temp;
		// index_tji
		index_tji.push_back((t * J + j) * I + i);
		// indep_tji
		for (int k = 0; k < K_tji; k++) {
			in_tji >> temp;
			indep_tji[(t * J + j) * I + i].push_back(temp);
		}
		Nrows++;
		// print percentage of completed input 
        if (((double)Nrows)/((double)SparseMax) >= ((double)percent)/100.0)
        {
            // cout << percent << "% complete" << endl;
            percent += 1;
        }
	}
	in_tji.close();	
	cout << "finishing table_tji\n";

	/* 3. Gibbs Sampler */

	// ------------ Unobservable Initialization ------------ //

	// gsl_r: random number generator
	gsl_rng_env_setup();
	const gsl_rng_type * gsl_T = gsl_rng_ranlxd2;
	gsl_rng * gsl_r = gsl_rng_alloc(gsl_T);
	
    cout << "\nInitialization\n";
	map<int, double> r;
	map<int, double> s;
	map<int, double> v;

	for (int j = 0; j < J; j++) {
		for (int t = born[j] + 1; t <= dead[j]; t++) {
			// s: regardless of the ultimate exit status
			s[t * J + j] = 0.0;
			// v: only for those not exit yet (in J(t))
			if (t != dead[j]) {
				for (int i = 0; i < I; i++) {
					v[(t * J + j) * I + i] = 0.0;
				}
			}
		}
	}	

	if (G1 > 0) {
		ifstream in_s("imputed_s.txt");
		while (in_s >> temp) {
			int j = (int) temp;
			in_s >> temp;
			int t = (int) temp;
			in_s >> temp;
			s[t * J + j] = temp;
		}
		in_s.close();
		ifstream in_v("imputed_v.txt");
		while (in_v >> temp) {
			int t = (int) temp;
			in_v >> temp;
			int j = (int) temp;
			in_v >> temp;
			int i = (int) temp;
			in_v >> temp;
			v[(t * J + j) * I + i] = temp;
		}
		in_v.close();
	}

	cout << "finishing initalizing v..." << endl;
	
	cout << "\nCalculating..." << endl;
	cout << "From iteration " << G1 << " to iteration " << G2 << endl;	

	for (int g = G1; g < G2; g++) {
        
        cout << "iter " << g << endl;
        
		// ------------ Impute Unobservables (r,s,v) ------------ //

		// (1) ------------------- Impute r -------------------
		cout << "impute r" << endl;
		#pragma omp parallel for
		for (int j = 0; j < J; j++) {
			
			int dur = dead[j] - born[j] + 1;
			double EX, VX, EXfw, VXfw, dE, dV;
			double EXtt[dur];
			double VXtt[dur];
			double EXfwtt[dur];
			
			// Step A) Filter Forward to get the mean and variance of v(t|t), stored
			//		   in EXtt and VXtt
			for (int t = born[j]; t <= dead[j]; t++) {
				int tau = t - born[j];
				
				// A1) The updating step
				// observable
				if (r_obs[t * J + j] != -999.0) {
					EX = r_obs[t * J + j];
					VX = 0.0;
				}
				else {
					gsl_vector * b = gsl_vector_calloc(2);
					gsl_matrix * p = gsl_matrix_calloc(2,2);
					gsl_vector * e = gsl_vector_calloc(2);
					gsl_vector * KF = gsl_vector_calloc(2);
					// b = (phi(s,-1), [phi(v,-1])
					double phi_s = match[(t - 1) * J + j] == -1 ? gsl_vector_get(phi_s_n, K_s_n - 1) : gsl_vector_get(phi_s_y, K_s_y - 1);
					gsl_vector_set(b, 0, phi_s);
					gsl_vector_set(b, 1, gsl_vector_get(phi_v, K_v - 1));
					// p = diag(sigma2_s, [sigma2_v]/I) 
					double sigma2_s = match[(t - 1) * J + j] == -1 ? sigma2_s_n : sigma2_s_y;
					gsl_matrix_set(p, 0, 0, sigma2_s);
					gsl_matrix_set(p, 1, 1, 1.0/I);
					// e = (s-s(EXfw), mean(v-v(EXfw))
					double ds;					
					if (match[(t - 1) * J + j] == -1) {
						gsl_vector * tempv = gsl_vector_calloc(K_s_n);
						gsl_vector_set(tempv, 0, 1.0);
						gsl_vector_set(tempv, K_s_n - 1, EXfw);
						gsl_vector_from_to(tempv, indep_s[t * J + j], 1, K_s_n - 1);
						gsl_blas_ddot(tempv, phi_s_n, &ds);
						gsl_vector_free(tempv);
					}
					else {
						int i = match[(t - 1) * J + j];
						gsl_vector * tempv = gsl_vector_calloc(K_s_y);
						gsl_vector_set(tempv, 0, 1.0);
						gsl_vector_set(tempv, K_s_y - 1, EXfw);
						gsl_vector_from_to(tempv, indep_s[t * J + j], 1, K_s_y - 1);
						gsl_blas_ddot(tempv, phi_s_y, &ds);
						gsl_vector_free(tempv);				
					}

					gsl_vector_set(e, 0, s[t * J + j] - ds);
					double e_v = 0.0;
					for (int i = 0; i < I; i++) {
						double dv;
						gsl_vector * tempv = gsl_vector_calloc(K_v);
						gsl_vector_set (tempv, K_v - 1, EXfw);	
						gsl_vector_from(tempv, indep_j[j], 0);
						gsl_vector_from(tempv, indep_tj[t * J + j], K_j);
						gsl_vector_from(tempv, indep_i[i], K_j + K_tj);
						gsl_vector_from(tempv, indep_ti[t * I + i], K_j + K_tj + K_i);
						gsl_vector_from(tempv, indep_ji[j * I + i], K_j + K_tj + K_i + K_ti);
						if (indep_tji.find((t * J + j) * I + i)!=indep_tji.end()) {
							gsl_vector_from(tempv, indep_tji[(t * J + j) * I + i], K_j + K_tj + K_i + K_ti + K_ji);
		                }
						gsl_blas_ddot(tempv, phi_v, &dv);						
						gsl_vector_free(tempv);
						e_v += v[(t * J + j) * I + i] - dv;
					}
		            gsl_vector_set (e, 1, e_v/I);			
					// KF = inv(b*VXfw*b' + p) * (VXfw*b)
					gsl_blas_dsyr(CblasUpper, VXfw, b, p);
					gsl_matrix_set (p, 1, 0, gsl_matrix_get(p, 0, 1));
					gsl_linalg_cholesky_decomp(p);
					gsl_linalg_cholesky_solve(p, b, KF);
					gsl_vector_scale(KF, VXfw);			
					// EX = EXfw + KF'*e
					double KFe;
					gsl_blas_ddot(KF, e, &KFe);
					EX = EXfw + KFe;
					// EV = VXfw - VXfw*KF*b
					double KFb;
					gsl_blas_ddot(KF, b, &KFb);
					VX = (1 - min(KFb,0.99)) * VXfw;	
					// free
					gsl_vector_free (b);
					gsl_matrix_free (p);
					gsl_vector_free (e);
					gsl_vector_free (KF);
				}
				// A2) The forecasting step
				if (t < dead[j]) {
					// previously unmatched
					int i = match[t * J + j];
					if (i == -1) {
						gsl_vector * tempv = gsl_vector_calloc(K_r_n);
						gsl_vector_set(tempv, 0, 1.0);
						gsl_vector_from_to(tempv, indep_r[(t + 1) * J + j], 1, K_r_n);
						gsl_blas_ddot(tempv, phi_r_n, &dE);
						gsl_vector_free(tempv);
						dV = sigma2_r_n;
					}
					// previously matched
					else {
						gsl_vector * tempv = gsl_vector_calloc(K_r_y);
						gsl_vector_set(tempv, 0, 1.0);
						gsl_vector_from_to(tempv, indep_r[(t + 1) * J + j], 1, K_r_y);
						gsl_blas_ddot(tempv, phi_r_y, &dE);
						gsl_vector_free(tempv);
						dV = sigma2_r_y;
					}
				}
				// A3) Update parameters
				EXfw = EX + dE;
				VXfw = VX + dV;
				EXtt[tau] = EX;
				EXfwtt[tau] = EXfw;	
				VXtt[tau] = VX;
			}

			// B. Backward
			#pragma omp critical(dataupdate)
			r[dead[j] * J + j] = EXtt[dead[j]-born[j]] + sqrt(VXtt[dead[j]-born[j]]) * gsl_ran_ugaussian_ratio_method (gsl_r);
			for (int t = dead[j]-1; t >= born[j]; t--) {
				int tau = t - born[j];
				if (VXtt[tau] == 0) {
					r[t * J + j] = EXtt[tau];
				}
				else {
					double GF;
					// currently unmatched
					if (match[t * J + j] == -1) {
						GF = VXtt[tau] / (VXtt[tau] + sigma2_r_n);
					}
					// currently matched
					else {
						GF = VXtt[tau] / (VXtt[tau] + sigma2_r_y);
					}
					double mu = EXtt[tau] + GF * (r[(t + 1) * J + j] - EXfwtt[tau]);
					double sigma2 = VXtt[tau] * (1 - GF);
					#pragma omp critical(dataupdate)
					r[t * J + j] = mu + sqrt(sigma2) * gsl_ran_ugaussian_ratio_method(gsl_r);
				}
			}
			
		}
		
		// (2) ------------------- Impute s -------------------
		cout << "impute s" << endl;
		#pragma omp parallel for
		for (int j = 0; j < J; j++) {
			for (int t = born[j] + 1; t <= dead[j]; t++) {
				// (mu, sigma2) for truncated normal
				double mu, sigma2;
				// previsouly unmatch
				int i = match[(t - 1) * J + j];
				if (i == -1) {
					gsl_vector * tempv = gsl_vector_calloc(K_s_n);
					gsl_vector_set(tempv, 0, 1.0);
					gsl_vector_set(tempv, K_s_n - 1, r[t * J + j]);
					gsl_vector_from_to(tempv, indep_s[t * J + j], 1, K_s_n - 1);
					gsl_blas_ddot(tempv, phi_s_n, &mu);
					gsl_vector_free(tempv);
					sigma2 = sigma2_s_n;
				} 
                // previously matched
                else {
					gsl_vector * tempv = gsl_vector_calloc(K_s_y);
					gsl_vector_set(tempv, 0, 1.0);
					gsl_vector_set(tempv, K_s_y - 1, r[t * J + j]);
					gsl_vector_from_to(tempv, indep_s[t * J + j], 1, K_s_y - 1);
					gsl_blas_ddot(tempv, phi_s_y, &mu);
					gsl_vector_free(tempv);
					sigma2 = sigma2_s_y;
				}
				// draw s
				if (status[t * J + j] == 1) {
					#pragma omp critical(dataupdate)
					s[t * J + j] = rn_tl(delta, mu, sigma2, gsl_r);
				} else if (status[t * J + j] == -1) {
					#pragma omp critical(dataupdate)
					s[t * J + j] = rn_tu(-delta, mu, sigma2, gsl_r);
				} else {
					#pragma omp critical(dataupdate)
					s[t * J + j] = rn_tlu(-delta, delta, mu, sigma2, gsl_r);
				}
			}
			
		}		
		
		// (3) ------------------- Impute v -------------------
		cout << "impute v" << endl;

		// create a bound
		double infty = 3.0;
		
		// inner paralell for each j
		for (int t = 1; t < T; t++) {
			//cout << "t: " << t << endl;			

			// setJ
			vector<int> SetJ;
			SetJ.insert(SetJ.end(), SetJ_U[t].begin(), SetJ_U[t].end());
			SetJ.insert(SetJ.end(), SetJ_M[t].begin(), SetJ_M[t].end());
			// (a) inverse_match: i's match
			map<int, vector<int> > i_match;
			for (int iter = 0; iter < SetJ_M[t].size(); iter++) {
				int j = SetJ_M[t][iter];
				int i = match[t * J + j];
				if (i != -1) {
					i_match[i].push_back(j);
				}
			}
			// (b) min{v_i: (i,j) in mu}
			map<int, double> min_v_i;
			for (int i = 0; i < I; i++) {
				if (i_match[i].size() > 0) {
					min_v_i[i] = v[(t * J + i_match[i][0]) * I + i];
					for (int k = 1; k < i_match[i].size(); k++) {
						min_v_i[i] = min(min_v_i[i], v[(t * J + i_match[i][k]) * I + i]);
					}					
				}
			}
			
			// (1) draw impute matched pair next (at most J pairs)
			//cout << "SetJ_M[t].size(): " << SetJ_M[t].size() << endl;
			//#pragma omp parallel for
			for (int iter = 0; iter < SetJ_M[t].size(); iter++) {
				int j = SetJ[iter];
				int i = match[t * J + j];
				// S_i: unpaired j's that prefers i
				vector<int> S_i;
				for (int it = 0; it < SetJ.size(); it++) {
					int k = SetJ[it];
					int i_k = match[t * J + k];
					if (i_k != i) {
						bool add = (i_k == -1) ? true : (v[(t * J + k) * I + i] > v[(t * J + k) * I + i_k] ? true : false); 
						if (add == true) S_i.push_back(k);
					}
				}
				// S_j: unpaired i's that prefers j
				vector<int> S_j;
				for (int k = 0; k < I; k++) {
					if (k != i) {
						bool add = (i_match[k].size() == 0) ? true: (v[(t * J + j) * I + k] > min_v_i[k] ? true : false);
						if (add == true) S_j.push_back(k);
					}
				}
				// v_i = max S_i
				double v_i = S_i.size() > 0 ? v[(t * J + S_i[0]) * I + i] : -infty; 
				for (int k = 1; k < S_i.size(); k++) {
					v_i = max(v_i, v[(t * J + S_i[k]) * I + i]);
				}
				// v_j = max S_j
				double v_j = S_j.size() > 0 ? v[(t * J + j) * I + S_j[0]] : -infty;
				for (int k = 1; k < S_j.size(); k++) {
					v_j = max(v_j, v[(t * J + j) * I + S_j[k]]);
				}
				// v_inf
				double v_inf = max(v_i, v_j);
				// mu
				double mu;
				gsl_vector * tempv = gsl_vector_calloc(K_v);
				gsl_vector_set(tempv, K_v - 1, r[t * J + j]);
				gsl_vector_from(tempv, indep_j[j], 0);
				gsl_vector_from(tempv, indep_tj[t * J + j], K_j);
				gsl_vector_from(tempv, indep_i[i], K_j + K_tj);
				gsl_vector_from(tempv, indep_ti[t * I + i], K_j + K_tj + K_i);
				gsl_vector_from(tempv, indep_ji[j * I + i], K_j + K_tj + K_i + K_ti);
				if (indep_tji.find((t * J + j) * I + i)!=indep_tji.end()) {
					gsl_vector_from(tempv, indep_tji[(t * J + j) * I + i], K_j + K_tj + K_i + K_ti + K_ji);
                }
				gsl_blas_ddot(tempv, phi_v, &mu);
				gsl_vector_free(tempv);
				// draw v
				//#pragma omp critical(dataupdate)
				v[(t * J + j) * I + i] = rn_tl(v_inf, mu, 1.0, gsl_r);
			}
			
			// (2) draw impute unmatched pair first (about J * I pairs)
			//cout << "SetJ_U[t].size(): " << SetJ_U[t].size() << endl;
			#pragma omp parallel for
			for (int iter = 0; iter < SetJ.size(); iter++) {
				int j = SetJ[iter];
				for (int i = 0; i < I; i++) {
					int mu_j = match[t * J + j];
					if (mu_j != i) {
						// v_j
						double v_j = match[t * J + j] > -1 ? v[(t * J + j) * I + mu_j] : -infty; // always prefer to be matched than unmatched
						// v_i
						double v_i = i_match[i].size() > 0 ? min_v_i[i] : -infty; // always prefer to be matched than unmatched
						// v_sup
						double v_sup = max(v_i, v_j);
						// mu
						double mu;
						gsl_vector * tempv = gsl_vector_calloc(K_v);
						gsl_vector_set(tempv, K_v - 1, r[t * J + j]);
						gsl_vector_from(tempv, indep_j[j], 0);
						gsl_vector_from(tempv, indep_tj[t * J + j], K_j);
						gsl_vector_from(tempv, indep_i[i], K_j + K_tj);
						gsl_vector_from(tempv, indep_ti[t * I + i], K_j + K_tj + K_i);
						gsl_vector_from(tempv, indep_ji[j * I + i], K_j + K_tj + K_i + K_ti);
						if (indep_tji.find((t * J + j) * I + i)!=indep_tji.end()) {
							gsl_vector_from(tempv, indep_tji[(t * J + j) * I + i], K_j + K_tj + K_i + K_ti + K_ji);
	                    }
						gsl_blas_ddot(tempv, phi_v, &mu);
						gsl_vector_free(tempv);
						// draw v
						#pragma omp critical(dataupdate)
						v[(t * J + j) * I + i] = rn_tu(v_sup, mu, 1.0, gsl_r);
					}
				}
			}
						
		}
		
		
		// ------------ Update Parameters ------------ //
		
		int count;
		int count_y = 0;
		int count_n = 0;
		int count_v = 0;
		for (int t = 1; t < T; t++) {
			count_y += SetJ_M[t - 1].size();
			count_n += SetJ_U[t - 1].size() + SetN[t - 1].size();
			count_v += SetJ_M[t].size();
			count_v += SetJ_U[t].size();
		}
		count_v *= I;
        
		// (1) (phi_r_y, sigma2_r_y)    
		cout << "update r_y" << endl;
		if (count_y > 0) {
			// initialize
			gsl_matrix * A_r_y_temp = gsl_matrix_calloc(K_r_y, K_r_y);
			gsl_matrix_memcpy (A_r_y_temp , A_r_y);
			gsl_vector * B_r_y_temp = gsl_vector_calloc(K_r_y);
			gsl_blas_dgemv(CblasNoTrans, 1.0, A_r_y, mu_r_y, 1.0, B_r_y_temp);	
			double a_r_y_temp = a_r_y; 
			double b_r_y_temp = b_r_y;
			// update
			#pragma omp parallel for
			for (int t = 1; t < T; t++) {
				// SetE_M[t - 1] = SetJ_M[t - 1]
				int size = SetJ_M[t-1].size();
				if (size > 0) {
					gsl_vector * regy_r_y = gsl_vector_calloc(size);
					gsl_matrix * regX_r_y = gsl_matrix_calloc(size, K_r_y);
					for (unsigned int c = 0; c < size; c++) {
						int j = SetJ_M[t - 1][c];
						int i = match[(t - 1) * J + j];
						gsl_vector_set(regy_r_y, c, r[t * J + j] - r[(t - 1) * J + j]);
						gsl_matrix_set(regX_r_y, c, 0, 1.0);
						gsl_matrix_from_to(regX_r_y, indep_r[t * J + j], c, 1, K_r_y);
					}
					#pragma omp critical(dataupdate)
					Bayesupdate(regX_r_y, regy_r_y, A_r_y_temp, B_r_y_temp, a_r_y_temp, b_r_y_temp);
					gsl_vector_free(regy_r_y);
					gsl_matrix_free(regX_r_y);
				}
			}
			// draw
			Bayesdraw(A_r_y_temp, B_r_y_temp, a_r_y_temp, b_r_y_temp, K_r_y, sigma2_r_y, phi_r_y, gsl_r);
			gsl_matrix_free (A_r_y_temp);
			gsl_vector_free (B_r_y_temp);
		}   
	  
		// (2) (phi_r_n, sigma2_r_n) 
		cout << "update r_n" << endl;	  
		if (count_n > 0) {
			// initialize
			gsl_matrix * A_r_n_temp = gsl_matrix_calloc(K_r_n, K_r_n);
			gsl_matrix_memcpy (A_r_n_temp , A_r_n);
			gsl_vector * B_r_n_temp = gsl_vector_calloc(K_r_n);
			gsl_blas_dgemv(CblasNoTrans, 1.0, A_r_n, mu_r_n, 1.0, B_r_n_temp);	
			double a_r_n_temp = a_r_n; 
			double b_r_n_temp = b_r_n;		
			// update
			#pragma omp parallel for
			for (int t = 1; t < T; t++) {
				// SetE_U[t - 1] = SetM_U[t - 1] + SetN[t - 1]
				vector<int> SetE_U_1;
				SetE_U_1.insert(SetE_U_1.end(), SetJ_U[t-1].begin(), SetJ_U[t-1].end());
				SetE_U_1.insert(SetE_U_1.end(), SetN[t-1].begin(), SetN[t-1].end());
				int size = SetE_U_1.size();
				if (size > 0) {
					gsl_vector * regy_r_n = gsl_vector_calloc(size);
					gsl_matrix * regX_r_n = gsl_matrix_calloc(size, K_r_n);
					for (unsigned int c = 0; c < size; c++) {
						int j = SetE_U_1[c];
						gsl_vector_set(regy_r_n, c, r[t * J + j] - r[(t - 1) * J + j]);
						gsl_matrix_set(regX_r_n, c, 0, 1.0);
						gsl_matrix_from_to(regX_r_n, indep_r[t * J + j], c, 1, K_r_n);				
					}
					#pragma omp critical(dataupdate)
					Bayesupdate(regX_r_n, regy_r_n, A_r_n_temp, B_r_n_temp, a_r_n_temp, b_r_n_temp);
					gsl_vector_free(regy_r_n);
					gsl_matrix_free(regX_r_n);
				}
			}
			// draw
			Bayesdraw(A_r_n_temp, B_r_n_temp, a_r_n_temp, b_r_n_temp, K_r_n, sigma2_r_n, phi_r_n, gsl_r);
			gsl_matrix_free (A_r_n_temp);
			gsl_vector_free (B_r_n_temp);				
		}
			    
		// (3) (phi_s_y, sigma2_s_y)
		cout << "update s_y" << endl;
		if (count_y > 0) {
			// initialize
			gsl_matrix * A_s_y_temp = gsl_matrix_calloc(K_s_y, K_s_y);
			gsl_matrix_memcpy (A_s_y_temp , A_s_y);
			gsl_vector * B_s_y_temp = gsl_vector_calloc(K_s_y);
			gsl_blas_dgemv(CblasNoTrans, 1.0, A_s_y, mu_s_y, 1.0, B_s_y_temp);	
			double a_s_y_temp = a_s_y; 
			double b_s_y_temp = b_s_y;		
			// update
			#pragma omp parallel for
			for (int t = 1; t < T; t++) {
				// SetE_M[t - 1] = SetJ_M[t - 1]
				int size = SetJ_M[t-1].size();
				if (size > 0) {
					gsl_vector * regy_s_y = gsl_vector_calloc(size);
					gsl_matrix * regX_s_y = gsl_matrix_calloc(size, K_s_y);
					for (unsigned int c = 0; c < size; c++) {
						int j = SetJ_M[t - 1][c];
						int i = match[(t - 1) * J + j];
						gsl_vector_set(regy_s_y, c, s[t * J + j]);
						gsl_matrix_set(regX_s_y, c, 0, 1.0);
						gsl_matrix_set(regX_s_y, c, K_s_y - 1, r[t * J + j]);
						gsl_matrix_from_to(regX_s_y, indep_s[t * J + j], c, 1, K_s_y - 1);	
					}
					#pragma omp critical(dataupdate)
					Bayesupdate(regX_s_y, regy_s_y, A_s_y_temp, B_s_y_temp, a_s_y_temp, b_s_y_temp);
					gsl_vector_free(regy_s_y);
					gsl_matrix_free(regX_s_y);
				}
			}				
			// draw
			Bayesdraw(A_s_y_temp, B_s_y_temp, a_s_y_temp, b_s_y_temp, K_s_y, sigma2_s_y, phi_s_y, gsl_r);
			gsl_matrix_free (A_s_y_temp);
			gsl_vector_free (B_s_y_temp);
		}
		
		// (4) (phi_s_n, sigma2_s_n)
		cout << "update s_n" << endl;
		if (count_n > 0) {
			// initialize
			gsl_matrix * A_s_n_temp = gsl_matrix_calloc(K_s_n, K_s_n);
			gsl_matrix_memcpy (A_s_n_temp , A_s_n);
			gsl_vector * B_s_n_temp = gsl_vector_calloc(K_s_n);
			gsl_blas_dgemv(CblasNoTrans, 1.0, A_s_n, mu_s_n, 1.0, B_s_n_temp);	
			double a_s_n_temp = a_s_n; 
			double b_s_n_temp = b_s_n;	
			// update
			#pragma omp parallel for
			for (int t = 1; t < T; t++) {
				// SetE_U[t - 1] = SetJ_U[t - 1] + SetJ_N[t - 1]
				vector<int> SetE_U_1;
				SetE_U_1.insert(SetE_U_1.end(), SetJ_U[t-1].begin(), SetJ_U[t-1].end());
				SetE_U_1.insert(SetE_U_1.end(), SetN[t-1].begin(), SetN[t-1].end());
				int size = SetE_U_1.size();
				if (size > 0) {
					gsl_vector * regy_s_n = gsl_vector_calloc(size);
					gsl_matrix * regX_s_n = gsl_matrix_calloc(size, K_s_n);
					for (unsigned int c = 0; c < size; c++) {
						int j = SetE_U_1[c];
						gsl_vector_set(regy_s_n, c, s[t * J + j]);
						gsl_matrix_set(regX_s_n, c, 0, 1.0);
						gsl_matrix_set(regX_s_n, c, K_s_n - 1, r[t * J + j]);
						gsl_matrix_from_to(regX_s_n, indep_s[t * J + j], c, 1, K_s_n - 1);
						gsl_matrix_from(regX_s_n, indep_t[t], c, 1);				
					}
					#pragma omp critical(dataupdate)
					Bayesupdate(regX_s_n, regy_s_n, A_s_n_temp, B_s_n_temp, a_s_n_temp, b_s_n_temp);
					gsl_vector_free(regy_s_n);
					gsl_matrix_free(regX_s_n);
				}
			}		
			// draw
			Bayesdraw(A_s_n_temp, B_s_n_temp, a_s_n_temp, b_s_n_temp, K_s_n, sigma2_s_n, phi_s_n, gsl_r);
			gsl_matrix_free (A_s_n_temp);
			gsl_vector_free (B_s_n_temp);	
		}		
		
		// (5) (phi_v)
		cout << "update v" << endl;
		if (count_v > 0) {
			// initialize
			gsl_matrix * A_v_temp = gsl_matrix_calloc(K_v, K_v);
			gsl_matrix_memcpy (A_v_temp, A_v);
			gsl_vector * B_v_temp  = gsl_vector_calloc(K_v);
			gsl_blas_dgemv(CblasNoTrans, 1.0, A_v, mu_v, 1.0, B_v_temp);
			// update
			#pragma omp parallel for
			for (int t = 1; t < T; t++) {
				vector<int> SetJ;
				SetJ.insert(SetJ.end(), SetJ_U[t].begin(), SetJ_U[t].end());
				SetJ.insert(SetJ.end(), SetJ_M[t].begin(), SetJ_M[t].end());
				int size = SetJ.size()*I;
				if (size > 0) {
					gsl_vector * regy_v = gsl_vector_calloc(size);
					gsl_matrix * regX_v = gsl_matrix_calloc(size, K_v);
					int count = 0;				
					for (int i = 0; i < I; i++) {
						for (int iter = 0; iter < SetJ.size(); iter++) {
							int j = SetJ[iter];
							gsl_vector_set(regy_v, count, v[(t * J + j) * I + i]);
							gsl_matrix_set(regX_v, count, K_v - 1, r[t * J + j]);
							gsl_matrix_from(regX_v, indep_j[j], count, 0);
							gsl_matrix_from(regX_v, indep_tj[t * J + j], count, K_j);
							gsl_matrix_from(regX_v, indep_i[i], count, K_j + K_tj);
							gsl_matrix_from(regX_v, indep_ti[t * I + i], count, K_j + K_tj + K_i);
							gsl_matrix_from(regX_v, indep_ji[j * I + i], count, K_j + K_tj + K_i + K_ti);
							if (indep_tji.find((t * J + j) * I + i)!=indep_tji.end()) {
	    						gsl_matrix_from(regX_v, indep_tji[(t * J + j) * I + i], count, K_j + K_tj + K_i + K_ti + K_ji);
	                        }
	                        // finish debugging
	                        count++;
						}
					}
					#pragma omp critical(dataupdate)
					Bayesupdate_sigknown(regX_v, regy_v, A_v_temp, B_v_temp);
					gsl_vector_free (regy_v);
					gsl_matrix_free (regX_v);
				}
			}
			// draw
			Bayesdraw_sigknown(A_v_temp, B_v_temp, K_v, 1.0, phi_v, gsl_r);
			gsl_matrix_free (A_v_temp);
			gsl_vector_free (B_v_temp);
		}
				
		// (5) Store updates in save

		// phi
		gsl_vector_view phi_r_y_save_g = gsl_matrix_row(phi_r_y_save, g);
		gsl_vector_memcpy(&phi_r_y_save_g.vector, phi_r_y);
		gsl_vector_view phi_r_n_save_g = gsl_matrix_row(phi_r_n_save, g);
		gsl_vector_memcpy(&phi_r_n_save_g.vector, phi_r_n);
		gsl_vector_view phi_s_y_save_g = gsl_matrix_row(phi_s_y_save, g);
		gsl_vector_memcpy(&phi_s_y_save_g.vector, phi_s_y);
		gsl_vector_view phi_s_n_save_g = gsl_matrix_row(phi_s_n_save, g);
		gsl_vector_memcpy(&phi_s_n_save_g.vector, phi_s_n);
		gsl_vector_view phi_v_save_g = gsl_matrix_row(phi_v_save, g);
		gsl_vector_memcpy(&phi_v_save_g.vector, phi_v);

		// sigma2
		sigma2_r_y_save[g] = sigma2_r_y;
		sigma2_r_n_save[g] = sigma2_r_n;
		sigma2_s_y_save[g] = sigma2_s_y;
		sigma2_s_n_save[g] = sigma2_s_n;  
				        
	}

	/* 4. Output Data */
	
	cout << "Starting Output..." << endl;

	// phi
	const char * f_r_y = "phi_r_y.txt";
	gsl_matrix_output(phi_r_y_save, f_r_y, G2, K_r_y);
	const char * f_r_n = "phi_r_n.txt";
	gsl_matrix_output(phi_r_n_save, f_r_n, G2, K_r_n);
	const char * f_s_y = "phi_s_y.txt";
	gsl_matrix_output(phi_s_y_save, f_s_y, G2, K_s_y);
	const char * f_s_n = "phi_s_n.txt";
	gsl_matrix_output(phi_s_n_save, f_s_n, G2, K_s_n);
	const char * f_v = "phi_v.txt";
	gsl_matrix_output(phi_v_save, f_v, G2, K_v);

	// sigma2
	const char * f2_r_y = "sigma2_r_y.txt";
	vector_output(sigma2_r_y_save, f2_r_y, G2);
	const char * f2_r_n = "sigma2_r_n.txt";
	vector_output(sigma2_r_n_save, f2_r_n, G2);
	const char * f2_s_y = "sigma2_s_y.txt";
	vector_output(sigma2_s_y_save, f2_s_y, G2);
	const char * f2_s_n = "sigma2_s_n.txt";
	vector_output(sigma2_s_n_save, f2_s_n, G2);

	// imputed values
	ofstream out;
	out.open("imputed_r.txt");
	for (int j = 0; j < J; j++) {
		for (int t = born[j]; t<= dead[j]; t++) {
			out << j << "\t" << t << "\t" << r[t * J + j] << endl;
		}
	}
	out.close();
	
	out.open("imputed_s.txt");
	for (int j = 0; j < J; j++) {
		for (int t = born[j]+1; t<= dead[j]; t++) {
			out << j << "\t" << t << "\t" << s[t * J + j] << endl;
		}
	}
	out.close();
	
	out.open("imputed_v.txt");
	for (int t = 1; t < T; t++) {
		for (int i = 0; i < I; i++) {
			for (int iter = 0; iter < SetJ_M[t].size(); iter++) {
				int j = SetJ_M[t][iter];
				out << t << "\t" << j << "\t" << i << "\t" << v[(t * J + j) * I + i] << endl;
			}
			for (int iter = 0; iter < SetJ_U[t].size(); iter++) {
				int j = SetJ_U[t][iter];
				out << t << "\t" << j << "\t" << i << "\t" << v[(t * J + j) * I + i] << endl;
			}
		}
	}
	out.close();
	
	/* 5. Clear */
	gsl_rng_free(gsl_r);

	// A
	gsl_matrix_free(A_r_y);
	gsl_matrix_free(A_r_n);
	gsl_matrix_free(A_s_y);
	gsl_matrix_free(A_s_n);
	gsl_matrix_free(A_v);

	// phi_save
	gsl_matrix_free(phi_r_y_save);
	gsl_matrix_free(phi_r_n_save);
	gsl_matrix_free(phi_s_y_save);
	gsl_matrix_free(phi_s_n_save);
	gsl_matrix_free(phi_v_save);

	// phi
	gsl_vector_free(phi_r_y);
	gsl_vector_free(phi_r_n);
	gsl_vector_free(phi_s_y);
	gsl_vector_free(phi_s_n);
	gsl_vector_free(phi_v);

	// mu
	gsl_vector_free(mu_r_y);
	gsl_vector_free(mu_r_n);
	gsl_vector_free(mu_s_y);
	gsl_vector_free(mu_s_n);
	gsl_vector_free(mu_v);
	
	return 0;	
	
}

