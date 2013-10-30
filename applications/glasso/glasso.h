#ifndef GLASSO_H
#define GLASSO_H

extern "C" {
	/**
	 * This is a wrapper around the glasso fortran function
	 * implemented in the glasso R package:
	 *
	 * http://cran.r-project.org/web/packages/glasso/index.html
	 *
	 * All credit goes to the authors of the package
	 */
	void glasso_(int* n, float* ss, float* rho, int* ap_flag,
	             int* init_flag, int* trace_flag, int* dpen_flag,
	             float* thr, int* maxit, float* www, float* wwwi,
	             int* nniter, float* ddel, int* jerr);
}

#endif //GLASSO_H
