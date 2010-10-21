#ifndef NNLS_H
#define NNLS_H

extern "C" {

int nnls_(double *a, int *mda, int *m, int *n, 
	  double *b, double *x, double *rnorm, 
	  double *w, double *zz, int *index, int *mode);

#define NNLS nnls_

}

#endif

