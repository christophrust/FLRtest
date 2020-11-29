void internal_crossprod(double *x, int nrx, int ncx,
                        double *y, int nry, int ncy, double *z)
{
  double sum;
#define CROSSPROD_BODY                          \
  int NRX = nrx, NRY = nry, NCX = ncx;          \
  for (int i = 0; i < ncx; i++)                 \
    for (int k = 0; k < ncy; k++) {             \
	    sum = 0.0;                                \
	    for (int j = 0; j < nrx; j++)             \
        sum += x[j + i * NRX] * y[j + k * NRY];	\
	    z[i + k * NCX] = (double) sum;            \
    }
  CROSSPROD_BODY;
}
