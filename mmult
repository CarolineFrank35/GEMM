#define __GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#ifndef SIZE
#define SIZE 1000
#endif

#ifndef SIZE_M
#define SIZE_M SIZE
#endif

#ifndef SIZE_N
#define SIZE_N SIZE
#endif

#ifndef SIZE_K
#define SIZE_K SIZE
#endif

#define long_long_t long long

static inline long_long_t timestamp();

//
// C = A*B
// (M x N) * (N x K) -> (M x K)
//
void mmult(double A[SIZE_M][SIZE_N],
           double B[SIZE_N][SIZE_K],
           double C[SIZE_M][SIZE_K])
{
    int jj, kk, ii, j, i, k, j0;
    const int mr = 4;
    const int nr = 4;
    const int mc = 100;
    const int kc = 250;
    const int nc = 500;
    
    double acc0, acc1, acc2, acc3;
    // Splitting Matrices
    for (jj = 0; jj < SIZE_K; jj += nc){
        
        for (kk = 0; kk < SIZE_N; kk += kc) {
            
            for (ii = 0; ii < SIZE_M; ii += mc) {
            
                for (j0 = jj; j0 < jj + nc; j0 += nr) {

                    // Kernel GEBP
                    for (j = j0; j < j0 + nr; ++j){
                        for (i = ii; i < ii + mc; i += mr) { // unrolled
                                // new row, restart accumulators
                                if (kk == 0) {
                                    acc0 = acc1 = acc2 = acc3 = 0;
                                } else {
                                    acc0 = C[i + 0][j];
                                    acc1 = C[i + 1][j];
                                    acc2 = C[i + 2][j];
                                    acc3 = C[i + 3][j];
                                }
                            
                            for (k = kk; k < kk + kc; ++k) {
                                acc0 += A[i + 0][k] * B[k][j];
                                acc1 += A[i + 1][k] * B[k][j];
                                acc2 += A[i + 2][k] * B[k][j];
                                acc3 += A[i + 3][k] * B[k][j];
                            }
                            C[i + 0][j] = acc0;
                            C[i + 1][j] = acc1;
                            C[i + 2][j] = acc2;
                            C[i + 3][j] = acc3;
                        }
                    }
                }
            }
        }
    }
}


double A[SIZE_M][SIZE_N];
double B[SIZE_N][SIZE_K];
double C[SIZE_M][SIZE_K];

int main(int argc, char* argv[])
{
    int         i, j;
    double      sum;
    /* Using double for nflop to prevent overflow on 32 bit archs */
    double      nflop;
    long_long_t tstart, tstop;
    double      tmmult;
    
    printf("Problem size: %d x %d\n", SIZE, SIZE);
    
    for (i = 0; i < SIZE_M; i++) {
        for (j = 0; j < SIZE_N; j++) {
            A[i][j] = (double)(i) + (double)(j);
        }
    }
    
    for (i = 0; i < SIZE_N; i++) {
        for (j = 0; j < SIZE_K; j++) {
            B[i][j] = (double)(i) + (double)(j);
        }
    }
    
    /* Two FLOP in inner loop: add and mul */
    nflop = 2.0 * (double)SIZE_M * (double)SIZE_N * (double)SIZE_K;
    
    tstart = timestamp();
    mmult(A, B, C);
    tstop  = timestamp();
    
    /* Duration in nanoseconds.
     * FLOP/ns = GFLOP/s
     */
    tmmult = (double)(tstop - tstart);
    
    /* Sum matrix elements as correctness hint */
    sum = 0.0;
    for (i = 0; i < SIZE_M && i < SIZE_K; i++) {
        sum += C[i][i];
    }
    printf("Trace mmult: %12.12g\n", sum);
    printf("M, N, K, tmmult_s, gflops_mmult\n");
    printf("%d, %d, %d, %f, %f \n",
           SIZE_M, SIZE_N, SIZE_K,
           tmmult, nflop / tmmult);
    
    return 0;
}

static inline long_long_t timestamp()
{
    struct timespec ts;
    long_long_t timestamp;
    clock_gettime(CLOCK_REALTIME, &ts);
    timestamp = ts.tv_sec * 1000000000LL + ts.tv_nsec;
    return timestamp;
}

