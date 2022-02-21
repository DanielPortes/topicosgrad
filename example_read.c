/*
*   Matrix Market I/O example program
*
*   Read a real (non-complex) sparse matrix from a Matrix Market (v. 2.0) file.
*   and copies it to stdout.  This porgram does nothing useful, but
*   illustrates common usage of the Matrix Matrix I/O routines.
*   (See http://math.nist.gov/MatrixMarket for details.)
*
*   Usage:  a.out [filename] > output
*
*
*   NOTES:
*
*   1) Matrix Market files are always 1-based, i.e. the index of the first
*      element of a matrix is (1,1), not (0,0) as in C.  ADJUST THESE
*      OFFSETS ACCORDINGLY offsets accordingly when reading and writing
*      to files.
*
*   2) ANSI C requires one to use the "l" format modifier when reading
*      double precision floating point numbers in scanf() and
*      its variants.  For example, use "%lf", "%lg", or "%le"
*      when reading doubles, otherwise errors will occur.
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mmio.h"

#define epsilon 1.0e-7

/*
 *
*/

void read_general(int *IA, double *AA, int *JA, int *J, int *I, double *val, int N, int nz, FILE *f);

void read_sym_general(int *IA, double *AA, int *JA, int *J, int *I, double *val, int N, int nz, FILE *f);

void csr_matrix_vector(int *IA, double *AA, int *JA, int *x, double *y, int N);

void csr_matrix_symmetric_vector(int *IA, double *AA, int *JA, double *x, double *y, int N);

void csr_matrix_symmetric_vector_cond(int *IA, double *AA, int *JA, double *x, double *y, int N);

void grad_conj(int *IA, double *AA, int *JA, double *b, long kmax, int N);

void grad_conj_precond(int *IA, double *AA, int *JA, double *b, long kmax, int N);

int main(int argc, char *argv[])
{
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M, N, nz;
    int i, *I, *J;
    double *val;
    // CSR
    double *AA;
    int *JA, *IA;

    if (argc < 2)
    {
        fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
        exit(1);
    }
    else
    {
        if ((f = fopen(argv[1], "r")) == NULL)
        {
            exit(1);
        }
    }

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }


    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
        mm_is_sparse(matcode))
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) != 0)
    {
        exit(1);
    }

    /*
     * LEITURA PARA CSR
     * */
    I = (int *) malloc(nz * sizeof(int));
    J = (int *) malloc(nz * sizeof(int));
    val = (double *) malloc(nz * sizeof(double));
    AA = (double *) malloc(nz * sizeof(double));
    JA = (int *) malloc(nz * sizeof(int));
    IA = (int *) malloc((N + 1) * sizeof(int));

    read_general(IA, AA, JA, J, I, val, N, nz, f);

    /*
     * matrix vector
     * */
    double *x = (double *) malloc((N) * sizeof(double));
    double *y = (double *) malloc((N) * sizeof(double));
    for (int j = 0; j < N; ++j) // preenche vetor x
    {
        x[j] = 1.0;
    }


    csr_matrix_symmetric_vector(IA, AA, JA, x, y, N);

    printf("\ny: ");
    for (int j = 0; j < N; ++j)
    {
        printf("%f ", y[j]);
    }

     grad_conj(IA, AA, JA, y, 1000, N);
//    grad_conj_precond(IA, AA, JA, y, 1000, N);


    free(I);
    free(J);
    free(val);
    free(AA);
    free(JA);
    free(IA);
    free(x);
    free(y);

    if (f != stdin)
    { fclose(f); }

    return 0;

}

void read_general(int *IA, double *AA, int *JA, int *J, int *I, double *val, int N, int nz, FILE *f)
{

    for (int i = 0; i < nz; i++)
    {
        fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
        I[i]--;  /* adjust from 1-based to 0-based */
        J[i]--;
    }

    // init IA
    for (int j = 0; j < N + 1; ++j)
    {
        IA[j] = 0;
    }

    int mont = 0;
    IA[0] = mont;
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < nz; ++j)
        {
            if (I[j] == i)
            {
                AA[mont] = val[j];
                JA[mont] = J[j];
                mont++;
            }
        }
        IA[i + 1] = mont;
    }
}

// todo: nao terminado, mas read_general funciona bem pra matrizes simetricas
void read_sym_general(int *IA, double *AA, int *JA, int *J, int *I, double *val, int N, int nz, FILE *f)
{
    /* reseve memory for matrices */
    I = (int *) malloc((nz * 2) * sizeof(int));
    J = (int *) malloc((nz * 2) * sizeof(int));
    val = (double *) malloc((nz * 2) * sizeof(double));

    for (int i = 0; i < (nz * 2); i++)
    {
        fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
        I[i]--;  /* adjust from 1-based to 0-based */
        J[i]--;

        if (I[i] != J[i])
        {
            I[i + 1] = J[i];
            J[i + 1] = I[i];
            printf("%d ", I[i + 1]);
            printf("%d ", J[i + 1]);
            i++;
        }
    }
    if (f != stdin)
    { fclose(f); }

    // init IA
    AA = (double *) malloc((nz * 2) * sizeof(double));
    JA = (int *) malloc((nz * 2) * sizeof(int));
    IA = (int *) malloc((N + 1) * sizeof(int));

    // init IA
    for (int j = 0; j < N + 1; ++j)
    {
        IA[j] = 0;
    }

    int mont = 0;
    IA[0] = mont;
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < nz; ++j)
        {
            if (I[j] == i)
            {
                AA[mont] = val[j];
                JA[mont] = J[j];
                mont++;
            }
        }
        IA[i + 1] = mont;
    }
}

void csr_matrix_vector(int *IA, double *AA, int *JA, int *x, double *y, int N)
{
    for (int j = 0; j < N; ++j) // preenche vetor x
    {
        x[j] = 1;
    }

    for (int i = 0; i < N; ++i)
    {
        int jstart = IA[i];
        int jend = IA[i + 1];
        y[i] = 0.0;
        for (int k = jstart; k < jend; ++k)
        {
            int j = JA[k];
            y[i] += AA[k] * x[j];
        }
    }
}

void csr_matrix_symmetric_vector(int *IA, double *AA, int *JA, double *x, double *y, int N)
{
    for (int i = 0; i < N; ++i)
    {
        int jstart = IA[i];
        int jend = IA[i + 1];
        y[i] = 0.0;
        for (int k = jstart; k < jend; ++k)
        {
            int j = JA[k];
            y[i] += AA[k] * x[j];
            if (i != j)
            {
                y[j] += AA[k] * x[i];
            }
        }
    }
}

void grad_conj(int *IA, double *AA, int *JA, double *b, long kmax, int N)
{
    double *r = (double *) malloc((N) * sizeof(double));
    double *x = (double *) malloc((N) * sizeof(double));
    double *p = (double *) malloc((N) * sizeof(double));
    double *v = (double *) malloc((N) * sizeof(double));

    double rsold = 0.0, rsnew = 0.0, alpha = 0.0, beta = 0.0;
    long k = 0;
    for (int i = 0; i < N; i++)
    {
        x[i] = 0.0;              // aproximação inicial
        r[i] = b[i];             // r = b - A*x;
        p[i] = r[i];             // p = r;
        rsold += r[i] * r[i];    // rsold = r'*r;
    }
    double e = sqrt(rsold);
    while ((e > epsilon) && (k < kmax))
    {
        k++;
        // v = A*p;
        csr_matrix_symmetric_vector(IA, AA, JA, p, v, N);

        for (int i = 0; i < N; ++i)
        {
            printf("\n v : %f \n", v[i]);
        }

        double vp = 0.0;
        for (int i = 0; i < N; ++i)
        {
            vp += v[i] * p[i];       //    alpha = rsold/(v'*p);
        }
        alpha = rsold / vp;          //   alpha = rsold/(v'*p);
        printf("\n alpha : %f \n", alpha);

        rsnew = 0.0;
        for (int i = 0; i < N; ++i)
        {
            x[i] += alpha * p[i];   //  x     = x + alpha*p;
            printf(" xi: %f ", x[i]);
            r[i] -= alpha * v[i];   //  r     = r - alpha*v;
            rsnew += r[i] * r[i];   //  rsnew = r'*r;
        }
        e = sqrt(rsnew);
        printf("\n e: %f ", e);

        if (e < epsilon)
        {
            break;
        }
        beta = rsnew / rsold;
        for (int i = 0; i < N; ++i)
        {
            p[i] = r[i] + beta * p[i];
        }
        rsold = rsnew;
    }

    free(v);
    free(p);
    free(r);
}

void csr_matrix_symmetric_vector_cond(int *IA, double *AA, int *JA, double *x, double *y, int N)
{
    for (int i = 0; i < N; ++i)
    {
        int jstart = IA[i];
        int jend = IA[i + 1];
        y[i] = 0.0;
        for (int k = jstart; k < jend; ++k)
        {
            int j = JA[k];
            y[i] += (1 / AA[k]) * x[j];
            if (i != j)
            {
                y[j] += 0.0;
            }
        }
    }
}

void grad_conj_precond(int *IA, double *AA, int *JA, double *b, long kmax, int N)
{
    double *r = (double *) malloc((N) * sizeof(double));
    double *x = (double *) malloc((N) * sizeof(double));
    double *p = (double *) malloc((N) * sizeof(double));
    double *v = (double *) malloc((N) * sizeof(double));
    double *z = (double *) malloc((N) * sizeof(double));

    double rsold = 0.0, rsnew = 0.0, alpha = 0.0, beta = 0.0;
    long k = 0;
    for (int i = 0; i < N; i++)
    {
        x[i] = 0.0;              // aproximação inicial
        r[i] = b[i];             // r = b - A*x;
    }
    csr_matrix_symmetric_vector_cond(IA, AA, JA, r, z, N); // z = m^-1 * r
    for (int i = 0; i < N; i++)
    {
        p[i] = z[i];             // p = z;
        rsold += r[i] * z[i];    // rsold = r'*z;
    }
    double e = sqrt(rsold);
    while ((e > epsilon) && (k < kmax))
    {
        k++;
        // v = A*p;
        csr_matrix_symmetric_vector(IA, AA, JA, p, v, N);

        double vp = 0.0;
        for (int i = 0; i < N; ++i)
        {
            vp += v[i] * p[i];       //    alpha = rsold/(v'*p);
        }
        alpha = rsold / vp;          //   alpha = rsold/(v'*p);

        rsnew = 0.0;
        for (int i = 0; i < N; ++i)
        {
            x[i] += alpha * p[i];   //  x     = x + alpha*p;
            r[i] -= alpha * v[i];   //  r     = r - alpha*v;
            rsnew += r[i] * z[i];   //  rsnew = r'* z;
        }
        e = sqrt(rsnew);

        if (e < epsilon)
        {
            break;
        }
        csr_matrix_symmetric_vector_cond(IA, AA, JA, r, z, N); // z = m^-1 * r
        beta = rsnew / rsold;
        for (int i = 0; i < N; ++i)
        {
            p[i] = z[i] + beta * p[i];
        }
        rsold = rsnew;
    }

    free(v);
    free(z);
    free(p);
    free(r);
}

