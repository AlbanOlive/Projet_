#include <stdlib.h>
#include "matrice.h"
#include "jgs.h"

double *gauss_seidel_tridiag(double *A, double *b, double *x0, double eps, int itmax, int n)
{
    double *x_old, *x, err, S;
    int it, i;
    it = 0;
    err = norm_2(soustraction(b, produit(A, x0, n), n), n);
    x_old = x0;
    x = (double *)calloc(n, sizeof(double));
    while (err > eps && it < itmax)
    {
        S = *(A + 1) * *(x_old + 1);
        *x = (*b - S) / *A;
        for (i = 0; i < n - 1; i++)
        {
            S = *(A + i * n + i - 1) * *(x + i - 1) + *(A + i * n + i + 1) * *(x_old + i + 1);
            *(x + i) = (*(b + i) - S) / *(A + i * n + i);
        }
        S = *(A + i * n + i - 1) * *(x + i - 1);
        *(x + i) = (*(b + i) - S) / *(A + i * n + i);
        x_old = x;
        err = norm_2(soustraction(b, produit(A, x, n), n), n);
        it++;
    }
    return x;
}

double *gauss_seidel_(double *A, double *b, double *x0, double eps, int itmax, int n)
{
    double *x_old, *x, err, S;
    int it, i, j;
    it = 0;
    err = norm_2(soustraction(b, produit(A, x0, n), n), n);
    x_old = x0;
    x = (double *)calloc(n, sizeof(double));
    while (err > eps && it < itmax)
    {
        for (i = 0; i < n; i++)
        {
            S = 0;
            for (j = 0; j < i; j++)
                S += *(A + i * n + j) * *(x_old + j);
            for (j = i + 1; j < n; j++)
                S += *(A + i * n + j) * *(x_old + j);
            *(x + i) = (*(b + i) - S) / *(A + i * n + i);
        }
        x_old = x;
        err = norm_2(soustraction(b, produit(A, x, n), n), n);
        it++;
    }
    return x;
}

double *gauss_seidel(double *A, double *b, double *x0, double eps, int itmax, int n)
{
    return estTridiagonale(A, n) ? gauss_seidel_tridiag(A, b, x0, eps, itmax, n) : gauss_seidel_(A, b, x0, eps, itmax, n);
}
