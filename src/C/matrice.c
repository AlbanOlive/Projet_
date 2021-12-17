#include <math.h>
#include <stdlib.h>
#include "matrice.h"

double norm_2(double *vect, int n)
{
    double S;
    int i;
    S = 0.0;
    for (i = 0; i < n; i++)
        S += *(vect + i) * *(vect + i);
    return sqrt(S);
}

double *produit(double *M, double *v, int n)
{
    double *P, S;
    int i, j;
    P = (double *)malloc(n * sizeof(double));
    for (i = 0; i < n; i++)
    {
        S = 0.0;
        for (j = 0; j < n; j++)
            S += *(M + i * n + j) * *(v + j);
        *(P + i) = S;
    }
    return P;
}

double *soustraction(double *v1, double *v2, int n)
{
    double *S;
    int i;
    S = (double *)malloc(n * sizeof(double));
    for (i = 0; i < n; i++)
        *(S + i) = *(v1 + i) - *(v2 + i);
    return S;
}

int estTridiagonale(double *M, int n)
{
    int i, j;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
        {
            if (i == j || i + 1 == j || j + 1 == i)
            {
                if (*(M + i * n + j) == 0)
                    return 0;
            }
            else if (*(M + i * n + j) != 0)
                return 0;
        }
    return 1;
}