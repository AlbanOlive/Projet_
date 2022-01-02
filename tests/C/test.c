#include "../../src/C/jgs.h"
#include <stdio.h>

#define EPS 0.00001
#define ITMAX 1000

void aff(double *R, int nbLignes, int nbColonnes)
{
    int i, j;
    for (i = 0; i < nbLignes; i++)
    {
        for (j = 0; j < nbColonnes; j++)
            printf("%lf ", *(R + i * nbColonnes + j));
        printf("\n");
    }
}

int main()
{
    double x0[] = {0.0, 0.0, 0.0, 0.0};
    double A1[] = {5.0, -1.0, 1.0, 2.0, 8.0, -1.0, -1.0, 1.0, 4.0};
    double A2[] = {2.0, 8.0, -1.0, 5.0, -1.0, -1.0, -1.0, 1.0, 4.0};
    double A3[] = {1.0, 2.0, -2.0, 1.0, 1.0, 1.0, 2.0, 2.0, 1.0};
    double A4[] = {2.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0};
    double b1[] = {10.0, 11.0, 3.0};
    double b2[] = {2.0, 1.0, 1.0};
    double b3[] = {3.0, 0.0, 1.0};
    double b4[] = {4.0, 4.0, 4.0};

    printf("\nD1\n");
    aff(A1, 3, 3);
    printf("Jacobi : ");
    aff(jacobi(A1, b1, x0, EPS, ITMAX, 3), 1, 3);
    printf("Gauss-Seidel : ");
    aff(gauss_seidel(A1, b1, x0, EPS, ITMAX, 3), 1, 3);
    printf("Resultat théorique : 2 1 1\n");

    printf("\nD2\n");
    aff(A2, 3, 3);
    printf("Jacobi : ");
    aff(jacobi(A2, b2, x0, EPS, ITMAX, 3), 1, 3);
    printf("Gauss-Seidel : ");
    aff(gauss_seidel(A2, b2, x0, EPS, ITMAX, 3), 1, 3);
    printf("Resultat théorique : Divergence\n");

    printf("\nD3\n");
    aff(A3, 3, 3);
    printf("Jacobi : ");
    aff(jacobi(A3, b3, x0, EPS, ITMAX, 3), 1, 3);
    printf("Gauss-Seidel : ");
    aff(gauss_seidel(A3, b3, x0, EPS, ITMAX, 3), 1, 3);
    printf("Resultat théorique : Divergence\n");

    printf("D4\n");
    aff(A4, 3, 3);
    printf("Jacobi : ");
    aff(jacobi(A4, b4, x0, EPS, ITMAX, 3), 1, 3);
    printf("Gauss-Seidel : ");
    aff(gauss_seidel(A4, b4, x0, EPS, ITMAX, 3), 1, 3);
    printf("Resultat théorique : 1 1 1\n");

    return 0;
}