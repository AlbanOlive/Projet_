import numpy as np
from tridiag import estTridiagonale


def gauss_seidel(A, b, x0, eps, itmax):
    return gauss_seidel_tridiag(A, b, x0, eps, itmax) if estTridiagonale(A) else gauss_seidel_(A, b, x0, eps, itmax)


def gauss_seidel_tridiag(A, b, x0, eps, itmax):
    it = 0
    err = np.linalg.norm(b - A.dot(x0))
    x_old = x0
    n = len(A)
    x = np.zeros(n)
    while err > eps and it < itmax:
        S = A[0][1] * x_old[1]
        x[0] = (b[0] - S) / A[0][0]
        for i in range(1, n-1):
            S = A[i][i-1] * x[i-1] + A[i][i+1] * x_old[i+1]
            x[i] = (b[i] - S) / A[i][i]
        i += 1
        S = A[i][i-1] * x[i-1]
        x[i] = (b[i] - S) / A[i][i]
        x_old = x
        err = np.linalg.norm(b - A.dot(x))
        it += 1
    return x, it


def gauss_seidel_(A, b, x0, eps, itmax):
    it = 0
    err = np.linalg.norm(b - A.dot(x0))
    x_old = x0
    n = len(A)
    x = np.zeros(n)
    while err > eps and it < itmax:
        for i in range(n):
            S = 0
            for j in range(i):
                S += A[i][j] * x[j]
            for j in range(i + 1, n):
                S += A[i][j] * x_old[j]
            x[i] = (b[i] - S) / A[i][i]
        x_old = x
        err = np.linalg.norm(b - A.dot(x))
        it += 1
    return x, it
