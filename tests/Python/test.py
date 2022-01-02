import numpy as np
import os,sys
sys.path.append(os.path.realpath(os.path.join(os.path.dirname(__file__), '../../src/Python/')))
from jacobi import *
from gauss_seidel import *


eps = 10**-5
itmax = 10**3
x0 = np.array([0, 0, 0]).T
A1 = np.array([[5, -1, 1], [2, 8, -1], [-1, 1, 4]])
A2 = np.array([[2, 8, -1], [5, -1, -1], [-1, 1, 4]])
A3 = np.array([[1, 2, -2], [1, 1, 1], [2, 2, 1]])
A4 = np.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]])
b1 = np.array([10, 11, 3]).T
b2 = np.array([2, 1, 1]).T
b3 = np.array([3, 0, 1]).T
b4 = np.array([4, 4, 4]).T

print("A1")
print(A1)
print("b1")
print(b1)
print("Jacobi :", jacobi(A1, b1, x0, eps, itmax)[0])
print("Gauss-Seidel :", gauss_seidel(A1, b1, x0, eps, itmax)[0])
print("Resultat théorique : [2, 1, 1]")

print("A2")
print(A2)
print("b2")
print(b2)
print("Jacobi :", jacobi(A2, b2, x0, eps, itmax)[0])
print("Gauss-Seidel :", gauss_seidel(A2, b2, x0, eps, itmax)[0])
print("Resultat théorique : Divergence")

print("A3")
print(A3)
print("b3")
print(b3)
print("Jacobi :", jacobi(A3, b3, x0, eps, itmax)[0])
print("Gauss-Seidel :", gauss_seidel(A3, b3, x0, eps, itmax)[0])
print("Resultat théorique : Jacobi : [1, 0, -1]")
print("Gauss-Seidel : Divergence")

print("A4")
print(A4)
print("b4")
print(b4)
print("Jacobi :", jacobi(A4, b4, x0, eps, itmax)[0])
print("Gauss-Seidel :", gauss_seidel(A4, b4, x0, eps, itmax)[0])
print("Resultat théorique : Jacobi : Divergence")
print("Gauss-Seidel : [1, 1, 1]")
