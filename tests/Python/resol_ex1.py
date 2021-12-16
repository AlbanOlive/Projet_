import numpy as np
from jacobi import *
from gauss_seidel import *

eps = float(input("Saisir epsilon : "))
itmax = int(input("Nombre d'iteration max : "))
N = int(input("Saisir N : "))
t0 = int(input("Saisir t0 : "))
T = int(input("Saisir T : "))
alpha = float(input("Saisir alpha : "))
beta = float(input("Saisir beta : "))
h = (T-t0)/N
v = h*h                                        # v = h² pour alleger les calculs

x0 = np.zeros(N).T
A = 1/v * (np.eye(N) + np.eye(N) - np.diag(np.ones(N-1),1) - np.diag(np.ones(N-1),-1))

def S(u) :          
    y = 2
    return y


def mat_b (t0, T, N, alpha, beta) :       # a = u(t0)       b = u(T)
    b = [0.0]*(N)
    h = (T - t0) / N
    b[0] = S(t0 + h) +  alpha / v      # b[0] = S_1 + a/h²   
    print (b[0])
    for i in range (1,N-1) : 
        b[i] =  S(t0 + h*(i+1))    # b[i] = S(t0 + h*(i+1)) = S(x_i+1) = S_i+1
        print(b[i])
    
    b[N-1] = S(T) + beta / v           # b[N-1] = S_N + b/h²
    
    return b

b = mat_b(t0, T, N, alpha, beta)

print("A :")
print(A)
print("b :")
print(b)
print("Jacobi :", jacobi(A, b, x0, eps, itmax))
print("Gauss-Seidel :", gauss_seidel(A, b, x0, eps, itmax))
