import numpy as np
from jacobi import *
from gauss_seidel import *

# On definie les paramètres
eps = float(input("Saisir epsilon : "))
itmax = int(input("Nombre d'iteration max : "))
N = int(input("Saisir N : "))
t0 = int(input("Saisir t0 : "))
T = int(input("Saisir T : "))
alpha = float(input("Saisir alpha : "))
beta = float(input("Saisir beta : "))
h = (T-t0) / (N+1)
v = h*h                                                  # v = h² pour alléger les calculs
x0 = np.zeros(N).T


# On créer le vecteur x des points d'évaluations
x = np.zeros(N)
for i in range(1,N+1) : 
    x[i-1] = t0 + i*h


# On créer la matrice A
A = 1/v * (np.eye(N) + np.eye(N) - np.diag(np.ones(N-1),1) - np.diag(np.ones(N-1),-1))

# On definie la fonction S selon le sujet
def S(u) :          
    y = 0
    return y

# On créer la matrice b
def mat_b (t0, T, N, alpha, beta) :       # α = u(t0)       β = u(T)
    b = [0.0]*(N)
    b[0] = S(x[1]) +  alpha / v           # b[0] = S_1 + α/h²   
    for i in range (1,N-1) : 
        b[i] =  S(x[i+1])                 # b[i] = S(t0 + h*(i+1)) = S(x_i+1) = S_i+1
    
    b[N-1] = S(T) + beta / v              # b[N-1] = S_N + β/h²
    return b

b = mat_b(t0, T, N, alpha, beta)

# Affichage des matrice A et b
print("A :")
print(A)
print("b :")
print(b)

# On retrouve u grâce aux méthodes de Jacobi et Gauss-Seidel
print("x =",x)

print("Jacobi : u(x) = ", jacobi(A, b, x0, eps, itmax)[0])

print("Gauss-Seidel : u(x) = ", gauss_seidel(A, b, x0, eps, itmax)[0])

