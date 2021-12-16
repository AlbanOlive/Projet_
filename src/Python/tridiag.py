import numpy as np


def estTridiagonale(A):
    nbLignes, nbColonnes = A.shape
    for i in range(nbLignes) :
        for j in range(nbColonnes) :
            if i == j or i + 1 == j or j + 1 == i :
                if A[i][j] == 0 :
                    return False
            elif A[i][j] != 0 :
                return False
    return True