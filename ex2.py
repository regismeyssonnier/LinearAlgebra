# coding: latin
import numpy as np

A = np.array([[0, 1],
			  [1, 0]])

eigenvalues, eigenvectors = np.linalg.eig(A)

v1 = eigenvectors[:, 0]
v2 = eigenvectors[:, 1]

print(eigenvalues)
print(v1)
print(v2)
print(v1 * eigenvalues[0])
print(v2 * eigenvalues[1])

print("--------------------------")
AX = np.dot(A,v1)
print(AX)
lambda_X = eigenvalues[0] * v1
print(lambda_X)

if np.allclose(AX, lambda_X):
    print("AX = lambda*X. Donc X est un vecteur propre avec valeur propre lambda.")
else:
    print("AX n'est pas égal à lambda*X. X n'est pas un vecteur propre avec valeur propre lambda.")