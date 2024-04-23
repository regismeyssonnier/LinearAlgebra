# coding: latin
import numpy as np

# Exemple de données (n observations, p variables)
data = np.array([[1, 2, 3],
                 [4, 5, 6],
                 [7, 8, 9],
                 [10, 8, 12]])

# Calculer la moyenne de chaque variable
means = np.mean(data, axis=0)

# Centrer les données (soustraire la moyenne de chaque variable)
centered_data = data - means

# Calculer la matrice de covariance
cov_matrix = np.cov(centered_data, rowvar=False, bias=True)

# Afficher la matrice de covariance
print("Matrice de covariance :\n", cov_matrix)