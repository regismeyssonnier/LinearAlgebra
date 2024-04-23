# coding: latin
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Importer les outils 3D

# Exemple de données (n observations, p variables)
# data = np.array([[1, 2, 3],
#                  [4, 5, 6],
#                  [7, 8, 9],
#                  [10, 8, 12]])


# Définir les dimensions des données (n observations, p variables)
n = 100  # Nombre d'observations
p = 3    # Nombre de variables
data = np.random.randn(n, p)*20.0

# Calculer la moyenne de chaque variable
means = np.mean(data, axis=0)

# Centrer les données (soustraire la moyenne de chaque variable)
centered_data = data - means

# Calculer la matrice de covariance
cov_matrix = np.cov(centered_data, rowvar=False, bias=True)

# Afficher la matrice de covariance
print("Matrice de covariance :\n", cov_matrix)

eigenvalues, eigenvectors = np.linalg.eig(cov_matrix)

print("Eigen Values:\n", eigenvalues)
print("Eigen Vector:\n", eigenvectors)


# Créer une figure 3D pour afficher les données et les vecteurs propres
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Afficher les points de données d'origine
ax.scatter(data[:, 0], data[:, 2], data[:, 1], color='b', label='Donn\xe9es')

# Afficher les moyennes de chaque variable
ax.scatter(means[0], means[2], means[1], color='r', s=100, label='Moyennes', marker='o')

min_length = 1
max_length = 10
col = ['green','orange', 'black']
axis=['X','Y','Z']
# Afficher les vecteurs propres à partir des moyennes
num = [0,2,1]
for i in range(len(eigenvalues)):
    length = np.sqrt(eigenvalues[num[i]])
    scaled_length = min_length + (max_length - min_length) * (length - min(eigenvalues)) / (max(eigenvalues) - min(eigenvalues))
    ax.quiver(means[0], means[2], means[1],
              eigenvectors[0, i], eigenvectors[2, i], eigenvectors[1, i],
              length=scaled_length, normalize=True,
              color=col[i], label=axis[num[i]])


# Ajouter des étiquettes et un titre au graphique
ax.set_xlabel('Variable X')
ax.set_ylabel('Variable Z')
ax.set_zlabel('Variable Y')
ax.set_title('Donn\xe9es et vecteurs propres de la matrice de covariance')

# Afficher la légende
ax.legend()

# Afficher le graphique 3D
plt.show()