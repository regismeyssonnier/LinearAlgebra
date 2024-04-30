# coding: latin
import numpy as np
import sympy as sp
import scipy as sc
import matplotlib.pyplot as plt
import random
from mpl_toolkits.mplot3d import Axes3D  # Importer les outils 3D

def norme5(s, t, x, y, z):
	return s**2 + t**2 + x**2 + y**2 + z**2

def SEP(vecteur, dim):

	# Définir les symboles
	s, t, x, y, z = sp.symbols('s t x y z')

	print(vecteur)

	dc = {s:0, t:1, x:2, y:3, z:4}
	sb = [s, t, x, y, z] 
	vis = [0] * 5
	eq = []
	ind_eq = 1
	ind_src = 0



	n = dim
	dsu = DisjointSetUnion(n)

	for id, v in enumerate(vecteur):
		group = []
		for s in v.free_symbols:
			group += [dc[s]]

		if len(group) > 1:
			num = group[0]
			for i in range(1, len(group)):
				dsu.union(num, group[i])
		elif len(group) == 1:
			stop = False
			for i in range(id, -1, -1):
				if i != id:
					for s in vecteur[i].free_symbols:
						if sb[dc[s]].compare(sb[group[0]]) == 0:
							print("add ", sb[dc[s]], sb[group[0]] )
							dsu.union(dc[s], group[0])
							stop = True
							break
					if stop: break

		print()

	print("solution")

	for i in range(0, 5):
		eq += [dsu.find(i)]
	print(eq)
	print(dsu.parent)

	seq = set()
	eqf = []
	for id, v in enumerate(vecteur):
		if len(v.free_symbols) == 0:
			eqf += [-1]
			continue

		for s in v.free_symbols:
			ds = dsu.find(dc[s]) 
			eqf += [ds]
			seq.add(ds)
			break

	print(eqf)
	print(seq)

	res = []

	for s in seq:
		i = 0
		dval={}
		for c in eqf:
			for sy in vecteur[i].free_symbols:
				if c == s:
					dval[sy] = 1
				else:
					dval[sy] = 0
			i+=1
		r = vecteur.subs(dval)
		res += [r]

	print( res )

	return res


class DisjointSetUnion:
	def __init__(self, n):
		self.parent = list(range(n))
		self.rank = [0] * n

	def find(self, u):
		if self.parent[u] != u:
			self.parent[u] = self.find(self.parent[u])  # Path compression
		return self.parent[u]

	def union(self, u, v):
		root_u = self.find(u)
		root_v = self.find(v)

		if root_u != root_v:
			# Union by rank
			if self.rank[root_u] > self.rank[root_v]:
				self.parent[root_v] = root_u
			elif self.rank[root_u] < self.rank[root_v]:
				self.parent[root_u] = root_v
			else:
				self.parent[root_v] = root_u
				self.rank[root_u] += 1

	def are_connected(self, u, v):
		return self.find(u) == self.find(v)





# Exemple de donnÃƒÂƒÃ‚Â©es (n observations, p variables)
#data = np.array([[4, 3, -2],
#                 [-3, -1, 3],
#                 [2, 3, 0]])

#data = np.array([[0, 3, -1],
#                 [2, -1, 1],
#                 [0, 0, 2]])

sz_MX = 5
symbols = 's t x y z'
data = np.array([[5, 0, 4, -2, -3],
			  [-2, 3, -3, 2, 4],
			  [0, 0, 3, 0, 0],
			  [0, 0, 0, 3, 1],
			  [1, 0, 2, -1, 1]])

#data = np.array([[5, 0, 4, -2, -3],
#			  [-2, 3, -3, 2, 4],
#			  [7, 8, 3, 6, 0],
#			  [1, 6, 5, 3, 1],
#			  [1, 10, 2, -1, 1]])

#data = np.random.randn(5, 5) * 10.0

#data = np.array([[4, 3, -2],
#                 [-3, -1, 3],
#                 [2, 3, 0]])

# data = []
# for i in range(0, 5):
# 	vd = []
# 	for j in range(0, 5):
# 		vd += [random.randint(0, 10)]
# 	data += [vd]

# data = np.array(data)


# data = np.array([[5, 1],
# 				 [0, 5]])

# Calculer les valeurs propres
eigenvalues, _ = np.linalg.eig(data)

eigenvalues = np.real(eigenvalues)
eigenvalues = np.round(eigenvalues, decimals=5)

# Tolerance pour la proximitÃ© des valeurs propres
tolerance = 1e-5

# Initialiser un dictionnaire pour stocker les valeurs propres uniques et leurs multiplicitÃ©s
eigenvalue_multiplicities = {}

# Parcourir toutes les valeurs propres calculÃ©es
for eigenvalue in eigenvalues:
	found = False
	
	# VÃ©rifier si cette valeur propre est dÃ©jÃ  dans le dictionnaire avec une tolÃ©rance
	for stored_eigenvalue in eigenvalue_multiplicities:
		if np.isclose(eigenvalue, stored_eigenvalue, atol=tolerance):
			eigenvalue_multiplicities[stored_eigenvalue] += 1
			found = True
			break
	
	# Si la valeur propre n'a pas Ã©tÃ© trouvÃ©e, l'ajouter au dictionnaire avec une multiplicitÃ© de 1
	if not found:
		eigenvalue_multiplicities[eigenvalue] = 1


eigenvectors = np.zeros((sz_MX, sz_MX))
ind_ev =0
P = np.zeros((sz_MX, sz_MX))
ind_P = 0

# Afficher les valeurs propres uniques avec leurs multiplicitÃ©s
for eigenvalue, multiplicity in eigenvalue_multiplicities.items():
	print(f"Valeur propre : {eigenvalue}, Multiplicité : {multiplicity}")

	# Séparer A en système d'équations A * x = b
	variables = sp.symbols(symbols)  # Définir les variables
	s, t, x, y, z = sp.symbols("s, t, x, y, z")
	_A = data - eigenvalue * np.identity(sz_MX)
	v = np.array([s, t, x, y, z])
	A = sp.Matrix(_A)
	print(_A)
	b = sp.Matrix([0]*sz_MX)
		
	Av = _A @ v
		
	system=(A, b)

	sl = sp.solve([Av[0], Av[1], Av[2], Av[3], Av[4]] , [s, t, x, y, z], dict=True)
	print(sl)

	#print("_________________________")
	#solution = sp.linsolve(system, variables)
	#print(solution)
	#print("_________________________")

	# Créer une liste de symboles
	symboles = [s, t, x, y, z]
	# Créer une liste de valeurs correspondantes avec une valeur par défaut
	valeurs = [sl[0].get(symbole, symbole) for symbole in symboles]
	print(valeurs)
	
	# Résoudre le système d'équations linéaires A * x = b
	#solution = sp.linsolve(valeurs , s, t, x, y, z)
	
	print("rank2 " , A.rank())
	# Afficher la solution
	print("Solution du système linéaire :")
	#print(type(solution))
	#print("s ", solution)
	print("ind_P ", ind_P)

	zero = True
	for i in range(sz_MX):
		if valeurs[i] != 0.0:
			zero = False
			break

	dim_A = sz_MX - A.rank()
	print("dim A ", dim_A)

	indpa = []
	if zero :
		b = sp.Matrix([0.00001]*sz_MX)
		system=(A, b)

		solution = sp.linsolve(system, variables)
		print(solution)
		sol = solution.args[0]

		ind = 0
		for s in sol:
			eigenvectors[ind, ind_ev] = s
			ind+= 1

		P[:, ind_P] = eigenvectors[:,ind_ev]
		ind_P += 1

		indpa += [ind_ev]
		ind_ev+=1

		print("EV ", eigenvectors)

	else:
		
		sol = sp.FiniteSet((valeurs[0], valeurs[1], valeurs[2], valeurs[3], valeurs[4])).args[0]
		print("sol ", sol)

		
	
		print("symbol present :", sol.free_symbols)
		print("symbol present :", sol.atoms())

		SVP = SEP(sol, sz_MX)

		

		for v in SVP:
			ind = 0
			for s in v:
				eigenvectors[ind, ind_ev] = s
				ind+= 1

			indpa += [ind_ev]
			ind_ev+=1
		print("EV ", eigenvectors)

	#dim restante

	if dim_A == multiplicity or (zero and multiplicity == 1) or multiplicity == 1 :continue

	already_copied = False
	if zero and multiplicity > 1:already_copied = True

	start_num = 1
	
	print("ind_ev ", ind_ev)
	for ind_VS in indpa:

		if not already_copied:
			P[:, ind_P] = eigenvectors[:,ind_VS]
			ind_P += 1
						
		ind_evc = 0

		variables = sp.symbols(symbols)  # Définir les variables
		b = []#sp.Matrix(eigenvectors[:,ind_VS])
		for indln in range(sz_MX):
			b += [eigenvectors[indln,ind_VS]]
		print("b= ", b)
		b = sp.Matrix(b)
					
		system=(A, b)
		#------------------

		# Résoudre le système d'équations linéaires A * x = b
		solution = sp.linsolve(system, variables)

		sol = list(solution)[0]
		print(solution)

		#start_num = 0
		dval = {}
		val = []
		for i in range(start_num, start_num+sz_MX):
			val += [i] 


		
		indv = 0
		for s in sol.free_symbols:
			dval[s] = val[indv]

			indv+=1

		start_num += indv

		print(dval)

		solf = sol.subs(dval)
		#n5 = norme5(solf[0], solf[1],solf[2], solf[3],solf[4])
		print(solf)

		ind = 0
		for s in solf:
			eigenvectors[ind, ind_ev] = s 
			ind+= 1

		print(ind_P, ind_ev)
		P[:, ind_P] = eigenvectors[:,ind_ev]
		ind_P += 1

		ind_evc = ind_ev
		ind_ev += 1
	
		print(eigenvectors)

		mult = 2
		while ind_ev < sz_MX and mult < multiplicity:
			print("ind_evc", ind_evc)
			variables = sp.symbols(symbols, real=False)  # Définir les variables
			b = []#sp.Matrix(eigenvectors[:,ind_VS])
			for indln in range(sz_MX):
				b += [eigenvectors[indln,ind_evc]]
			b = sp.Matrix(b)
			system = A, b

			print("b ", b)

			# Résoudre le système d'équations linéaires A * x = b
			solution = sp.linsolve(system, variables)
			print(solution)
			if solution is sp.EmptySet:
				break

			sol = list(solution)[0]
			
			#start_num = 0
			dval = {}
			val = []
			for i in range(start_num, start_num+sz_MX):
				val += [i] 


			indv = 0
			for s in sol.free_symbols:
				dval[s] = val[indv]

			indv+=1

			start_num += indv

			print(dval)

			solf = sol.subs(dval)
			#n5 = norme5(solf[0], solf[1],solf[2], solf[3],solf[4])
			print(solf)

			ind = 0
			for s in solf:
				eigenvectors[ind, ind_ev] = s
				ind+= 1

			P[:, ind_P] = eigenvectors[:,ind_ev]
			ind_P += 1

			ind_evc += 1
			ind_ev += 1

			mult +=1
							
			print(eigenvectors)

print ("v", eigenvalues) 
print("P", P)
P1 = np.linalg.inv(P)
print(P1)
J = P1 @ data @ P
print("La matrice de Jordan est : ")
print(np.round(J, decimals=2))
print("A-------------------------------------------")
#J[4, 0]=J[4, 1]=J[4, 2]=J[4, 3]=0 
AB = P @ J@P1
print(np.round(AB, decimals=0))
print(np.round(data, decimals=0))
print("---------------------------------------------")
eigenvectors1 = np.linalg.inv(eigenvectors)
print(eigenvectors1)
print(np.round(np.dot(np.dot(eigenvectors1, data), eigenvectors), decimals=4))
print("A3-----------------------------------")
print(np.round(data@data@data, decimals=4))
J = J@J@J
print("-A32")
print(np.round((P@J)@P1, decimals=4))
