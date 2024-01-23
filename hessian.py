import numpy as np

# Given matrix
matrix_values = [
    [0.098537, 0.213269, -0.109665, -0.089204, -0.202469, 0.067879, 0,0,0],
    [0.213269, 0.512136, -0.070125, -0.195327, -0.448948, 0.134033, 0,0,0],
    [-0.109665, -0.070125, 0.712756, 0.092903, 0.196423, -0.123554,0,0,0 ],
    [-0.089204, -0.195327, 0.092903, 0.090416, 0.199241, -0.091000,0,0,0],
    [-0.202469, -0.448948, 0.196423, 0.199241, 0.443863, -0.185553,0,0,0],
    [0.067879, 0.134033, -0.123554, -0.091000, -0.185553, 0.144876, 0,0,0],
    [-0.009334, -0.017941, 0.016762, -0.001212, 0.003228, 0.023121,0.010546, 0,0],
    [-0.010799, -0.063188, -0.126297, -0.003913, 0.005086, 0.051519,0.014713, 0.058102, 0],
    [0.041787, -0.063909, -0.589202, -0.001904, -0.010870, -0.021322, -0.039883, 0.074778, 0.610524]
]
np.set_printoptions(formatter={'float': '{: 0.6f}'.format}, linewidth=150)

lower_triangular_matrix = np.tril(matrix_values)
symmetric_matrix = lower_triangular_matrix + lower_triangular_matrix.T - np.diag(lower_triangular_matrix.diagonal())
print(symmetric_matrix)
# print(np.array(matrix_values))
# # Convert the matrix to a NumPy array
matrix_array = np.array(symmetric_matrix)

# # Diagonalize the matrix
eigenvalues, eigenvectors = np.linalg.eigh(matrix_array)

# # Print the eigenvalues and eigenvectors
print("Eigenvalues:")
print(np.sqrt(abs(eigenvalues)))

print("\nEigenvectors:")
print(eigenvectors)
