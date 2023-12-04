import sympy
from sympy import Matrix

E, m1, m2, S11, S12, S21, S22 = sympy.symbols("E m1 m2 S11 S12 S21 S22")

# Define the matrix
A = Matrix([[E-m1-S11, -S21], [-S12, E-m2-S22]])

# Compute the determinant of A
det_A = A.det()

# Print the determinant
print("The determinant of the matrix A is:", det_A)
