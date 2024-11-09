import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

# Variable definition
x, y = sp.symbols('x y')
L = float(input("Enter the length of the beam L (m): "))  # Unit: meters (m)
b = float(input("Enter the width of the beam section b (m): "))  # Unit: meters (m)
h = float(input("Enter the height of the beam section h (m): "))  # Unit: meters (m)
I = b * h**3 / 12  # Second moment of area (unit: m^4)

# Support reaction variables
R_A, R_B = sp.symbols('R_A R_B')  # Unit: Newton (N)

# Point moment load input
num_point_moments = int(input("Enter the number of point moments: "))
moment_positions = []  # Position unit: meters (m)
moment_magnitudes = []  # Magnitude unit: Newton-meter (Nm)

for i in range(num_point_moments):
    pos = float(input(f"Enter the x-position for moment load {i+1} (m): "))
    mag = float(input(f"Enter the magnitude of moment load {i+1} (Nm, counterclockwise positive): "))
    moment_positions.append(pos)
    moment_magnitudes.append(mag)

# Point vertical load input (downward load is positive)
num_point_loads = int(input("Enter the number of point vertical loads: "))
point_positions = []  # Position unit: meters (m)
point_magnitudes = []  # Magnitude unit: Newton (N)

for i in range(num_point_loads):
    pos = float(input(f"Enter the x-position for vertical load {i+1} (m): "))
    mag = float(input(f"Enter the magnitude of vertical load {i+1} (N, downward positive): "))
    point_positions.append(pos)
    point_magnitudes.append(mag)

# Uniform distributed load input (downward load is positive)
num_continuous_loads = int(input("Enter the number of uniform distributed loads: "))
continuous_loads = []  # Position unit: meters (m), load magnitude unit: Newton/meter (N/m)

for i in range(num_continuous_loads):
    start_pos = float(input(f"Enter the start x-position for uniform load {i+1} (m): "))
    end_pos = float(input(f"Enter the end x-position for uniform load {i+1} (m): "))
    magnitude = float(input(f"Enter the magnitude of uniform load {i+1} (N/m, downward positive): "))
    continuous_loads.append((start_pos, end_pos, magnitude))

# Linear distributed load input (downward load is positive)
num_linear_loads = int(input("Enter the number of linear distributed loads: "))
linear_loads = []  # Position unit: meters (m), start and end load magnitude units: Newton/meter (N/m)

for i in range(num_linear_loads):
    start_pos = float(input(f"Enter the start x-position for linear load {i+1} (m): "))
    end_pos = float(input(f"Enter the end x-position for linear load {i+1} (m): "))
    start_mag = float(input(f"Enter the start magnitude of linear load {i+1} (N/m): "))
    end_mag = float(input(f"Enter the end magnitude of linear load {i+1} (N/m): "))
    linear_loads.append((start_pos, end_pos, start_mag, end_mag))

# Shear force function definition (according to Hibbeler's sign convention)
V = R_A  # Starting with support reaction R_A (unit: Newton, N)
M = sp.Integer(0)  # Initialize bending moment symbolically at zero (unit: Newton-meter, Nm)

# Shear force due to point loads
for pos, mag in zip(point_positions, point_magnitudes):
    V += sp.Piecewise((-mag, x >= pos), (0, x < pos))

# Shear force due to uniform distributed loads
for start_pos, end_pos, mag in continuous_loads:
    V += sp.Piecewise(
        (0, x < start_pos),
        (-mag * (x - start_pos), (x >= start_pos) & (x <= end_pos)),
        (-mag * (end_pos - start_pos), x > end_pos)
    )

# Shear force due to linear distributed loads
for start_pos, end_pos, start_mag, end_mag in linear_loads:
    # Slope of linear load
    slope = (end_mag - start_mag) / (end_pos - start_pos)  # Unit: N/m^2
    a = start_pos
    b_pos = end_pos
    w0 = start_mag
    V += sp.Piecewise(
        (0, x < a),
        (-(w0 * (x - a) + 0.5 * slope * (x - a)**2), (x >= a) & (x <= b_pos)),
        (-(w0 * (b_pos - a) + 0.5 * slope * (b_pos - a)**2), x > b_pos)
    )

# Moment due to point moments
for pos, mag in zip(moment_positions, moment_magnitudes):
    M += sp.Piecewise((-mag, x >= pos), (0, x < pos))

# Integrate shear force to find the moment
M += sp.integrate(V, x) + sp.Symbol('C1')  # Add integration constant C1 (unit: Nm)

# Use boundary conditions to solve for integration constant and reactions
boundary_conditions = [M.subs(x, 0) - 0]

# Calculate total load
total_point_load = sum(point_magnitudes)  # Unit: N
total_cont_load = sum(mag * (end_pos - start_pos) for start_pos, end_pos, mag in continuous_loads)  # Unit: N
total_linear_load = sum(0.5 * (start_mag + end_mag) * (end_pos - start_pos) for start_pos, end_pos, start_mag, end_mag in linear_loads)  # Unit: N
total_load = total_point_load + total_cont_load + total_linear_load  # Unit: N

# Force equilibrium equation
force_eq = sp.Eq(R_A + R_B, total_load)

# Moment equilibrium equation (sum of moments along the beam length)
moment_eq = sp.Eq(M.subs(x, L), 0)

# Solve equations
solutions = sp.solve([force_eq] + boundary_conditions + [moment_eq], (R_A, R_B, sp.Symbol('C1')))
R_A_sol = solutions[R_A]  # Unit: N
R_B_sol = solutions[R_B]  # Unit: N
C1_sol = solutions[sp.Symbol('C1')]  # Unit: Nm

# Substitute R_A, R_B, and C1 in V and M
V = V.subs(R_A, R_A_sol)  # Unit: N
M = M.subs(R_A, R_A_sol).subs('C1', C1_sol)  # Unit: Nm

# Stress distribution function (unit: Pa, N/m^2)
sigma = -M * y / I

# Output results
print("\nShear force function V(x) (unit: N):")
sp.pprint(V)
print("\nMoment function M(x) (unit: Nm):")
sp.pprint(M)
print("\nStress distribution function σ(x, y) (unit: N/m^2):")
sp.pprint(sigma)
print("\nLeft support reaction R_A (unit: N):", R_A_sol)
print("Right support reaction R_B (unit: N):", R_B_sol)

# Function to calculate first moment of area Q(y) (unit: m^3)
def Q(y_val):
    return b * (h**2 / 4 - y_val**2)

# Convert shear force function for numerical computation
V_func = sp.lambdify(x, V, 'numpy')

# Calculate stress components at a specific location
x_val = float(input("\nEnter the x-position to calculate stress (m): "))
y_val = float(input("Enter the y-position to calculate stress (m, y=0 at section center): "))
z_val = float(input("Enter the z-position to calculate stress (m): "))

# Small epsilon value for discontinuity handling
epsilon = 1e-6

# Calculate left and right x values around the specified position
x_left = x_val - epsilon
x_right = x_val + epsilon

# Calculate shear force at specified positions (unit: N)
V_at_x_left = V_func(x_left)
V_at_x_right = V_func(x_right)

# Calculate Q(y) (unit: m^3)
Q_at_y = Q(y_val)

# Calculate shear stress σ_xy (unit: N/m^2)
sigma_xy_left = V_at_x_left * Q_at_y / (I * b)
sigma_xy_right = V_at_x_right * Q_at_y / (I * b)

# Calculate normal stress σ_xx (unit: N/m^2)
sigma_xx_left = sigma.subs({x: x_left, y: y_val})
sigma_xx_left = float(sigma_xx_left)

sigma_xx_right = sigma.subs({x: x_right, y: y_val})
sigma_xx_right = float(sigma_xx_right)

# Stress tensor components
sigma_xz = 0  # Shear stress in z-direction (assumed zero, unit: N/m^2)
sigma_yy = 0  # Normal stress in y-direction (unit: N/m^2)
sigma_zz = 0  # Normal stress in z-direction (unit: N/m^2)

# Left stress tensor (unit: N/m^2)
stress_tensor_left = sp.Matrix([
    [sigma_xx_left, sigma_xy_left, sigma_xz],
    [sigma_xy_left, sigma_yy, 0],
    [sigma_xz, 0, sigma_zz]
])

# Right stress tensor (unit: N/m^2)
stress_tensor_right = sp.Matrix([
    [sigma_xx_right, sigma_xy_right, sigma_xz],
    [sigma_xy_right, sigma_yy, 0],
    [sigma_xz, 0, sigma_zz]
])

# Calculate norm (absolute sum of stress components)
norm_left = abs(sigma_xx_left) + abs(sigma_xy_left)
norm_right = abs(sigma_xx_right) + abs(sigma_xy_right)

# Choose larger stress tensor
if norm_left >= norm_right:
    stress_tensor = stress_tensor_left
    print(f"\nStress tensor at ({x_val}, {y_val}, {z_val}) (left value used, unit: N/m^2):")
else:
    stress_tensor = stress_tensor_right
    print(f"\nStress tensor at ({x_val}, {y_val}, {z_val}) (right value used, unit: N/m^2):")

# Print stress tensor
sp.pprint(stress_tensor)
