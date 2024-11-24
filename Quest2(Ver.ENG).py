import sympy as sp
import numpy as np
import matplotlib
matplotlib.use('TkAgg')  # Enable interactive backend
import matplotlib.pyplot as plt

# Variable definitions
x, y = sp.symbols('x y')
L = float(input("Enter the length of the beam L (m): "))
b = float(input("Enter the width of the beam section b (m): "))
h = float(input("Enter the height of the beam section h (m): "))
I_original = b * h**3 / 12  # Original second moment of area (unit: m^4)

# Support reaction variables
R_A, R_B = sp.symbols('R_A R_B')  # Unit: Newton (N)

# Input for point moment loads
num_point_moments = int(input("Enter the number of point moment loads: "))
moment_positions = []  # Unit: meters (m)
moment_magnitudes = []  # Unit: Newton-meter (Nm)

for i in range(num_point_moments):
    pos = float(input(f"Enter the position x for moment load {i+1} (m): "))
    mag = float(input(f"Enter the magnitude of moment load {i+1} (Nm, counterclockwise positive): "))
    moment_positions.append(pos)
    moment_magnitudes.append(mag)

# Input for point vertical loads (downward loads are positive)
num_point_loads = int(input("Enter the number of point vertical loads: "))
point_positions = []  # Unit: meters (m)
point_magnitudes = []  # Unit: Newton (N)

for i in range(num_point_loads):
    pos = float(input(f"Enter the position x for vertical load {i+1} (m): "))
    mag = float(input(f"Enter the magnitude of vertical load {i+1} (N, downward positive): "))
    point_positions.append(pos)
    point_magnitudes.append(mag)

# Input for uniformly distributed loads (downward loads are positive)
num_continuous_loads = int(input("Enter the number of uniformly distributed loads: "))
continuous_loads = []  # Position unit: meters (m), load unit: Newton/meter (N/m)

for i in range(num_continuous_loads):
    start_pos = float(input(f"Enter the start position x of load {i+1} (m): "))
    end_pos = float(input(f"Enter the end position x of load {i+1} (m): "))
    magnitude = float(input(f"Enter the magnitude of load {i+1} (N/m, downward positive): "))
    continuous_loads.append((start_pos, end_pos, magnitude))

# Input for linearly varying loads (downward loads are positive)
num_linear_loads = int(input("Enter the number of linearly varying loads: "))
linear_loads = []  # Position unit: meters (m), load unit: Newton/meter (N/m)

for i in range(num_linear_loads):
    start_pos = float(input(f"Enter the start position x of load {i+1} (m): "))
    end_pos = float(input(f"Enter the end position x of load {i+1} (m): "))
    start_mag = float(input(f"Enter the starting magnitude of load {i+1} (N/m): "))
    end_mag = float(input(f"Enter the ending magnitude of load {i+1} (N/m): "))
    linear_loads.append((start_pos, end_pos, start_mag, end_mag))

# Input for crack information
num_cracks = int(input("Enter the number of cracks: "))
crack_positions = []  # Crack positions along x (unit: m)
crack_depths = []     # Crack depths (unit: m, from the bottom)

for i in range(num_cracks):
    crack_x = float(input(f"Enter the x position of crack {i+1} (m): "))
    crack_d = float(input(f"Enter the y-direction depth of crack {i+1} (m, from the bottom): "))
    crack_positions.append(crack_x)
    crack_depths.append(crack_d)

# Define shear force function (Hibbeler's sign convention)
V = R_A  # Shear force, starting with R_A (unit: Newton, N)
M = sp.Integer(0)  # Initialize bending moment as symbolic 0 (unit: Newton-meter, Nm)

# Shear force due to point vertical loads
for pos, mag in zip(point_positions, point_magnitudes):
    V += sp.Piecewise((-mag, x >= pos), (0, x < pos))

# Shear force due to uniformly distributed loads
for start_pos, end_pos, mag in continuous_loads:
    V += sp.Piecewise(
        (0, x < start_pos),
        (-mag * (x - start_pos), (x >= start_pos) & (x <= end_pos)),
        (-mag * (end_pos - start_pos), x > end_pos)
    )

# Shear force due to linearly varying loads
for start_pos, end_pos, start_mag, end_mag in linear_loads:
    slope = (end_mag - start_mag) / (end_pos - start_pos)  # Slope (N/m²)
    a = start_pos
    b_pos = end_pos
    w0 = start_mag
    V += sp.Piecewise(
        (0, x < a),
        (- (w0 * (x - a) + 0.5 * slope * (x - a)**2), (x >= a) & (x <= b_pos)),
        (- (w0 * (b_pos - a) + 0.5 * slope * (b_pos - a)**2), x > b_pos)
    )

# Bending moment calculation due to point moments
for pos, mag in zip(moment_positions, moment_magnitudes):
    M += sp.Piecewise((-mag, x >= pos), (0, x < pos))

# Integrate shear force to find bending moment
C1 = sp.symbols('C1')
M += sp.integrate(V, x) + C1  # Add integration constant C1 (unit: Nm)

# Apply boundary conditions to solve for constants
boundary_conditions = [M.subs(x, 0)]  # Moment is 0 at the left end
total_load = sum(point_magnitudes) + sum(mag * (end_pos - start_pos) for start_pos, end_pos, mag in continuous_loads)
force_eq = sp.Eq(R_A + R_B - total_load, 0)

# Solve equations
moment_eq = sp.Eq(M.subs(x, L), 0)  # Moment equilibrium at right end
solutions = sp.solve([force_eq] + boundary_conditions + [moment_eq], (R_A, R_B, C1))
R_A_sol = solutions[R_A]
R_B_sol = solutions[R_B]
C1_sol = solutions[C1]

# Substitute constants into V and M
V = V.subs([(R_A, R_A_sol), (R_B, R_B_sol)])
M = M.subs([(R_A, R_A_sol), (R_B, R_B_sol), (C1, C1_sol)])

# Evaluate shear force and moment for plotting
x_vals = np.linspace(0, L, 500)
V_vals = np.array([V.subs(x, val).evalf() for val in x_vals], dtype=float)
M_vals = np.array([M.subs(x, val).evalf() for val in x_vals], dtype=float)

# Plot shear force diagram
plt.figure(figsize=(10, 4))
plt.plot(x_vals, V_vals, label='Shear Force V(x)')
plt.title('Shear Force Diagram (SFD)')
plt.xlabel('x (m)')
plt.ylabel('Shear Force V(x) (N)')
plt.grid(True)
plt.legend()
plt.show()

# Plot bending moment diagram
plt.figure(figsize=(10, 4))
plt.plot(x_vals, M_vals, label='Bending Moment M(x)', color='orange')
plt.title('Bending Moment Diagram (BMD)')
plt.xlabel('x (m)')
plt.ylabel('Moment M(x) (Nm)')
plt.grid(True)
plt.legend()
plt.show()

# Calculate stresses
while True:
    x_val = float(input(f"\nEnter the x position to calculate stress (0 <= x <= {L} m): "))
    y_val = float(input(f"Enter the y position to calculate stress (m, -{h/2} <= y <= {h/2}): "))
    z_val = float(input("Enter the z position to calculate stress (m): "))

    crack_depth = next((d for p, d in zip(crack_positions, crack_depths) if abs(p - x_val) < 1e-6), None)
    if crack_depth and y_val < (-h/2 + crack_depth):
        print("Stress tensor is zero due to crack.")
        continue

    Q_y = b * (h**2 / 4 - y_val**2)  # First moment of area
    I = I_original if crack_depth is None else b * (h - crack_depth)**3 / 12
    V_x = float(V.subs(x, x_val).evalf())
    M_x = float(M.subs(x, x_val).evalf())
    sigma_xx = -M_x * y_val / I
    sigma_xy = V_x * Q_y / (I * b)
    stress_tensor = sp.Matrix([
        [sigma_xx, sigma_xy, 0],
        [sigma_xy, 0, 0],
        [0, 0, 0]
    ])
    print(f"\nStress tensor at ({x_val}, {y_val}, {z_val}) (N/m²):")
    sp.pprint(stress_tensor)

    repeat = input("\nDo you want to calculate another position? (y/n): ").strip().lower()
    if repeat != 'y':
        print("Exiting calculation.")
        break
