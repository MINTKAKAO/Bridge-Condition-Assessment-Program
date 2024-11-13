import sympy as sp
import numpy as np
import matplotlib
matplotlib.use('TkAgg')  # Setting to interactive backend
import matplotlib.pyplot as plt

# Define variables
x, y = sp.symbols('x y')
L = float(input("Enter the length of the beam L (m): "))  # Unit: meters (m)
b = float(input("Enter the width of the beam cross-section b (m): "))  # Unit: meters (m)
h = float(input("Enter the height of the beam cross-section h (m): "))  # Unit: meters (m)
I = b * h**3 / 12  # Second moment of area (unit: m^4)

# Define support reaction variables
R_A, R_B = sp.symbols('R_A R_B')  # Unit: Newtons (N)

# Point moment load input
num_point_moments = int(input("Enter the number of point moments: "))
moment_positions = []  # Position unit: meters (m)
moment_magnitudes = []  # Magnitude unit: Newton-meters (Nm)

for i in range(num_point_moments):
    pos = float(input(f"Enter the x position of moment load {i+1} (m): "))
    mag = float(input(f"Enter the magnitude of moment load {i+1} (Nm, counterclockwise positive): "))
    moment_positions.append(pos)
    moment_magnitudes.append(mag)

# Point vertical load input (downward force is positive)
num_point_loads = int(input("Enter the number of point vertical loads: "))
point_positions = []  # Position unit: meters (m)
point_magnitudes = []  # Magnitude unit: Newtons (N)

for i in range(num_point_loads):
    pos = float(input(f"Enter the x position of vertical load {i+1} (m): "))
    mag = float(input(f"Enter the magnitude of vertical load {i+1} (N, downward positive): "))
    point_positions.append(pos)
    point_magnitudes.append(mag)

# Uniform distributed load input (downward force is positive)
num_continuous_loads = int(input("Enter the number of uniform distributed loads: "))
continuous_loads = []  # Position unit: meters (m), load magnitude unit: Newton/meter (N/m)

for i in range(num_continuous_loads):
    start_pos = float(input(f"Enter the start x position of uniform load {i+1} (m): "))
    end_pos = float(input(f"Enter the end x position of uniform load {i+1} (m): "))
    magnitude = float(input(f"Enter the magnitude of uniform load {i+1} (N/m, downward positive): "))
    continuous_loads.append((start_pos, end_pos, magnitude))

# Linear distributed load input (downward force is positive)
num_linear_loads = int(input("Enter the number of linear distributed loads: "))
linear_loads = []  # Position unit: meters (m), start and end load magnitude units: Newton/meter (N/m)

for i in range(num_linear_loads):
    start_pos = float(input(f"Enter the start x position of linear load {i+1} (m): "))
    end_pos = float(input(f"Enter the end x position of linear load {i+1} (m): "))
    start_mag = float(input(f"Enter the start magnitude of linear load {i+1} (N/m): "))
    end_mag = float(input(f"Enter the end magnitude of linear load {i+1} (N/m): "))
    linear_loads.append((start_pos, end_pos, start_mag, end_mag))

# Define shear force function (according to Hibbeler's sign convention)
V = R_A  # Start with reaction force R_A (unit: Newtons, N)
M = sp.Integer(0)  # Initialize bending moment symbolically at zero (unit: Newton-meter, Nm)

# Calculate shear force due to point loads
for pos, mag in zip(point_positions, point_magnitudes):
    V += sp.Piecewise((-mag, x >= pos), (0, x < pos))

# Calculate shear force due to uniform distributed loads
for start_pos, end_pos, mag in continuous_loads:
    V += sp.Piecewise(
        (0, x < start_pos),
        (-mag * (x - start_pos), (x >= start_pos) & (x <= end_pos)),
        (-mag * (end_pos - start_pos), x > end_pos)
    )

# Calculate shear force due to linear distributed loads
for start_pos, end_pos, start_mag, end_mag in linear_loads:
    slope = (end_mag - start_mag) / (end_pos - start_pos)  # Slope (N/m²)
    a = start_pos
    b_pos = end_pos
    w0 = start_mag
    V += sp.Piecewise(
        (0, x < a),
        (-(w0 * (x - a) + 0.5 * slope * (x - a)**2), (x >= a) & (x <= b_pos)),
        (-(w0 * (b_pos - a) + 0.5 * slope * (b_pos - a)**2), x > b_pos)
    )

# Calculate moment due to point moments
for pos, mag in zip(moment_positions, moment_magnitudes):
    M += sp.Piecewise((-mag, x >= pos), (0, x < pos))

# Integrate shear force to calculate the moment
C1 = sp.symbols('C1')
M += sp.integrate(V, x) + C1  # Add integration constant C1 (unit: Nm)

# Solve for constants using boundary conditions and reaction forces
boundary_conditions = [M.subs(x, 0) - 0]

# Total load calculations
total_point_load = sum(point_magnitudes)  # Unit: N
total_cont_load = sum(mag * (end_pos - start_pos) for start_pos, end_pos, mag in continuous_loads)  # Unit: N
total_linear_load = sum(0.5 * (start_mag + end_mag) * (end_pos - start_pos) for start_pos, end_pos, start_mag, end_mag in linear_loads)  # Unit: N
total_load = total_point_load + total_cont_load + total_linear_load  # Unit: N

# Equilibrium equations
force_eq = sp.Eq(R_A + R_B, total_load)
moment_eq = sp.Eq(M.subs(x, L), 0)

# Solve equations
solutions = sp.solve([force_eq] + boundary_conditions + [moment_eq], (R_A, R_B, C1))
R_A_sol = solutions[R_A]  # Unit: N
R_B_sol = solutions[R_B]  # Unit: N
C1_sol = solutions[C1]    # Unit: Nm

# Substitute R_A, R_B, and C1 into V and M
V = V.subs(R_A, R_A_sol)  # Unit: N
M = M.subs(R_A, R_A_sol).subs(C1, C1_sol)  # Unit: Nm

# Stress distribution function (unit: Pa, N/m^2)
sigma = -M * y / I

# Output reactions
print("\nReaction force at left support R_A (unit: N):", R_A_sol)
print("Reaction force at right support R_B (unit: N):", R_B_sol)

# First moment of area Q(y) calculation function (unit: m^3)
def Q(y_val, b, h):
    return (b / 2) * ((h**2 / 4) - y_val**2)

# Numerical evaluation of shear force and moment
x_vals = np.linspace(0, L, 500)
V_vals = np.array([V.subs(x, val).evalf() for val in x_vals], dtype=float)
M_vals = np.array([M.subs(x, val).evalf() for val in x_vals], dtype=float)

# Plot diagrams
plt.figure(figsize=(10, 4))
plt.plot(x_vals, V_vals, label='Shear Force V(x)')
plt.title('Shear Force Diagram (SFD)')
plt.xlabel('x (m)')
plt.ylabel('Shear Force V(x) (N)')
plt.grid(True)
plt.legend()
plt.show()  # Blocking mode

plt.figure(figsize=(10, 4))
plt.plot(x_vals, M_vals, label='Moment M(x)', color='orange')
plt.title('Bending Moment Diagram (BMD)')
plt.xlabel('x (m)')
plt.ylabel('Moment M(x) (Nm)')
plt.grid(True)
plt.legend()
plt.show()  # Blocking mode

# List of discontinuities
discontinuities = sorted(
    set(
        moment_positions +
        point_positions +
        [pos for load in continuous_loads for pos in (load[0], load[1])] +
        [pos for load in linear_loads for pos in (load[0], load[1])]
    )
)

# Check for discontinuity
def is_discontinuity(x_val, discontinuities, epsilon=1e-6):
    return any(abs(x_val - pos) < epsilon for pos in discontinuities)

# Stress calculation loop
while True:
    try:
        x_val = float(input(f"\nEnter x position to calculate stress (0 <= x <= {L} m): "))
        if not (0 <= x_val <= L):
            print(f"x position must be between 0 and {L}.")
            continue
        y_val = float(input(f"Enter y position to calculate stress (m, y=0 at the section center, -{h/2} <= y <= {h/2} m): "))
        if not (-h/2 <= y_val <= h/2):
            print(f"y position must be between -{h/2} and {h/2}.")
            continue
        z_val = float(input("Enter z position to calculate stress (m): "))

        Q_at_y = Q(y_val, b, h)

        if is_discontinuity(x_val, discontinuities):
            epsilon = 1e-6
            x_left = x_val - epsilon
            x_right = x_val + epsilon

            V_at_x_left = float(V.subs(x, x_left).evalf())
            V_at_x_right = float(V.subs(x, x_right).evalf())

            sigma_xy_left = V_at_x_left * Q_at_y / (I * b)
            sigma_xy_right = V_at_x_right * Q_at_y / (I * b)

            sigma_xx_left = float(sigma.subs({x: x_left, y: y_val}).evalf())
            sigma_xx_right = float(sigma.subs({x: x_right, y: y_val}).evalf())

            sigma_xz = 0  # Assume shear stress in z direction is zero (unit: N/m²)
            sigma_yy = 0  # Unit: N/m²
            sigma_zz = 0  # Unit: N/m²

            stress_tensor_left = sp.Matrix([
                [sigma_xx_left, sigma_xy_left, sigma_xz],
                [sigma_xy_left, sigma_yy, 0],
                [sigma_xz, 0, sigma_zz]
            ])

            stress_tensor_right = sp.Matrix([
                [sigma_xx_right, sigma_xy_right, sigma_xz],
                [sigma_xy_right, sigma_yy, 0],
                [sigma_xz, 0, sigma_zz]
            ])

            norm_left = abs(sigma_xx_left) + abs(sigma_xy_left)
            norm_right = abs(sigma_xx_right) + abs(sigma_xy_right)

            if norm_left >= norm_right:
                stress_tensor = stress_tensor_left
                print(f"\nLeft stress tensor at discontinuity {x_val} m (unit: N/m²):")
            else:
                stress_tensor = stress_tensor_right
                print(f"\nRight stress tensor at discontinuity {x_val} m (unit: N/m²):")

            sp.pprint(stress_tensor)

        else:
            V_at_x = float(V.subs(x, x_val).evalf())
            sigma_xx = float(sigma.subs({x: x_val, y: y_val}).evalf())
            sigma_xy = V_at_x * Q_at_y / (I * b)

            sigma_xz = 0
            sigma_yy = 0
            sigma_zz = 0

            stress_tensor = sp.Matrix([
                [sigma_xx, sigma_xy, sigma_xz],
                [sigma_xy, sigma_yy, 0],
                [sigma_xz, 0, sigma_zz]
            ])

            print(f"\nStress tensor at ({x_val} m, {y_val} m, {z_val} m) (unit: N/m²):")
            sp.pprint(stress_tensor)

    except Exception as e:
        print(f"Error occurred: {e}")

    repeat = input("\nWould you like to calculate stress at another position? (y/n): ").strip().lower()
    if repeat != 'y':
        print("Exiting the calculation.")
        break
