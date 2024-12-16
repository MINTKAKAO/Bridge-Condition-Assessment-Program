import sympy as sp
import numpy as np
import matplotlib
matplotlib.use('TkAgg')  # Set interactive backend
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

# Variable definitions
x, y = sp.symbols('x y')
L = float(input("Enter the beam length L (m): "))  # Unit: m
b = float(input("Enter the beam width b (m): "))  # Unit: m
h = float(input("Enter the beam height h (m): "))  # Unit: m
I_original = b * h**3 / 12  # Original second moment of area (Unit: m^4)

# Ultimate stress values input
ultimate_tensile_stress = float(input("Enter the ultimate tensile stress (N/m²): "))
ultimate_compressive_stress = float(input("Enter the ultimate compressive stress (N/m²): "))
ultimate_shear_stress = float(input("Enter the ultimate shear stress (N/m²): "))

# Reaction force variables
R_A, R_B = sp.symbols('R_A R_B')  # Unit: N
epsilon = 1e-6  # Tolerance

# Point moment loads input
num_point_moments = int(input("Enter the number of point moment loads: "))
moment_positions = []  # Unit: m
moment_magnitudes = []  # Unit: Nm

for i in range(num_point_moments):
    pos = float(input(f"Enter the position x of moment load {i+1} (m, 0 <= x <= {L}): "))
    mag = float(input(f"Enter the magnitude of moment load {i+1} (Nm, counter-clockwise positive): "))
    moment_positions.append(pos)
    moment_magnitudes.append(mag)

# Point vertical loads input (downward loads positive)
num_point_loads = int(input("Enter the number of point vertical loads: "))
point_positions = []  # Unit: m
point_magnitudes = []  # Unit: N

for i in range(num_point_loads):
    pos = float(input(f"Enter the position x of vertical load {i+1} (m, 0 <= x <= {L}): "))
    mag = float(input(f"Enter the magnitude of vertical load {i+1} (N, downward positive): "))
    point_positions.append(pos)
    point_magnitudes.append(mag)

# Uniform distributed loads input
num_continuous_loads = int(input("Enter the number of uniform distributed loads: "))
continuous_loads = []  # (start_pos, end_pos, magnitude), N/m

for i in range(num_continuous_loads):
    start_pos = float(input(f"Enter the start position x of uniform load {i+1} (m, 0 <= x <= {L}): "))
    end_pos = float(input(f"Enter the end position x of uniform load {i+1} (m, {start_pos} <= x <= {L}): "))
    magnitude = float(input(f"Enter the magnitude of uniform load {i+1} (N/m, downward positive): "))
    continuous_loads.append((start_pos, end_pos, magnitude))

# Linear distributed loads input
num_linear_loads = int(input("Enter the number of linear distributed loads: "))
linear_loads = []  # (start_pos, end_pos, start_mag, end_mag)

for i in range(num_linear_loads):
    start_pos = float(input(f"Enter the start position x of linear load {i+1} (m, 0 <= x <= {L}): "))
    end_pos = float(input(f"Enter the end position x of linear load {i+1} (m, {start_pos} <= x <= {L}): "))
    start_mag = float(input(f"Enter the start magnitude of linear load {i+1} (N/m): "))
    end_mag = float(input(f"Enter the end magnitude of linear load {i+1} (N/m): "))
    linear_loads.append((start_pos, end_pos, start_mag, end_mag))

# Crack information input
num_cracks = int(input("Enter the number of cracks: "))
crack_positions = []  # Unit: m
crack_depths = []     # Unit: m

for i in range(num_cracks):
    crack_x = float(input(f"Enter the x position of crack {i+1} (m): "))
    crack_d = float(input(f"Enter the depth of crack {i+1} (m, measured from bottom): "))
    crack_positions.append(crack_x)
    crack_depths.append(crack_d)

# Shear force function definition
V = R_A  # Shear force starts with R_A
M = sp.Integer(0)  # Initial moment is 0

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
    slope = (end_mag - start_mag) / (end_pos - start_pos)
    a = start_pos
    b_pos = end_pos
    w0 = start_mag
    V += sp.Piecewise(
        (0, x < a),
        (-(w0 * (x - a) + 0.5 * slope * (x - a)**2), (x >= a) & (x <= b_pos)),
        (-(w0 * (b_pos - a) + 0.5 * slope * (x - a)**2), x > b_pos)
    )

# Moment due to point moments
for pos, mag in zip(moment_positions, moment_magnitudes):
    M += sp.Piecewise((-mag, x >= pos), (0, x < pos))

# Integrate shear force to get bending moment
C1 = sp.symbols('C1')
M += sp.integrate(V, (x)) + C1

# Boundary conditions
boundary_conditions = []
boundary_conditions.append(M.subs(x, 0))  # At x=0, M=0

# Total loads
total_point_load = sum(point_magnitudes)
total_cont_load = sum(mag * (end_pos - start_pos) for start_pos, end_pos, mag in continuous_loads)
total_linear_load = sum(0.5*(start_mag+end_mag)*(end_pos - start_pos) for start_pos, end_pos, start_mag, end_mag)
total_load = total_point_load + total_cont_load + total_linear_load

# Force equilibrium
force_eq = sp.Eq(R_A + R_B - total_load, 0)

# Moment equilibrium
moment_total = sp.Integer(0)
for pos, mag in zip(point_positions, point_magnitudes):
    moment_total += mag * pos
for start_pos, end_pos, mag in continuous_loads:
    moment_total += mag*(end_pos - start_pos)*(start_pos+(end_pos - start_pos)/2)
for start_pos, end_pos, start_mag, end_mag in linear_loads:
    w_avg = (start_mag+end_mag)/2
    moment_total += w_avg*(end_pos - start_pos)*(start_pos+(end_pos - start_pos)/2)
for pos, mag in zip(moment_positions, moment_magnitudes):
    moment_total -= mag

moment_eq = sp.Eq(R_B * L - moment_total, 0)

solutions = sp.solve([force_eq]+boundary_conditions+[moment_eq], (R_A, R_B, C1))
R_A_sol = solutions[R_A]
R_B_sol = solutions[R_B]
C1_sol = solutions[C1]

V = V.subs(R_A, R_A_sol)
M = M.subs(R_A, R_A_sol).subs(C1, C1_sol)

print("\nLeft reaction force R_A (N):", float(R_A_sol))
print("Right reaction force R_B (N):", float(R_B_sol))

discontinuities = sorted(
    set(
        moment_positions +
        point_positions +
        [pos for load in continuous_loads for pos in (load[0], load[1])] +
        [pos for load in linear_loads for pos in (load[0], load[1])]
    )
)
discontinuities = [pos for pos in discontinuities if not (abs(pos) < epsilon or abs(pos - L) < epsilon)]

def is_discontinuity(x_val, discontinuities, epsilon=1e-6):
    return any(abs(x_val - pos) < epsilon for pos in discontinuities)

def get_crack_depth_at_x(x_val, crack_positions, crack_depths, epsilon=1e-6):
    for pos, d in zip(crack_positions, crack_depths):
        if abs(x_val - pos) < epsilon:
            return d
    return None

def Q(y_val, b, h, d=None):
    if d is None:
        y_neutral = 0
        y_top = h/2
        y_bottom = -h/2
    else:
        y_bottom = -h/2 + d
        y_top = h/2
        y_neutral = (y_bottom + y_top)/2
        if y_val < y_bottom:
            return 0

    y_prime = y_val - y_neutral
    h_effective = y_top - y_bottom

    return (b/2)*((h_effective/2)**2 - y_prime**2)

def I_effective(b, h, d=None):
    if d is None:
        return b*h**3/12
    else:
        h_effective = h - d
        return b*h_effective**3/12

def calculate_neutral_axis(h, d=None):
    if d is None:
        return 0
    else:
        y_bottom = -h/2 + d
        y_top = h/2
        return (y_bottom + y_top)/2

x_vals = np.linspace(0, L, 500)
V_vals = np.array([V.subs(x, val).evalf() for val in x_vals], dtype=float)
M_vals = np.array([M.subs(x, val).evalf() for val in x_vals], dtype=float)

plt.figure(figsize=(10,4))
plt.plot(x_vals, V_vals, label='Shear Force V(x)')
plt.title('Shear Force Diagram (SFD)')
plt.xlabel('x (m)')
plt.ylabel('Shear Force V(x) (N)')
plt.grid(True)
plt.legend()
plt.show()

plt.figure(figsize=(10,4))
plt.plot(x_vals, M_vals, label='Moment M(x)', color='orange')
plt.title('Bending Moment Diagram (BMD)')
plt.xlabel('x (m)')
plt.ylabel('Moment M(x) (Nm)')
plt.grid(True)
plt.legend()
plt.show()

        except Exception as e:
        print(f"Error occurred: {e}")

    repeat = input("\nWould you like to continue calculations at another location? (y/n): ").strip().lower()
    if repeat != 'y':
        print("Terminating the calculation process.")
        break

