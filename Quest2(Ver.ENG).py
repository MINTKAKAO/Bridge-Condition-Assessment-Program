import sympy as sp
import numpy as np
import matplotlib
matplotlib.use('TkAgg')  # Set interactive backend
import matplotlib.pyplot as plt

# Define variables
x, y = sp.symbols('x y')
L = float(input("Enter the length L of the beam (m): "))  # Unit: meters (m)
b = float(input("Enter the width b of the beam's cross-section (m): "))  # Unit: meters (m)
h = float(input("Enter the height h of the beam's cross-section (m): "))  # Unit: meters (m)
I_original = b * h**3 / 12  # Original moment of inertia (Unit: m^4)

# Define reaction force variables
R_A, R_B = sp.symbols('R_A R_B')  # Unit: Newton (N)

epsilon = 1e-6  # Set tolerance

# Input point moment loads
num_point_moments = int(input("Enter the number of point moment loads: "))
moment_positions = []  # Position unit: meters (m)
moment_magnitudes = []  # Magnitude unit: Newton-meters (Nm)

for i in range(num_point_moments):
    pos = float(input(f"Enter the x-position of point moment load {i+1} (m, 0 <= x <= {L}): "))
    mag = float(input(f"Enter the magnitude of point moment load {i+1} (Nm, positive for counterclockwise): "))
    moment_positions.append(pos)
    moment_magnitudes.append(mag)

# Input point vertical loads (enter positive values for downward loads)
num_point_loads = int(input("Enter the number of point vertical loads: "))
point_positions = []  # Position unit: meters (m)
point_magnitudes = []  # Magnitude unit: Newton (N)

for i in range(num_point_loads):
    pos = float(input(f"Enter the x-position of point load {i+1} (m, 0 <= x <= {L}): "))
    mag = float(input(f"Enter the magnitude of point load {i+1} (N, positive for downward): "))
    point_positions.append(pos)
    point_magnitudes.append(mag)

# Input uniformly distributed loads (enter positive values for downward loads)
num_continuous_loads = int(input("Enter the number of uniformly distributed loads: "))
continuous_loads = []  # Position unit: meters (m), load magnitude unit: Newton/meter (N/m)

for i in range(num_continuous_loads):
    start_pos = float(input(f"Enter the start x-position of uniformly distributed load {i+1} (m, 0 <= x <= {L}): "))
    end_pos = float(input(f"Enter the end x-position of uniformly distributed load {i+1} (m, {start_pos} <= x <= {L}): "))
    magnitude = float(input(f"Enter the magnitude of uniformly distributed load {i+1} (N/m, positive for downward): "))
    continuous_loads.append((start_pos, end_pos, magnitude))

# Input linearly distributed loads (enter positive values for downward loads)
num_linear_loads = int(input("Enter the number of linearly distributed loads: "))
linear_loads = []  # Position unit: meters (m), start and end load magnitudes unit: Newton/meter (N/m)

for i in range(num_linear_loads):
    start_pos = float(input(f"Enter the start x-position of linear load {i+1} (m, 0 <= x <= {L}): "))
    end_pos = float(input(f"Enter the end x-position of linear load {i+1} (m, {start_pos} <= x <= {L}): "))
    start_mag = float(input(f"Enter the start magnitude of linear load {i+1} (N/m): "))
    end_mag = float(input(f"Enter the end magnitude of linear load {i+1} (N/m): "))
    linear_loads.append((start_pos, end_pos, start_mag, end_mag))

# Input crack information
num_cracks = int(input("Enter the number of cracks: "))
crack_positions = []  # Crack x-positions (Unit: m)
crack_depths = []     # Crack depths d (Unit: m, measured from the bottom)

for i in range(num_cracks):
    crack_x = float(input(f"Enter the x-position of crack {i+1} (m): "))
    crack_d = float(input(f"Enter the y-direction depth d of crack {i+1} (m, measured from the bottom): "))
    crack_positions.append(crack_x)
    crack_depths.append(crack_d)

# Define shear force function (following Hibbeler's sign convention)
V = R_A  # Shear force starts with reaction R_A (Unit: N)
M = sp.Integer(0)  # Initialize bending moment as symbolic zero (Unit: Nm)

# Calculate shear force due to point vertical loads
for pos, mag in zip(point_positions, point_magnitudes):
    V += sp.Piecewise((-mag, x >= pos), (0, x < pos))

# Calculate shear force due to uniformly distributed loads
for start_pos, end_pos, mag in continuous_loads:
    V += sp.Piecewise(
        (0, x < start_pos),
        (-mag * (x - start_pos), (x >= start_pos) & (x <= end_pos)),
        (-mag * (end_pos - start_pos), x > end_pos)
    )

# Calculate shear force due to linearly distributed loads
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

# Calculate bending moment due to point moment loads
for pos, mag in zip(moment_positions, moment_magnitudes):
    M += sp.Piecewise((-mag, x >= pos), (0, x < pos))

# Integrate shear force to get bending moment
C1 = sp.symbols('C1')
M += sp.integrate(V, x) + C1  # Add integration constant C1 (Unit: Nm)

# Use boundary conditions to solve for constants and reaction forces
boundary_conditions = []

# Moment boundary condition (moment is zero at the left end)
boundary_conditions.append(M.subs(x, 0))

# Total load calculation (positive for downward loads)
total_point_load = sum(point_magnitudes)  # Unit: N
total_cont_load = sum(mag * (end_pos - start_pos) for start_pos, end_pos, mag in continuous_loads)  # Unit: N
total_linear_load = sum(0.5 * (start_mag + end_mag) * (end_pos - start_pos)
                        for start_pos, end_pos, start_mag, end_mag in linear_loads)  # Unit: N
total_load = total_point_load + total_cont_load + total_linear_load  # Unit: N

# Force equilibrium equation (positive upwards)
force_eq = sp.Eq(R_A + R_B - total_load, 0)

# Moment equilibrium equation (sum of moments about the left end)
moment_total = sp.Integer(0)

# Moment due to point loads
for pos, mag in zip(point_positions, point_magnitudes):
    moment_total += mag * pos

# Moment due to uniformly distributed loads
for start_pos, end_pos, mag in continuous_loads:
    moment_total += mag * (end_pos - start_pos) * (start_pos + (end_pos - start_pos)/2)

# Moment due to linearly distributed loads
for start_pos, end_pos, start_mag, end_mag in linear_loads:
    w_avg = (start_mag + end_mag) / 2
    moment_total += w_avg * (end_pos - start_pos) * (start_pos + (end_pos - start_pos)/2)

# **Subtract external moments from total moment**
for pos, mag in zip(moment_positions, moment_magnitudes):
    moment_total -= mag  # Subtract external moments with opposite sign

# Moment equilibrium equation
moment_eq = sp.Eq(R_B * L - moment_total, 0)

# Solve the equations
solutions = sp.solve([force_eq] + boundary_conditions + [moment_eq], (R_A, R_B, C1))
R_A_sol = solutions[R_A]  # Unit: N
R_B_sol = solutions[R_B]  # Unit: N
C1_sol = solutions[C1]    # Unit: Nm

# Substitute R_A, R_B, C1 into V and M
V = V.subs(R_A, R_A_sol)  # Unit: N
M = M.subs(R_A, R_A_sol).subs(C1, C1_sol)  # Unit: Nm

# Print results
print("\nLeft reaction force R_A (Unit: N):", float(R_A_sol))
print("Right reaction force R_B (Unit: N):", float(R_B_sol))

# Create list of discontinuities (excluding crack positions)
discontinuities = sorted(
    set(
        moment_positions +
        point_positions +
        [pos for load in continuous_loads for pos in (load[0], load[1])] +
        [pos for load in linear_loads for pos in (load[0], load[1])]
    )
)

# **Modify discontinuities list to exclude beam start (x=0) and end (x=L) (apply tolerance)**
discontinuities = [pos for pos in discontinuities if not (abs(pos) < epsilon or abs(pos - L) < epsilon)]

# Function to check if a point is a discontinuity
def is_discontinuity(x_val, discontinuities, epsilon=1e-6):
    return any(abs(x_val - pos) < epsilon for pos in discontinuities)

# Function to get crack depth at a specific x position
def get_crack_depth_at_x(x_val, crack_positions, crack_depths, epsilon=1e-6):
    for pos, d in zip(crack_positions, crack_depths):
        if abs(x_val - pos) < epsilon:
            return d
    return None

# Function to calculate first moment of area Q(y) from neutral axis to y position (Unit: m^3)
def Q(y_val, b, h, d=None):
    """
    y_val: y position to calculate (Unit: m)
    b: Cross-sectional width (Unit: m)
    h: Cross-sectional height (Unit: m)
    d: Crack depth (Unit: m), None if no crack
    """
    if d is None:
        # No crack
        y_neutral = 0  # Neutral axis at center
        y_top = h / 2
        y_bottom = -h / 2
    else:
        # With crack
        y_bottom = -h / 2 + d
        y_top = h / 2
        y_neutral = (y_bottom + y_top) / 2  # Neutral axis of effective cross-section

        if y_val < y_bottom:
            # Inside crack, so Q(y) = 0
            return 0

    y_prime = y_val - y_neutral  # Adjusted y value
    h_effective = y_top - y_bottom

    return (b / 2) * ((h_effective / 2)**2 - y_prime**2)

# Function to calculate effective moment of inertia I
def I_effective(b, h, d=None):
    if d is None:
        # No crack
        return b * h**3 / 12
    else:
        # With crack, calculate using effective height
        h_effective = h - d
        return b * h_effective**3 / 12

# Function to calculate neutral axis position
def calculate_neutral_axis(h, d=None):
    if d is None:
        return 0  # Neutral axis at center
    else:
        y_bottom = -h / 2 + d
        y_top = h / 2
        # Calculate neutral axis of effective cross-section
        return (y_bottom + y_top) / 2

# Evaluate shear force and moment numerically
x_vals = np.linspace(0, L, 500)
V_vals = np.array([V.subs(x, val).evalf() for val in x_vals], dtype=float)
M_vals = np.array([M.subs(x, val).evalf() for val in x_vals], dtype=float)

# Plot graphs
plt.figure(figsize=(10, 4))
plt.plot(x_vals, V_vals, label='Shear Force V(x)')
plt.title('Shear Force Diagram (SFD)')
plt.xlabel('x (m)')
plt.ylabel('Shear Force V(x) (N)')
plt.grid(True)
plt.legend()
plt.show()  # Display in blocking mode

plt.figure(figsize=(10, 4))
plt.plot(x_vals, M_vals, label='Moment M(x)', color='orange')
plt.title('Bending Moment Diagram (BMD)')
plt.xlabel('x (m)')
plt.ylabel('Moment M(x) (Nm)')
plt.grid(True)
plt.legend()
plt.show()  # Display in blocking mode

# Start stress calculation loop
while True:  # Start repeat loop
    try:
        # Calculate stress components at a specific position
        x_val = float(input(f"\nEnter the x-position to calculate stress (0 <= x <= {L} m): "))
        if not (0 <= x_val <= L):
            print(f"x must be between 0 and {L}.")
            continue
        y_val = float(input(f"Enter the y-position to calculate stress (m, y=0 at cross-section center, -{h/2} <= y <= {h/2} m): "))
        if not (-h/2 <= y_val <= h/2):
            print(f"y must be between -{h/2} and {h/2}.")
            continue
        z_val = float(input("Enter the z-position to calculate stress (m): "))

        # Get crack depth at x position
        crack_d = get_crack_depth_at_x(x_val, crack_positions, crack_depths)

        # Check if y position is inside crack
        if crack_d is not None:
            y_bottom = -h / 2 + crack_d
            if y_val < y_bottom:
                # Inside crack, so stress tensor is zero
                sigma_xx = 0
                sigma_xy = 0
                sigma_xz = 0
                sigma_yy = 0
                sigma_zz = 0

                stress_tensor = sp.Matrix([
                    [sigma_xx, sigma_xy, sigma_xz],
                    [sigma_xy, sigma_yy, 0],
                    [sigma_xz, 0, sigma_zz]
                ])

                print(f"\nAt position ({x_val} m, {y_val} m, {z_val} m), inside crack. Stress tensor is zero:")
                sp.pprint(stress_tensor)
                continue  # Move to next iteration

        # Calculate moment of inertia at the position
        I = I_effective(b, h, crack_d)

        # Calculate neutral axis position
        y_neutral = calculate_neutral_axis(h, crack_d)

        # Adjust y value from neutral axis
        y_corrected = y_val - y_neutral

        # Calculate Q(y) (Unit: m^3)
        Q_at_y = Q(y_val, b, h, crack_d)

        if is_discontinuity(x_val, discontinuities, epsilon):
            # At discontinuity, calculate left and right limits
            x_left = x_val - epsilon
            x_right = x_val + epsilon

            # Shear force at positions (Unit: N)
            V_at_x_left = float(V.subs(x, x_left).evalf())
            V_at_x_right = float(V.subs(x, x_right).evalf())

            # Moment at positions (Unit: Nm)
            M_at_x_left = float(M.subs(x, x_left).evalf())
            M_at_x_right = float(M.subs(x, x_right).evalf())

            # Calculate normal stress σ_xx (Unit: N/m^2)
            sigma_xx_left = -M_at_x_left * y_corrected / I
            sigma_xx_right = -M_at_x_right * y_corrected / I

            # Calculate shear stress σ_xy (Unit: N/m^2)
            sigma_xy_left = -V_at_x_left * Q_at_y / (I * b)
            sigma_xy_right = -V_at_x_right * Q_at_y / (I * b)

            # Construct stress tensor
            sigma_xz = 0  # Ignore shear stress in z-direction (Unit: N/m^2)
            sigma_yy = 0  # Unit: N/m^2
            sigma_zz = 0  # Unit: N/m^2

            # Modify stress tensor calculation
            # Choose the side with larger absolute value for each stress component
            if abs(sigma_xx_left) > abs(sigma_xx_right):
                sigma_xx = sigma_xx_left
            else:
                sigma_xx = sigma_xx_right

            if abs(sigma_xy_left) > abs(sigma_xy_right):
                sigma_xy = sigma_xy_left
            else:
                sigma_xy = sigma_xy_right

            # Construct stress tensor
            stress_tensor = sp.Matrix([
                [sigma_xx, sigma_xy, sigma_xz],
                [sigma_xy, sigma_yy, 0],
                [sigma_xz, 0, sigma_zz]
            ])

            print(f"\nStress tensor at discontinuity {x_val} m (Unit: N/m²):")
            sp.pprint(stress_tensor)

        else:
            # Not at discontinuity
            # Shear force at position (Unit: N)
            V_at_x = float(V.subs(x, x_val).evalf())

            # Moment at position (Unit: Nm)
            M_at_x = float(M.subs(x, x_val).evalf())

            # Calculate normal stress σ_xx (Unit: N/m^2)
            sigma_xx = -M_at_x * y_corrected / I

            # Calculate shear stress σ_xy (Unit: N/m^2)
            sigma_xy = -V_at_x * Q_at_y / (I * b)

            # Construct stress tensor
            sigma_xz = 0  # Ignore shear stress in z-direction (Unit: N/m^2)
            sigma_yy = 0  # Unit: N/m^2
            sigma_zz = 0  # Unit: N/m^2

            stress_tensor = sp.Matrix([
                [sigma_xx, sigma_xy, sigma_xz],
                [sigma_xy, sigma_yy, 0],
                [sigma_xz, 0, sigma_zz]
            ])

            print(f"\nStress tensor at position ({x_val} m, {y_val} m, {z_val} m) (Unit: N/m²):")
            sp.pprint(stress_tensor)

        # **Added section: Plot stress distribution**
        plot_stress = input("\nIf you want to see the stress distribution graph at this cross-section, enter 'p' (press any other key to skip): ").strip().lower()
        if plot_stress == 'p':
            # Calculate stress distribution at selected x position
            selected_x = x_val
            num_y_points = 100  # Number of points in y-direction
            y_distribution = np.linspace(-h/2, h/2, num_y_points)
            sigma_xx_distribution = []
            sigma_xy_distribution = []

            for y_pt in y_distribution:
                # Check for crack
                if crack_d is not None and y_pt < (-h/2 + crack_d):
                    # Inside crack, stress is zero
                    sigma_xx_distribution.append(0)
                    sigma_xy_distribution.append(0)
                    continue

                # Calculate effective moment of inertia I and neutral axis
                I_eff = I_effective(b, h, crack_d)
                y_neutral_axis = calculate_neutral_axis(h, crack_d)
                y_corr = y_pt - y_neutral_axis
                Q_y = Q(y_pt, b, h, crack_d)

                # Calculate moment and shear force
                M_at_selected_x = float(M.subs(x, selected_x).evalf())
                V_at_selected_x = float(V.subs(x, selected_x).evalf())

                # Calculate stresses
                sigma_xx_pt = -M_at_selected_x * y_corr / I_eff
                sigma_xy_pt = -V_at_selected_x * Q_y / (I_eff * b)

                sigma_xx_distribution.append(float(sigma_xx_pt))
                sigma_xy_distribution.append(float(sigma_xy_pt))

            # Plot graphs
            plt.figure(figsize=(8, 6))
            plt.plot(sigma_xx_distribution, y_distribution, label=r'$\sigma_{xx}$ (N/m²)')
            plt.plot(sigma_xy_distribution, y_distribution, label=r'$\sigma_{xy}$ (N/m²)')
            plt.title(f'Stress Distribution at x = {selected_x} m')
            plt.xlabel('Stress (N/m²)')
            plt.ylabel('y Position (m)')
            plt.legend()
            plt.grid(True)
            plt.show()

    except Exception as e:
        print(f"An error occurred: {e}")

    # Check if user wants to calculate again
    repeat = input("\nWould you like to calculate at another position? (y/n): ").strip().lower()
    if repeat != 'y':
        print("Ending calculations.")
        break
