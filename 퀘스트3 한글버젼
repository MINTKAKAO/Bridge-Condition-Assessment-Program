import sympy as sp
import numpy as np
import matplotlib
matplotlib.use('TkAgg')  # 대화형 백엔드로 설정
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

# 변수 정의
x, y = sp.symbols('x y')
L = float(input("보의 길이 L을 입력하세요 (m): "))  # 단위: m
b = float(input("보의 단면 폭 b을 입력하세요 (m): "))  # 단위: m
h = float(input("보의 단면 높이 h을 입력하세요 (m): "))  # 단위: m
I_original = b * h**3 / 12  # 원래 단면 2차 모멘트 (단위: m^4)

# 극한 응력값 입력
ultimate_tensile_stress = float(input("극한 인장응력을 입력하세요 (N/m²): "))
ultimate_compressive_stress = float(input("극한 압축응력을 입력하세요 (N/m²): "))
ultimate_shear_stress = float(input("극한 전단응력을 입력하세요 (N/m²): "))

# 지점 반력 변수 정의
R_A, R_B = sp.symbols('R_A R_B')  # 단위: N
epsilon = 1e-6  # 허용 오차 설정

# 포인트 모멘트 하중 입력
num_point_moments = int(input("포인트 모멘트 하중의 개수를 입력하세요: "))
moment_positions = []  # (단위: m)
moment_magnitudes = []  # (단위: Nm)

for i in range(num_point_moments):
    pos = float(input(f"{i+1}번째 모멘트 하중 위치 x 값을 입력하세요 (m, 0 <= x <= {L}): "))
    mag = float(input(f"{i+1}번째 모멘트 하중 크기 값을 입력하세요 (Nm, 시계 반대 방향 양수): "))
    moment_positions.append(pos)
    moment_magnitudes.append(mag)

# 포인트 수직하중 입력 (아래 방향 하중 양수)
num_point_loads = int(input("포인트 수직하중의 개수를 입력하세요: "))
point_positions = []  # (단위: m)
point_magnitudes = [] # (단위: N)

for i in range(num_point_loads):
    pos = float(input(f"{i+1}번째 포인트 하중 위치 x 값을 입력하세요 (m, 0 <= x <= {L}): "))
    mag = float(input(f"{i+1}번째 포인트 하중 크기 값을 입력하세요 (N, 아래 방향 양수): "))
    point_positions.append(pos)
    point_magnitudes.append(mag)

# 균등 분포하중 입력
num_continuous_loads = int(input("균등 분포하중의 개수를 입력하세요: "))
continuous_loads = []  # (start_pos, end_pos, magnitude), N/m

for i in range(num_continuous_loads):
    start_pos = float(input(f"{i+1}번째 균등 분포하중 시작 위치 x (m, 0 <= x <= {L}): "))
    end_pos = float(input(f"{i+1}번째 균등 분포하중 끝 위치 x (m, {start_pos} <= x <= {L}): "))
    magnitude = float(input(f"{i+1}번째 균등 분포하중 크기 (N/m, 아래 방향 양수): "))
    continuous_loads.append((start_pos, end_pos, magnitude))

# 선형 연속하중 입력
num_linear_loads = int(input("선형 연속하중의 개수를 입력하세요: "))
linear_loads = []  # (start_pos, end_pos, start_mag, end_mag)

for i in range(num_linear_loads):
    start_pos = float(input(f"{i+1}번째 선형 하중 시작 위치 x (m, 0 <= x <= {L}): "))
    end_pos = float(input(f"{i+1}번째 선형 하중 끝 위치 x (m, {start_pos} <= x <= {L}): "))
    start_mag = float(input(f"{i+1}번째 선형 하중 시작 크기 (N/m): "))
    end_mag = float(input(f"{i+1}번째 선형 하중 끝 크기 (N/m): "))
    linear_loads.append((start_pos, end_pos, start_mag, end_mag))

# 크랙 정보 입력
num_cracks = int(input("크랙의 개수를 입력하세요: "))
crack_positions = []  # (m)
crack_depths = []     # (m)

for i in range(num_cracks):
    crack_x = float(input(f"{i+1}번째 크랙 x 위치 (m): "))
    crack_d = float(input(f"{i+1}번째 크랙 깊이 d (m, 하부 기준): "))
    crack_positions.append(crack_x)
    crack_depths.append(crack_d)

# 전단력 함수 정의
V = R_A  # 전단력, R_A로 시작
M = sp.Integer(0)  # 굽힘모멘트 초기값 0

# 포인트 수직하중에 의한 전단력
for pos, mag in zip(point_positions, point_magnitudes):
    V += sp.Piecewise((-mag, x >= pos), (0, x < pos))

# 균등분포하중 전단력
for start_pos, end_pos, mag in continuous_loads:
    V += sp.Piecewise(
        (0, x < start_pos),
        (-mag * (x - start_pos), (x >= start_pos) & (x <= end_pos)),
        (-mag * (end_pos - start_pos), x > end_pos)
    )

# 선형하중 전단력
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

# 포인트 모멘트 하중에 의한 모멘트
for pos, mag in zip(moment_positions, moment_magnitudes):
    M += sp.Piecewise((-mag, x >= pos), (0, x < pos))

# 전단력 적분하여 모멘트
C1 = sp.symbols('C1')
M += sp.integrate(V, (x)) + C1

# 경계조건
boundary_conditions = []
boundary_conditions.append(M.subs(x, 0))  # x=0에서 M=0

# 전체 하중
total_point_load = sum(point_magnitudes)
total_cont_load = sum(mag * (end_pos - start_pos) for start_pos, end_pos, mag in continuous_loads)
total_linear_load = sum(0.5*(start_mag+end_mag)*(end_pos - start_pos) for start_pos, end_pos, start_mag, end_mag in linear_loads)
total_load = total_point_load + total_cont_load + total_linear_load

# 힘 평형
force_eq = sp.Eq(R_A + R_B - total_load, 0)

# 모멘트 평형
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

print("\n왼쪽 지점 반력 R_A (단위: N):", float(R_A_sol))
print("오른쪽 지점 반력 R_B (단위: N):", float(R_B_sol))

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

def check_failure_mohr(sigma_xx, sigma_xy, uts, ucs, uss, plot_mohr=False, failure_type="point"):
    # 주응력 계산
    sigma_1 = (sigma_xx/2.0) + np.sqrt((sigma_xx/2.0)**2 + sigma_xy**2)
    sigma_2 = (sigma_xx/2.0) - np.sqrt((sigma_xx/2.0)**2 + sigma_xy**2)

    # 최대 전단응력
    tau_max = np.sqrt((sigma_xx/2.0)**2 + sigma_xy**2)
    # 최대 인장응력
    sigma_max_tensile = max(sigma_1, sigma_2) if max(sigma_1, sigma_2) > 0 else 0
    # 최대 압축응력
    sigma_max_compressive = abs(min(sigma_1, sigma_2)) if min(sigma_1, sigma_2) < 0 else 0

    # failure_type에 따라 출력 조정
    # section일 경우 계산 과정을 보여주지 않는다.
    if failure_type == "point":
        print("\n--- 파괴 판정 계산 과정 ---")
        print(f"주응력 sigma_1: {sigma_1:.4e} N/m²")
        print(f"주응력 sigma_2: {sigma_2:.4e} N/m²")
        print(f"최대 전단응력 tau_max: {tau_max:.4e} N/m²")
        print(f"최대 인장응력 sigma_max_tensile: {sigma_max_tensile:.4e} N/m²")
        print(f"최대 압축응력 sigma_max_compressive: {sigma_max_compressive:.4e} N/m²")
        print(f"재료 강도: σ_ut={uts:.4e} N/m², σ_uc={ucs:.4e} N/m², τ_us={uss:.4e} N/m²")

    def plot_mohr_circle():
        # 여기부터 수정 시작
        fig, ax = plt.subplots()

        # Mohr's Criterion 내부 영역을 정의할 좌표
        # 인장 사각형과 압축 사각형을 결합한 다각형
        # 순서대로 인장 사각형, (-σ_uc,0) -> (0,σ_ut), 압축 사각형, (0,-σ_uc) -> (σ_ut,0)
        polygon_points = [
            (0, 0),
            (ultimate_tensile_stress, 0),
            (ultimate_tensile_stress, ultimate_tensile_stress),
            (0, ultimate_tensile_stress),
            (-ultimate_compressive_stress, 0),
            (-ultimate_compressive_stress, -ultimate_compressive_stress),
            (0, -ultimate_compressive_stress),
            (ultimate_tensile_stress, 0),
        ]

        # 내부 영역 색칠
        polygon = Polygon(polygon_points, closed=True, facecolor='green', edgecolor='black', alpha=0.3, label='Safe Region')
        ax.add_patch(polygon)

        # 경계선 그리기 (모든 선을 동일한 색상으로 설정)
        # 인장 사각형
        ax.plot([0, ultimate_tensile_stress, ultimate_tensile_stress, 0, 0],
                [0, 0, ultimate_tensile_stress, ultimate_tensile_stress, 0],
                color='black', linestyle='-', linewidth=2.0, label='Tensile Limit')

        # 압축 사각형
        ax.plot([0, -ultimate_compressive_stress, -ultimate_compressive_stress, 0, 0],
                [0, 0, -ultimate_compressive_stress, -ultimate_compressive_stress, 0],
                color='black', linestyle='-', linewidth=2.0, label='Compressive Limit')

        # (-σ_uc, 0)과 (0, σ_ut)을 잇는 직선
        ax.plot([-ultimate_compressive_stress, 0], [0, ultimate_tensile_stress],
                color='black', linestyle='--', linewidth=2.0, label='Failure Line 1')

        # (0, -σ_uc)과 (σ_ut, 0)을 잇는 직선
        ax.plot([0, ultimate_tensile_stress], [-ultimate_compressive_stress, 0],
                color='black', linestyle='--', linewidth=2.0, label='Failure Line 2')

        # 현재 응력 상태 점 표시
        ax.plot(sigma_1, sigma_2, 'ro', label='Current Stress State')

        ax.set_xlabel(r'$\sigma_1$ (N/m²)')
        ax.set_ylabel(r'$\sigma_2$ (N/m²)')
        ax.set_title("Mohr's Criterion (σ₁-σ₂ plane)")
        ax.grid(True)
        ax.legend()
        plt.show()
        # 여기까지 수정 끝

    # 파괴 판정 조건
    if sigma_max_tensile > uts:
        if failure_type == "point":
            print("판정: 최대 인장응력 > σ_ut → 파괴됨")
        if plot_mohr and failure_type == "point":
            plot_mohr_circle()
        return True, "Tensile failure"

    if tau_max > uss:
        if failure_type == "point":
            print("판정: τ_max > τ_us → 파괴됨")
        if plot_mohr and failure_type == "point":
            plot_mohr_circle()
        return True, "Shear failure"

    if sigma_max_compressive > ucs:
        if failure_type == "point":
            print("판정: 최대 압축응력 > σ_uc → 파괴됨")
        if plot_mohr and failure_type == "point":
            plot_mohr_circle()
        return True, "Compressive failure"

    if failure_type == "point":
        print("판정: 모든 조건 미충족 → 파괴되지 않음")
        if plot_mohr:
            plot_mohr_circle()

    return False, "No failure"

while True:
    try:
        x_val = float(input(f"\n응력을 계산할 x 위치를 입력하세요 (0 <= x <= {L} m): "))
        if not (0 <= x_val <= L):
            print(f"x는 0과 {L} 사이여야 합니다.")
            continue
        y_val = float(input(f"응력을 계산할 y 위치를 입력하세요 (m, 단면 중앙 y=0, -{h/2} <= y <= {h/2}): "))
        if not (-h/2 <= y_val <= h/2):
            print(f"y는 -{h/2}과 {h/2} 사이여야 합니다.")
            continue
        z_val = float(input("응력을 계산할 z 위치를 입력하세요 (m): "))

        crack_d = get_crack_depth_at_x(x_val, crack_positions, crack_depths)

        if crack_d is not None:
            y_bottom = -h/2 + crack_d
            if y_val < y_bottom:
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

                print(f"\n({x_val} m, {y_val} m, {z_val} m) 위치는 크랙 내부입니다. 응력 텐서는 0입니다:")
                sp.pprint(stress_tensor)
                continue

        I = I_effective(b, h, crack_d)
        y_neutral = calculate_neutral_axis(h, crack_d)
        y_corrected = y_val - y_neutral
        Q_at_y = Q(y_val, b, h, crack_d)

        if is_discontinuity(x_val, discontinuities, epsilon):
            x_left = x_val - epsilon
            x_right = x_val + epsilon

            V_at_x_left = float(V.subs(x, x_left).evalf())
            V_at_x_right = float(V.subs(x, x_right).evalf())
            M_at_x_left = float(M.subs(x, x_left).evalf())
            M_at_x_right = float(M.subs(x, x_right).evalf())

            sigma_xx_left = -M_at_x_left * y_corrected / I
            sigma_xx_right = -M_at_x_right * y_corrected / I
            sigma_xy_left = -V_at_x_left * Q_at_y / (I * b)
            sigma_xy_right = -V_at_x_right * Q_at_y / (I * b)

            if abs(sigma_xx_left) > abs(sigma_xx_right):
                sigma_xx = sigma_xx_left
            else:
                sigma_xx = sigma_xx_right

            if abs(sigma_xy_left) > abs(sigma_xy_right):
                sigma_xy = sigma_xy_left
            else:
                sigma_xy = sigma_xy_right

            sigma_xz = 0
            sigma_yy = 0
            sigma_zz = 0

            stress_tensor = sp.Matrix([
                [sigma_xx, sigma_xy, sigma_xz],
                [sigma_xy, sigma_yy, 0],
                [sigma_xz, 0, sigma_zz]
            ])

            print(f"\n불연속 지점 {x_val} m에서의 응력 텐서 (단위: N/m²):")
            sp.pprint(stress_tensor)

        else:
            V_at_x = float(V.subs(x, x_val).evalf())
            M_at_x = float(M.subs(x, x_val).evalf())

            sigma_xx = -M_at_x * y_corrected / I
            sigma_xy = -V_at_x * Q_at_y / (I * b)

            sigma_xz = 0
            sigma_yy = 0
            sigma_zz = 0

            stress_tensor = sp.Matrix([
                [sigma_xx, sigma_xy, sigma_xz],
                [sigma_xy, sigma_yy, 0],
                [sigma_xz, 0, sigma_zz]
            ])

            print(f"\n({x_val} m, {y_val} m, {z_val} m) 위치에서의 응력 텐서 (단위: N/m²):")
            sp.pprint(stress_tensor)

        # 부재(점) 파괴 판정
        failed, fail_mode = check_failure_mohr(
            sigma_xx,
            sigma_xy,
            ultimate_tensile_stress,
            ultimate_compressive_stress,
            ultimate_shear_stress,
            plot_mohr=True,
            failure_type="point"
        )
        if failed:
            print(f"파괴 발생(해당 점): {fail_mode}")
        else:
            print("파괴 없음(해당 점)")

        # 단면 전체 파괴 여부 판정 (계산 과정 숨김)
        num_y_samples = 50
        y_sample_points = np.linspace(-h/2, h/2, num_y_samples)
        cross_section_failed = True  # 모든 지점 파괴 가정 후, 하나라도 파괴 안되면 False

        # 여기는 failure_type="section" 이지만 어떤 출력도 하지 않도록 check_failure_mohr 내 로직 변경함
        for y_test in y_sample_points:
            if crack_d is not None:
                y_bottom = -h/2 + crack_d
                if y_test < y_bottom:
                    # 크랙 내부는 응력 0 → 파괴 아님
                    cross_section_failed = False
                    break

            I_eff = I_effective(b, h, crack_d)
            y_neutral_axis = calculate_neutral_axis(h, crack_d)
            y_corr = y_test - y_neutral_axis
            Q_y = Q(y_test, b, h, crack_d)

            V_at_x_test = float(V.subs(x, x_val).evalf())
            M_at_x_test = float(M.subs(x, x_val).evalf())

            sigma_xx_test = -M_at_x_test * y_corr / I_eff
            sigma_xy_test = -V_at_x_test * Q_y / (I_eff * b)

            # section 파괴판정: 결과만 판단
            test_failed, _ = check_failure_mohr(
                sigma_xx_test,
                sigma_xy_test,
                ultimate_tensile_stress,
                ultimate_compressive_stress,
                ultimate_shear_stress,
                plot_mohr=False,
                failure_type="section"
            )
            if not test_failed:
                cross_section_failed = False
                break

        if cross_section_failed:
            print(f"주의: x={x_val} m 단면에서 모든 지점 파괴 발생 → 단면 파괴")
        else:
            print(f"x={x_val} m 단면은 파괴되지 않음")

        # 단면 파괴 여부 판정 후 그래프 보고 싶은지 묻기
        plot_stress = input("\n이 지점(x={:.3f} m)에서의 단면 응력 분포 그래프를 보고 싶다면 'p'를 입력하세요 (그 외 키는 건너뜁니다): ".format(x_val)).strip().lower()
        if plot_stress == 'p':
            # 단면 응력분포 그래프 그리기
            selected_x = x_val
            num_y_points = 100
            y_distribution = np.linspace(-h/2, h/2, num_y_points)
            sigma_xx_distribution = []
            sigma_xy_distribution = []

            for y_pt in y_distribution:
                if crack_d is not None and y_pt < (-h/2 + crack_d):
                    sigma_xx_distribution.append(0)
                    sigma_xy_distribution.append(0)
                    continue

                I_eff = I_effective(b, h, crack_d)
                y_neutral_axis = calculate_neutral_axis(h, crack_d)
                y_corr = y_pt - y_neutral_axis
                Q_y = Q(y_pt, b, h, crack_d)

                M_at_selected_x = float(M.subs(x, selected_x).evalf())
                V_at_selected_x = float(V.subs(x, selected_x).evalf())

                sigma_xx_pt = -M_at_selected_x * y_corr / I_eff
                sigma_xy_pt = -V_at_selected_x * Q_y / (I_eff * b)

                sigma_xx_distribution.append(float(sigma_xx_pt))
                sigma_xy_distribution.append(float(sigma_xy_pt))

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
        print(f"오류 발생: {e}")

    repeat = input("\n다른 위치에서 계산을 진행하시겠습니까? (y/n): ").strip().lower()
    if repeat != 'y':
        print("계산을 종료합니다.")
        break
