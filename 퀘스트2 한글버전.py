import sympy as sp
import numpy as np
import matplotlib
matplotlib.use('TkAgg')  # 대화형 백엔드로 설정
import matplotlib.pyplot as plt

# 변수 정의
x, y = sp.symbols('x y')
L = float(input("보의 길이 L을 입력하세요 (m): "))  # 단위: 미터 (m)
b = float(input("보의 단면 폭 b을 입력하세요 (m): "))  # 단위: 미터 (m)
h = float(input("보의 단면 높이 h을 입력하세요 (m): "))  # 단위: 미터 (m)
I_original = b * h**3 / 12  # 원래 단면 2차 모멘트 (단위: m^4)

# 지점 반력 변수 정의
R_A, R_B = sp.symbols('R_A R_B')  # 단위: 뉴턴 (N)

epsilon = 1e-6  # 허용 오차 설정

# 포인트 모멘트 하중 입력
num_point_moments = int(input("포인트 모멘트 하중의 개수를 입력하세요: "))
moment_positions = []  # 위치 단위: 미터 (m)
moment_magnitudes = []  # 크기 단위: 뉴턴미터 (Nm)

for i in range(num_point_moments):
    pos = float(input(f"{i+1}번째 모멘트 하중 위치 x 값을 입력하세요 (m, 0 <= x <= {L}): "))
    mag = float(input(f"{i+1}번째 모멘트 하중 크기 값을 입력하세요 (Nm, 시계 반대 방향 양수): "))
    moment_positions.append(pos)
    moment_magnitudes.append(mag)

# 포인트 수직하중 입력 (아래 방향 하중은 양수로 입력)
num_point_loads = int(input("포인트 수직하중의 개수를 입력하세요: "))
point_positions = []  # 위치 단위: 미터 (m)
point_magnitudes = []  # 크기 단위: 뉴턴 (N)

for i in range(num_point_loads):
    pos = float(input(f"{i+1}번째 포인트 하중 위치 x 값을 입력하세요 (m, 0 <= x <= {L}): "))
    mag = float(input(f"{i+1}번째 포인트 하중 크기 값을 입력하세요 (N, 아래 방향 양수): "))
    point_positions.append(pos)
    point_magnitudes.append(mag)

# 연속 수직하중 입력 (아래 방향 하중은 양수로 입력)
num_continuous_loads = int(input("균등 분포하중의 개수를 입력하세요: "))
continuous_loads = []  # 위치 단위: 미터 (m), 하중 크기 단위: 뉴턴/미터 (N/m)

for i in range(num_continuous_loads):
    start_pos = float(input(f"{i+1}번째 균등 분포하중 시작 위치 x 값을 입력하세요 (m, 0 <= x <= {L}): "))
    end_pos = float(input(f"{i+1}번째 균등 분포하중 끝 위치 x 값을 입력하세요 (m, {start_pos} <= x <= {L}): "))
    magnitude = float(input(f"{i+1}번째 균등 분포하중 크기 값을 입력하세요 (N/m, 아래 방향 양수): "))
    continuous_loads.append((start_pos, end_pos, magnitude))

# 선형 연속하중 입력 (아래 방향 하중은 양수로 입력)
num_linear_loads = int(input("선형 연속하중의 개수를 입력하세요: "))
linear_loads = []  # 위치 단위: 미터 (m), 시작 및 끝 하중 크기 단위: 뉴턴/미터 (N/m)

for i in range(num_linear_loads):
    start_pos = float(input(f"{i+1}번째 선형 하중 시작 위치 x 값을 입력하세요 (m, 0 <= x <= {L}): "))
    end_pos = float(input(f"{i+1}번째 선형 하중 끝 위치 x 값을 입력하세요 (m, {start_pos} <= x <= {L}): "))
    start_mag = float(input(f"{i+1}번째 선형 하중 시작 크기 값을 입력하세요 (N/m): "))
    end_mag = float(input(f"{i+1}번째 선형 하중 끝 크기 값을 입력하세요 (N/m): "))
    linear_loads.append((start_pos, end_pos, start_mag, end_mag))

# 크랙 정보 입력
num_cracks = int(input("크랙의 개수를 입력하세요: "))
crack_positions = []  # 크랙 위치 x (단위: m)
crack_depths = []     # 크랙 깊이 d (단위: m, 하부에서부터의 거리)

for i in range(num_cracks):
    crack_x = float(input(f"{i+1}번째 크랙 위치 x 값을 입력하세요 (m): "))
    crack_d = float(input(f"{i+1}번째 크랙의 y 방향 깊이 d 값을 입력하세요 (m, 하부에서부터의 거리): "))
    crack_positions.append(crack_x)
    crack_depths.append(crack_d)

# 전단력 함수 정의 (히블러의 부호 규약에 맞게 설정)
V = R_A  # 전단력, 반력 R_A로 시작 (단위: N)
M = sp.Integer(0)  # 굽힘 모멘트를 심볼릭 0으로 초기화 (단위: Nm)

# 포인트 수직하중에 따른 전단력 계산
for pos, mag in zip(point_positions, point_magnitudes):
    V += sp.Piecewise((-mag, x >= pos), (0, x < pos))

# 균등 분포하중에 따른 전단력 계산
for start_pos, end_pos, mag in continuous_loads:
    V += sp.Piecewise(
        (0, x < start_pos),
        (-mag * (x - start_pos), (x >= start_pos) & (x <= end_pos)),
        (-mag * (end_pos - start_pos), x > end_pos)
    )

# 선형 연속하중에 따른 전단력 계산
for start_pos, end_pos, start_mag, end_mag in linear_loads:
    slope = (end_mag - start_mag) / (end_pos - start_pos)  # 기울기 (N/m²)
    a = start_pos
    b_pos = end_pos
    w0 = start_mag
    V += sp.Piecewise(
        (0, x < a),
        (- (w0 * (x - a) + 0.5 * slope * (x - a)**2), (x >= a) & (x <= b_pos)),
        (- (w0 * (b_pos - a) + 0.5 * slope * (b_pos - a)**2), x > b_pos)
    )

# 포인트 모멘트 하중에 따른 모멘트 계산
for pos, mag in zip(moment_positions, moment_magnitudes):
    M += sp.Piecewise((-mag, x >= pos), (0, x < pos))

# 전단력을 적분하여 모멘트 계산
C1 = sp.symbols('C1')
M += sp.integrate(V, x) + C1  # 적분 상수 C1 추가 (단위: Nm)

# 경계 조건을 이용하여 적분 상수 및 반력 계산
boundary_conditions = []

# 모멘트 경계 조건 (좌측 끝에서 모멘트는 0)
boundary_conditions.append(M.subs(x, 0))

# 전체 하중 계산 (아래 방향 하중은 양수, 위 방향은 음수로 입력)
total_point_load = sum(point_magnitudes)  # 단위: N
total_cont_load = sum(mag * (end_pos - start_pos) for start_pos, end_pos, mag in continuous_loads)  # 단위: N
total_linear_load = sum(0.5 * (start_mag + end_mag) * (end_pos - start_pos) for start_pos, end_pos, start_mag, end_mag in linear_loads)  # 단위: N
total_load = total_point_load + total_cont_load + total_linear_load  # 단위: N

# 힘의 평형 방정식 (위 방향을 양수로 취급)
force_eq = sp.Eq(R_A + R_B - total_load, 0)

# 모멘트 평형 방정식 (좌측 끝점에서 모멘트 합계)
moment_total = sp.Integer(0)

# 포인트 하중에 의한 모멘트
for pos, mag in zip(point_positions, point_magnitudes):
    moment_total += mag * pos

# 균등 분포하중에 의한 모멘트
for start_pos, end_pos, mag in continuous_loads:
    moment_total += mag * (end_pos - start_pos) * (start_pos + (end_pos - start_pos)/2)

# 선형 분포하중에 의한 모멘트
for start_pos, end_pos, start_mag, end_mag in linear_loads:
    w_avg = (start_mag + end_mag) / 2
    moment_total += w_avg * (end_pos - start_pos) * (start_pos + (end_pos - start_pos)/2)

# **모멘트 합계에서 외부 모멘트 빼기**
for pos, mag in zip(moment_positions, moment_magnitudes):
    moment_total -= mag  # 외부 모멘트는 반대 부호로 빼줌

# 모멘트 평형 방정식
moment_eq = sp.Eq(R_B * L - moment_total, 0)

# 방정식 풀기
solutions = sp.solve([force_eq] + boundary_conditions + [moment_eq], (R_A, R_B, C1))
R_A_sol = solutions[R_A]  # 단위: N
R_B_sol = solutions[R_B]  # 단위: N
C1_sol = solutions[C1]    # 단위: Nm

# R_A, R_B, C1 값을 V와 M에 대입
V = V.subs(R_A, R_A_sol)  # 단위: N
M = M.subs(R_A, R_A_sol).subs(C1, C1_sol)  # 단위: Nm

# 결과 출력
print("\n왼쪽 지점 반력 R_A (단위: N):", float(R_A_sol))
print("오른쪽 지점 반력 R_B (단위: N):", float(R_B_sol))

# 불연속점 목록 생성 (크랙 위치는 제외)
discontinuities = sorted(
    set(
        moment_positions +
        point_positions +
        [pos for load in continuous_loads for pos in (load[0], load[1])] +
        [pos for load in linear_loads for pos in (load[0], load[1])]
    )
)

# **불연속 지점 목록에서 보의 시작점(x=0)과 끝점(x=L)을 제외하도록 수정 (허용 오차 적용)**
discontinuities = [pos for pos in discontinuities if not (abs(pos) < epsilon or abs(pos - L) < epsilon)]

# 불연속점인지 확인하는 함수
def is_discontinuity(x_val, discontinuities, epsilon=1e-6):
    return any(abs(x_val - pos) < epsilon for pos in discontinuities)

# 특정 위치에서 크랙 깊이 가져오는 함수
def get_crack_depth_at_x(x_val, crack_positions, crack_depths, epsilon=1e-6):
    for pos, d in zip(crack_positions, crack_depths):
        if abs(x_val - pos) < epsilon:
            return d
    return None

# 중립축에서 y 위치까지의 1차 모멘트 Q(y) 계산 함수 (단위: m^3)
def Q(y_val, b, h, d=None):
    """
    y_val: 계산할 y 위치 (단위: m)
    b: 단면 폭 (단위: m)
    h: 단면 높이 (단위: m)
    d: 크랙 깊이 (단위: m), None이면 크랙 없음
    """
    if d is None:
        # 크랙이 없는 경우
        y_neutral = 0  # 중립축은 단면의 중앙
        y_top = h / 2
        y_bottom = -h / 2
    else:
        # 크랙이 있는 경우
        y_bottom = -h / 2 + d
        y_top = h / 2
        y_neutral = (y_bottom + y_top) / 2  # 유효 단면의 중립축 위치

        if y_val < y_bottom:
            # 크랙 내부이므로 Q(y) = 0
            return 0

    y_prime = y_val - y_neutral  # 수정된 y 값
    h_effective = y_top - y_bottom

    return (b / 2) * ((h_effective / 2)**2 - y_prime**2)

# 수정된 단면 2차 모멘트 I 계산 함수
def I_effective(b, h, d=None):
    if d is None:
        # 크랙이 없는 경우
        return b * h**3 / 12
    else:
        # 크랙이 있는 경우, 유효 높이로 계산
        h_effective = h - d
        return b * h_effective**3 / 12

# 중립축 위치 계산 함수 추가
def calculate_neutral_axis(h, d=None):
    if d is None:
        return 0  # 중립축은 단면의 중앙
    else:
        y_bottom = -h / 2 + d
        y_top = h / 2
        # 유효 단면의 중립축 위치 계산
        return (y_bottom + y_top) / 2

# 전단력과 모멘트를 수치적으로 평가
x_vals = np.linspace(0, L, 500)
V_vals = np.array([V.subs(x, val).evalf() for val in x_vals], dtype=float)
M_vals = np.array([M.subs(x, val).evalf() for val in x_vals], dtype=float)

# 그래프 그리기
plt.figure(figsize=(10, 4))
plt.plot(x_vals, V_vals, label='Shear Force V(x)')
plt.title('Shear Force Diagram (SFD)')
plt.xlabel('x (m)')
plt.ylabel('Shear Force V(x) (N)')
plt.grid(True)
plt.legend()
plt.show()  # 블로킹 모드로 표시

plt.figure(figsize=(10, 4))
plt.plot(x_vals, M_vals, label='Moment M(x)', color='orange')
plt.title('Bending Moment Diagram (BMD)')
plt.xlabel('x (m)')
plt.ylabel('Moment M(x) (Nm)')
plt.grid(True)
plt.legend()
plt.show()  # 블로킹 모드로 표시

# 응력 계산 루프 시작
while True:  # 반복 루프 시작
    try:
        # 특정 위치에서의 응력 성분을 계산
        x_val = float(input(f"\n응력을 계산할 x 위치를 입력하세요 (0 <= x <= {L} m): "))
        if not (0 <= x_val <= L):
            print(f"x 위치는 0과 {L} 사이여야 합니다.")
            continue
        y_val = float(input(f"응력을 계산할 y 위치를 입력하세요 (m, 단면 중앙이 y=0, -{h/2} <= y <= {h/2} m): "))
        if not (-h/2 <= y_val <= h/2):
            print(f"y 위치는 -{h/2}과 {h/2} 사이여야 합니다.")
            continue
        z_val = float(input("응력을 계산할 z 위치를 입력하세요 (m): "))

        # 해당 x 위치에서 크랙 깊이 가져오기
        crack_d = get_crack_depth_at_x(x_val, crack_positions, crack_depths)

        # y 위치가 크랙 내부에 있는지 확인
        if crack_d is not None:
            y_bottom = -h / 2 + crack_d
            if y_val < y_bottom:
                # 크랙 내부이므로 응력 텐서는 0
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
                continue  # 다음 반복으로 넘어감

        # 해당 위치에서 단면 2차 모멘트 I 계산
        I = I_effective(b, h, crack_d)

        # 중립축 위치 계산
        y_neutral = calculate_neutral_axis(h, crack_d)

        # 이동된 중립축으로부터의 y 값 계산
        y_corrected = y_val - y_neutral

        # Q(y) 계산 (단위: m^3)
        Q_at_y = Q(y_val, b, h, crack_d)

        if is_discontinuity(x_val, discontinuities, epsilon):
            # 불연속점인 경우 좌극한과 우극한 계산
            x_left = x_val - epsilon
            x_right = x_val + epsilon

            # 해당 위치에서의 전단력 계산 (단위: N)
            V_at_x_left = float(V.subs(x, x_left).evalf())
            V_at_x_right = float(V.subs(x, x_right).evalf())

            # 모멘트 계산 (단위: Nm)
            M_at_x_left = float(M.subs(x, x_left).evalf())
            M_at_x_right = float(M.subs(x, x_right).evalf())

            # 수직 응력 σ_xx 계산 (단위: N/m^2)
            sigma_xx_left = -M_at_x_left * y_corrected / I
            sigma_xx_right = -M_at_x_right * y_corrected / I

            # 전단 응력 σ_xy 계산 (단위: N/m^2)
            sigma_xy_left = -V_at_x_left * Q_at_y / (I * b)
            sigma_xy_right = -V_at_x_right * Q_at_y / (I * b)

            # 응력 텐서 구성
            sigma_xz = 0  # z 방향 전단 응력은 무시 (단위: N/m^2)
            sigma_yy = 0  # 단위: N/m^2
            sigma_zz = 0  # 단위: N/m^2

            # 응력 텐서 계산 부분 수정
            # 각 응력 성분별로 절댓값이 큰 쪽 선택
            if abs(sigma_xx_left) > abs(sigma_xx_right):
                sigma_xx = sigma_xx_left
            else:
                sigma_xx = sigma_xx_right

            if abs(sigma_xy_left) > abs(sigma_xy_right):
                sigma_xy = sigma_xy_left
            else:
                sigma_xy = sigma_xy_right

            # 응력 텐서 구성
            stress_tensor = sp.Matrix([
                [sigma_xx, sigma_xy, sigma_xz],
                [sigma_xy, sigma_yy, 0],
                [sigma_xz, 0, sigma_zz]
            ])

            print(f"\n불연속 지점 {x_val} m에서의 응력 텐서 (단위: N/m²):")
            sp.pprint(stress_tensor)

        else:
            # 불연속점이 아닌 경우
            # 해당 위치에서의 전단력 계산 (단위: N)
            V_at_x = float(V.subs(x, x_val).evalf())

            # 모멘트 계산 (단위: Nm)
            M_at_x = float(M.subs(x, x_val).evalf())

            # 수직 응력 σ_xx 계산 (단위: N/m^2)
            sigma_xx = -M_at_x * y_corrected / I

            # 전단 응력 σ_xy 계산 (단위: N/m^2)
            sigma_xy = -V_at_x * Q_at_y / (I * b)

            # 응력 텐서 구성
            sigma_xz = 0  # z 방향 전단 응력은 무시 (단위: N/m^2)
            sigma_yy = 0  # 단위: N/m^2
            sigma_zz = 0  # 단위: N/m^2

            stress_tensor = sp.Matrix([
                [sigma_xx, sigma_xy, sigma_xz],
                [sigma_xy, sigma_yy, 0],
                [sigma_xz, 0, sigma_zz]
            ])

            print(f"\n({x_val} m, {y_val} m, {z_val} m) 위치에서의 응력 텐서 (단위: N/m²):")
            sp.pprint(stress_tensor)

        # **추가된 부분: 응력 분포 그래프 그리기**
        plot_stress = input("\n특정 단면의 응력 분포를 그래프로 보고 싶으시면 'p'를 입력하세요 (그 외 입력 시 건너뜁니다): ").strip().lower()
        if plot_stress == 'p':
            # 선택된 x 위치에서의 응력 분포 계산
            selected_x = x_val
            num_y_points = 100  # y 방향으로 분할할 지점 수
            y_distribution = np.linspace(-h/2, h/2, num_y_points)
            sigma_xx_distribution = []
            sigma_xy_distribution = []

            for y_pt in y_distribution:
                # 크랙 유무 확인
                if crack_d is not None and y_pt < (-h/2 + crack_d):
                    # 크랙 내부이므로 응력 0
                    sigma_xx_distribution.append(0)
                    sigma_xy_distribution.append(0)
                    continue

                # 단면 2차 모멘트 I 및 중립축 계산
                I_eff = I_effective(b, h, crack_d)
                y_neutral_axis = calculate_neutral_axis(h, crack_d)
                y_corr = y_pt - y_neutral_axis
                Q_y = Q(y_pt, b, h, crack_d)

                # 모멘트 및 전단력 계산
                M_at_selected_x = float(M.subs(x, selected_x).evalf())
                V_at_selected_x = float(V.subs(x, selected_x).evalf())

                # 응력 계산
                sigma_xx_pt = -M_at_selected_x * y_corr / I_eff
                sigma_xy_pt = -V_at_selected_x * Q_y / (I_eff * b)

                sigma_xx_distribution.append(float(sigma_xx_pt))
                sigma_xy_distribution.append(float(sigma_xy_pt))

            # 그래프 그리기
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

    # 다시 계산 여부 확인
    repeat = input("\n다른 위치에서 계산을 진행하시겠습니까? (y/n): ").strip().lower()
    if repeat != 'y':
        print("계산을 종료합니다.")
        break
