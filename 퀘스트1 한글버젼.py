import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

# plt.ion()  # 대화형 모드 비활성화

# 변수 정의
x, y = sp.symbols('x y')
L = float(input("보의 길이 L을 입력하세요 (m): "))  # 단위: 미터 (m)
b = float(input("보의 단면 폭 b을 입력하세요 (m): "))  # 단위: 미터 (m)
h = float(input("보의 단면 높이 h을 입력하세요 (m): "))  # 단위: 미터 (m)
I = b * h**3 / 12  # 단면 2차 모멘트 (단위: m^4)

# 지점 반력 변수 정의
R_A, R_B = sp.symbols('R_A R_B')  # 단위: 뉴턴 (N)

# 포인트 모멘트 하중 입력
num_point_moments = int(input("포인트 모멘트 하중의 개수를 입력하세요: "))
moment_positions = []  # 위치 단위: 미터 (m)
moment_magnitudes = []  # 크기 단위: 뉴턴미터 (Nm)

for i in range(num_point_moments):
    pos = float(input(f"{i+1}번째 모멘트 하중 위치 x 값을 입력하세요 (m): "))
    mag = float(input(f"{i+1}번째 모멘트 하중 크기 값을 입력하세요 (Nm, 시계 반대 방향 양수): "))
    moment_positions.append(pos)
    moment_magnitudes.append(mag)

# 포인트 수직하중 입력 (아래 방향 하중은 양수로 입력)
num_point_loads = int(input("포인트 수직하중의 개수를 입력하세요: "))
point_positions = []  # 위치 단위: 미터 (m)
point_magnitudes = []  # 크기 단위: 뉴턴 (N)

for i in range(num_point_loads):
    pos = float(input(f"{i+1}번째 포인트 하중 위치 x 값을 입력하세요 (m): "))
    mag = float(input(f"{i+1}번째 포인트 하중 크기 값을 입력하세요 (N, 아래 방향 양수): "))
    point_positions.append(pos)
    point_magnitudes.append(mag)

# 연속 수직하중 입력 (아래 방향 하중은 양수로 입력)
num_continuous_loads = int(input("균등 분포하중의 개수를 입력하세요: "))
continuous_loads = []  # 위치 단위: 미터 (m), 하중 크기 단위: 뉴턴/미터 (N/m)

for i in range(num_continuous_loads):
    start_pos = float(input(f"{i+1}번째 균등 분포하중 시작 위치 x 값을 입력하세요 (m): "))
    end_pos = float(input(f"{i+1}번째 균등 분포하중 끝 위치 x 값을 입력하세요 (m): "))
    magnitude = float(input(f"{i+1}번째 균등 분포하중 크기 값을 입력하세요 (N/m, 아래 방향 양수): "))
    continuous_loads.append((start_pos, end_pos, magnitude))

# 선형 연속하중 입력 (아래 방향 하중은 양수로 입력)
num_linear_loads = int(input("선형 연속하중의 개수를 입력하세요: "))
linear_loads = []  # 위치 단위: 미터 (m), 시작 및 끝 하중 크기 단위: 뉴턴/미터 (N/m)

for i in range(num_linear_loads):
    start_pos = float(input(f"{i+1}번째 선형 하중 시작 위치 x 값을 입력하세요 (m): "))
    end_pos = float(input(f"{i+1}번째 선형 하중 끝 위치 x 값을 입력하세요 (m): "))
    start_mag = float(input(f"{i+1}번째 선형 하중 시작 크기 값을 입력하세요 (N/m): "))
    end_mag = float(input(f"{i+1}번째 선형 하중 끝 크기 값을 입력하세요 (N/m): "))
    linear_loads.append((start_pos, end_pos, start_mag, end_mag))

# 전단력 함수 정의 (히블러의 부호 규약에 맞게 설정)
V = R_A  # 전단력, 반력 R_A로 시작 (단위: 뉴턴, N)
M = sp.Integer(0)  # 굽힘 모멘트를 심볼릭 0으로 초기화 (단위: 뉴턴미터, Nm)

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
boundary_conditions = [M.subs(x, 0) - 0]

# 전체 하중 계산
total_point_load = sum(point_magnitudes)  # 단위: N
total_cont_load = sum(mag * (end_pos - start_pos) for start_pos, end_pos, mag in continuous_loads)  # 단위: N
total_linear_load = sum(0.5 * (start_mag + end_mag) * (end_pos - start_pos) for start_pos, end_pos, start_mag, end_mag in linear_loads)  # 단위: N
total_load = total_point_load + total_cont_load + total_linear_load  # 단위: N

# 힘의 평형 방정식
force_eq = sp.Eq(R_A + R_B, total_load)

# 모멘트 평형 방정식 (전체 보 길이에서 모멘트 합계)
moment_eq = sp.Eq(M.subs(x, L), 0)

# 방정식 풀기
solutions = sp.solve([force_eq] + boundary_conditions + [moment_eq], (R_A, R_B, C1))
R_A_sol = solutions[R_A]  # 단위: N
R_B_sol = solutions[R_B]  # 단위: N
C1_sol = solutions[C1]    # 단위: Nm

# R_A, R_B, C1 값을 V와 M에 대입
V = V.subs(R_A, R_A_sol)  # 단위: N
M = M.subs(R_A, R_A_sol).subs(C1, C1_sol)  # 단위: Nm

# 응력 분포 함수 (단위: Pa, N/m^2)
sigma = -M * y / I

# 결과 출력 (응력 함수 출력 제거)
print("\n왼쪽 지점 반력 R_A (단위: N):", R_A_sol)
print("오른쪽 지점 반력 R_B (단위: N):", R_B_sol)

# 중립축에서 y 위치까지의 1차 모멘트 Q(y) 계산 함수 (단위: m^3)
def Q(y_val, b, h):
    return (b / 2) * ((h**2 / 4) - y_val**2)

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
plt.show()  # 블로킹 모드으로 표시

# 불연속점 목록 생성
discontinuities = sorted(
    set(
        moment_positions +
        point_positions +
        [pos for load in continuous_loads for pos in (load[0], load[1])] +
        [pos for load in linear_loads for pos in (load[0], load[1])]
    )
)

# 불연속점인지 확인하는 함수
def is_discontinuity(x_val, discontinuities, epsilon=1e-6):
    return any(abs(x_val - pos) < epsilon for pos in discontinuities)

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

        # Q(y) 계산 (단위: m^3)
        Q_at_y = Q(y_val, b, h)

        if is_discontinuity(x_val, discontinuities):
            # 불연속점인 경우 좌극한과 우극한 계산
            epsilon = 1e-6
            x_left = x_val - epsilon
            x_right = x_val + epsilon

            # 해당 위치에서의 전단력 계산 (단위: N)
            V_at_x_left = float(V.subs(x, x_left).evalf())
            V_at_x_right = float(V.subs(x, x_right).evalf())

            # 전단 응력 σ_xy 계산 (단위: N/m^2)
            sigma_xy_left = V_at_x_left * Q_at_y / (I * b)
            sigma_xy_right = V_at_x_right * Q_at_y / (I * b)

            # 수직 응력 σ_xx 계산 (단위: N/m^2)
            sigma_xx_left = float(sigma.subs({x: x_left, y: y_val}).evalf())
            sigma_xx_right = float(sigma.subs({x: x_right, y: y_val}).evalf())

            # 응력 텐서 구성
            sigma_xz = 0  # z 방향 전단 응력은 무시 (단위: N/m^2)
            sigma_yy = 0  # 단위: N/m^2
            sigma_zz = 0  # 단위: N/m^2

            # 좌측 응력 텐서 구성 (단위: N/m^2)
            stress_tensor_left = sp.Matrix([
                [sigma_xx_left, sigma_xy_left, sigma_xz],
                [sigma_xy_left, sigma_yy, 0],
                [sigma_xz, 0, sigma_zz]
            ])

            # 우측 응력 텐서 구성 (단위: N/m^2)
            stress_tensor_right = sp.Matrix([
                [sigma_xx_right, sigma_xy_right, sigma_xz],
                [sigma_xy_right, sigma_yy, 0],
                [sigma_xz, 0, sigma_zz]
            ])

            # 응력 텐서의 크기 계산 (절댓값의 합)
            norm_left = abs(sigma_xx_left) + abs(sigma_xy_left)
            norm_right = abs(sigma_xx_right) + abs(sigma_xy_right)

            # 큰 응력 텐서 선택
            if norm_left >= norm_right:
                stress_tensor = stress_tensor_left
                print(f"\n({x_val} m, {y_val} m, {z_val} m) 위치에서의 응력 텐서 (단위: N/m^2):")
            else:
                stress_tensor = stress_tensor_right
                print(f"\n({x_val} m, {y_val} m, {z_val} m) 위치에서의 응력 텐서 (단위: N/m^2):")

            # 응력 텐서 출력
            sp.pprint(stress_tensor)

        else:
            # 불연속점이 아닌 경우
            # 해당 위치에서의 전단력 계산 (단위: N)
            V_at_x = float(V.subs(x, x_val).evalf())

            # 수직 응력 σ_xx 계산 (단위: N/m^2)
            sigma_xx = float(sigma.subs({x: x_val, y: y_val}).evalf())

            # 전단 응력 σ_xy 계산 (단위: N/m^2)
            sigma_xy = V_at_x * Q_at_y / (I * b)

            # 응력 텐서 구성
            sigma_xz = 0  # z 방향 전단 응력은 무시 (단위: N/m^2)
            sigma_yy = 0  # 단위: N/m^2
            sigma_zz = 0  # 단위: N/m^2

            stress_tensor = sp.Matrix([
                [sigma_xx, sigma_xy, sigma_xz],
                [sigma_xy, sigma_yy, 0],
                [sigma_xz, 0, sigma_zz]
            ])

            print(f"\n({x_val} m, {y_val} m, {z_val} m) 위치에서의 응력 텐서 (단위: N/m^2):")
            sp.pprint(stress_tensor)

    except Exception as e:
        print(f"오류 발생: {e}")

    # 다시 계산 여부 확인
    repeat = input("\n다른 위치에서 계산을 진행하시겠습니까? (y/n): ").strip().lower()
    if repeat != 'y':
        print("계산을 종료합니다.")
        break

# 인터랙티브 모드 비활성화 (필요 시)
# plt.ioff()
# plt.show()  # 최종적으로 한 번 더 호출할 필요 없음
