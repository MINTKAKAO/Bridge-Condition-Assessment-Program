import sympy as sp

# Load 클래스 정의
class Load:
    def __init__(self, load_type, parameters):
        self.load_type = load_type  # 1: 포인트 모멘트, 2: 포인트 수직하중, 3: 연속 수직하중
        self.parameters = parameters

# 사용자 입력 처리
def get_user_input():
    length = float(input("교량 길이 (m): "))
    height = float(input("단면 높이 (m): "))
    width = float(input("단면 너비 (m): "))
    loads = []

    num_loads = int(input("하중의 개수: "))

    for i in range(num_loads):
        print(f"\n하중 #{i+1}")
        print("하중의 종류를 선택하세요:")
        print("1: 포인트 모멘트 하중")
        print("2: 포인트 수직하중")
        print("3: 연속 수직하중")
        load_type = int(input("하중 종류 (1/2/3): "))

        if load_type == 1:
            position = float(input("하중의 위치 x (m): "))
            magnitude = float(input("하중의 크기 (Nm): "))
            parameters = {'position': position, 'magnitude': magnitude}
            loads.append(Load(1, parameters))
        elif load_type == 2:
            position = float(input("하중의 위치 x (m): "))
            magnitude = float(input("하중의 크기 (N): "))
            parameters = {'position': position, 'magnitude': magnitude}
            loads.append(Load(2, parameters))
        elif load_type == 3:
            start = float(input("하중의 시작 위치 x1 (m): "))
            end = float(input("하중의 끝 위치 x2 (m): "))
            magnitude = float(input("하중의 크기 (N/m): "))
            parameters = {'start': start, 'end': end, 'magnitude': magnitude}
            loads.append(Load(3, parameters))
        else:
            print("잘못된 입력입니다.")

    return length, height, width, loads

# 단면의 모멘트 관성 계산
def moment_of_inertia(height, width):
    return (width * height**3) / 12

# 단순 지지 보의 모멘트 계산 함수
def beam_moment(length, loads):
    x = sp.Symbol('x')
    M = 0

    for load in loads:
        if load.load_type == 1:  # 포인트 모멘트
            position = load.parameters['position']
            magnitude = load.parameters['magnitude']
            M += magnitude * sp.Piecewise((0, x < position), (1, x >= position))
        elif load.load_type == 2:  # 포인트 수직하중
            position = load.parameters['position']
            magnitude = load.parameters['magnitude']
            M += magnitude * sp.Piecewise(
                (0, x < position),
                (x - position, (x >= position) & (x <= length)),
                (0, True)
            )
        elif load.load_type == 3:  # 연속 수직하중
            start = load.parameters['start']
            end = load.parameters['end']
            magnitude = load.parameters['magnitude']
            M += sp.Piecewise(
                (0, x < start),
                (magnitude * (x - start)**2 / 2, (x >= start) & (x <= end)),
                (magnitude * (end - start) * (x - (end + start)/2), x > end)
            )

    return M

# 응력 계산 함수
def beam_stress(length, height, width, loads):
    x = sp.Symbol('x')
    y = sp.Symbol('y')
    I = moment_of_inertia(height, width)
    M = beam_moment(length, loads)

    sigma_x = -M * y / I
    return sigma_x

# 응력 텐서 계산 함수
def stress_tensor_at_point(sigma_x, x_val, y_val, z_val):
    sigma_xx = sigma_x.subs({'x': x_val, 'y': y_val})
    sigma_yy = 0  # 보의 축 방향 변형이 없다고 가정
    sigma_zz = 0  # 폭 방향 응력은 무시

    # 전단 응력 계산 (필요 시 추가)
    tau_xy = 0
    tau_xz = 0
    tau_yz = 0

    # 응력 텐서 구성
    stress_tensor = sp.Matrix([[sigma_xx, tau_xy, tau_xz],
                               [tau_xy, sigma_yy, tau_yz],
                               [tau_xz, tau_yz, sigma_zz]])
    return stress_tensor

# 메인 함수
def main():
    length, height, width, loads = get_user_input()
    sigma_x = beam_stress(length, height, width, loads)

    while True:
        x_val = float(input("x 위치 (m): "))
        y_val = float(input("y 위치 (m): "))
        z_val = float(input("z 위치 (m): "))

        stress_tensor = stress_tensor_at_point(sigma_x, x_val, y_val, z_val)
        print("응력 텐서:")
        sp.pprint(stress_tensor)

        cont = input("다른 위치에서 응력을 계산하시겠습니까? (y/n): ")
        if cont.lower() != 'y':
            break

if __name__ == "__main__":
    main()
