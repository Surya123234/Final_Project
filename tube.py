from math import pi, log


# Constants
DENS = 998
MU = 0.001
A_BOX = 0.32 * 0.26
D = 0.00794
A_TUBE = pi * D**2 / 4
EPSILON = 0.0000015
G = 9.81
L = (0.2, 0.3, 0.4, 0.6)
TOLERANCE = 0.00001


def reynolds(v):
    return DENS * v * D / MU


def is_laminar(r):
    return r <= 2300


def turbulent_f(r):
    bracket_1 = EPSILON / (3.7 * D)
    bracket_2 = 5.74 / r**0.9
    return 0.25 / (log(bracket_1 + bracket_2, 10)) ** 2


def laminar_f(r):
    return 64 / r


def f(r):
    return laminar_f(r) if is_laminar(r) else turbulent_f(r)


def alpha(r):
    return 2 if is_laminar(r) else 1


def get_V2(z1, a, f, l):
    z1 += 0.02 + l / 150
    denom_1 = a * A_TUBE**2 / A_BOX**2
    denom_2 = f * l / D
    denom = a + denom_1 + denom_2
    return (2 * G * z1 / denom) ** 0.5


def get_V1(v2):
    return v2 * A_TUBE / A_BOX


def main():
    for l in L:
        h = 0.08
        t = 0
        v2 = 0.6755

        while h > 0:
            while True:
                r = reynolds(v2)
                new_v2 = get_V2(h, alpha(r), f(r), l)
                if abs(new_v2 - v2) <= TOLERANCE:
                    break
                v2 = new_v2

            v1 = get_V1(new_v2)
            t += 0.1
            h -= 0.1 * v1

        print(f"Final time for length {l} metres: {int(t)} seconds")


if __name__ == "__main__":
    main()
