from math import pi, log


# Constants
DENS = 998
MU = 0.001
A_BOX = 0.32 * 0.26
D = 0.00794
A_TUBE = pi * D**2 / 4
EPSILON = 0.0000015
G = 9.81
L = (0.2, 0.6)
TOLERANCE = 0.00001

def reynolds(v):
    return DENS * v * D / MU


def is_laminar(r):
    return r <= 2300


def turbulent_f(r):
    bracket_1 = EPSILON / (3.7 * D)
    bracket_2 = 5.74 / r**0.9
    return 0.25 / (log(bracket_1 + bracket_2)) ** 2


def laminar_f(r):
    return 64 / r


def f(r):
    return laminar_f(r) if is_laminar(r) else turbulent_f(r)


def alpha(r):
    return 2 if is_laminar(r) else 1


def get_V2(z1, a, f):
    denom_3 = a * A_TUBE**2 / A_BOX**2
    denom_2 = f * L[0] / D
    return (2 * G * z1 / (a + denom_2 - denom_3 + 0.5)) ** 0.5


def get_V1(v2):
    return v2 * A_TUBE / A_BOX

def main():
    h = L[0]
    h_i = h
    t = 0
    # Variables
    v2 = 0.6755
    r = reynolds(v2)
    new_v2 = 0

    while h > 0:
        i = 0
        while abs(new_v2 - v2) <= TOLERANCE and i < 1000:
            v2 = new_v2
            r = reynolds(v2)
            new_v2 = get_V2(h+ 0.02 + (L / 150), alpha(r), f(r))
            i += 1

        v1 = get_V1(v2)
        h = h_i - t * v1
        t += 0.1

    print("Final time", t)


if __name__ == "__main__":
    main()
