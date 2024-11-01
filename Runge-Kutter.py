from cmath import cos, sin, sqrt, sinh, exp, cosh
from sympy import sin, cos, Matrix, diff, solve, N, symbols
import numpy as np
from sympy.abc import x, y, z, t
import matplotlib.pyplot as plt


def func():
    return Matrix([cos(x) * t, sin(y) + x])


def fx(x, y, z, t):
    global eps
    return -2 * x - 2 * y + sinh(z + y) ** 2 + sin(x) ** 2


def fy(x, y, z, t):
    global eps
    return 2 * x - y + cosh(x + y) * sinh(z + y) ** 2


def fz(x, y, z, t):
    global eps
    return 3 * z + exp((6 * x + y) ** 2) - 1 + x ** 2 * cos(z)


def find_dots(x_0_1, y_0_1, z_0_1):
    t_0 = 0
    h = 0.0001
    points_1 = []
    points_2 = []
    for i in range(200000):
        x1 = fx(x_0_1, y_0_1,z_0_1, t_0)
        y1 = fy(x_0_1, y_0_1,z_0_1, t_0)
        z1 = fy(x_0_1, y_0_1,z_0_1, t_0)
        x2 = fx(x_0_1 + h * x1 / 3, y_0_1 + h * y1 / 3,z_0_1 + h * z1 / 3, t_0 + h / 3)
        y2 = fy(x_0_1 + h * x1 / 3, y_0_1 + h * y1 / 3,z_0_1 + h * z1 / 3, t_0 + h / 3)
        z2 = fy(x_0_1 + h * x1 / 3, y_0_1 + h * y1 / 3,z_0_1 + h * z1 / 3, t_0 + h / 3)
        x3 = fx(x_0_1 - h * x1 / 3 + h * x2, y_0_1 - h * y1 / 3 + h * y2,z_0_1 - h * z1 / 3 + h * z2, t_0 + 2 / 3 * h)
        y3 = fy(x_0_1 - h * x1 / 3 + h * x2, y_0_1 - h * y1 / 3 + h * y2,z_0_1 - h * z1 / 3 + h * z2, t_0 + 2 / 3 * h)
        z3 = fy(x_0_1 - h * x1 / 3 + h * x2, y_0_1 - h * y1 / 3 + h * y2,z_0_1 - h * z1 / 3 + h * z2, t_0 + 2 / 3 * h)
        x4 = fx(x_0_1 + h * x1 - h * x2 + h * x3, y_0_1 + h + y1 - h * y2 + h * y3,z_0_1 + h + z1 - h * z2 + h * z3, t_0 + h)
        y4 = fy(x_0_1 + h * x1 - h * x2 + h * x3, y_0_1 + h + y1 - h * y2 + h * y3,z_0_1 + h + z1 - h * z2 + h * z3, t_0 + h)
        z4 = fy(x_0_1 + h * x1 - h * x2 + h * x3, y_0_1 + h + y1 - h * y2 + h * y3,z_0_1 + h + z1 - h * z2 + h * z3, t_0 + h)
        x_0_1 += h / 8 * (x1 + 3 * x2 + 3 * x3 + x4)
        y_0_1 += h / 8 * (y1 + 3 * y2 + 3 * y3 + y4)
        z_0_1 += h / 8 * (z1 + 3 * z2 + 3 * z3 + z4)
        # print(x_0_1, x_0_2)
        t_0 += h
        points_1.append(x_0_1)
        points_2.append(y_0_1)
        if i%1000==0:
            print(i, t_0, x_0_1  )
    res = []
    res.append(points_1)
    res.append(points_2)
    return res


def print_dots(res):
    plt.plot(res[0], res[1])
    plt.show()


eps = 0.1
x_0_1, y_0_1, z_0_1 = map(float, input('Введите начальное значение x,y,z\n').split())
print_dots(find_dots(x_0_1, y_0_1,z_0_1))
