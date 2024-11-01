# from cmath import e

import sympy
from sympy import Matrix, sinh, sin, cosh, cos, exp, diff, integrate, expand, pretty
from sympy.abc import x, y, z, e, a, b, c, f, d, g, h, s, q, u, r
from sympy.core.backend import sympify

A_plus = Matrix([[3]])
A_minus = Matrix([[-2, - 2], [2, - 1]])
all_variables = (x, y, z)


def method_gauss(matrix, b):
    n = len(b)
    for i in range(n):

        if matrix[i][i] == 0:
            fl = 0
            for j in range(i + 1, n):
                if matrix[j][i] != 0:
                    b[i], b[j] = b[j], b[i]
                    matrix[i], matrix[j] = matrix[j], matrix[i]
                    fl = 1
                    break
            if not fl:
                return None

        bb = matrix[i][i]
        for j in range(i, n):
            matrix[i][j] /= bb
        b[i] /= bb

        for j in range(i + 1, n):
            if matrix[j][i] != 0:
                kk = matrix[j][i]
                for h in range(i, n):
                    matrix[j][h] -= kk * matrix[i][h]
                b[j] -= kk * b[i]
    for i in range(n - 1, -1, -1):
        for j in range(i - 1, -1, -1):
            kk = matrix[j][i]
            matrix[j][i] = 0
            b[j] -= kk * b[i]
    return b


def Telor(f):
    f = expand(f)
    buf_f = f.series(x, 0, 4).removeO().series(y, 0, 4).removeO().series(z, 0, 4).removeO()
    buf_f = expand(buf_f)
    for i1 in all_variables:
        for i2 in all_variables:
            for i3 in all_variables:
                for i4 in all_variables:
                    buff_f = diff(buf_f, i1)
                    buff_f = diff(buff_f, i2)
                    buff_f = diff(buff_f, i3)
                    buff_f = diff(buff_f, i4)
                    buff_f = integrate(buff_f, i1)
                    buff_f = integrate(buff_f, i2)
                    buff_f = integrate(buff_f, i3)
                    buff_f = integrate(buff_f, i4)
                    buff_f = expand(buff_f)
                    buf_f = buf_f - buff_f
                    sympify(buf_f)
    return buf_f


def ff_u(x1, y1, z1):
    # print(x1, y1, z1)
    qq = f_u()
    # print(qq[0], qq[1])
    return Matrix([Telor(qq[0].subs(x, x1).subs(y, y1).subs(z, z1)), Telor(qq[1].subs(x, x1).subs(y, y1).subs(z, z1))])


def f_u():
    return Matrix(
        [[Telor(sinh(z + y) ** 2 + sin(x) ** 2)], [Telor(cosh(x + y) * sinh(z + y) ** 2)]])


def gg_u(x1, y1, z1):
    return Matrix([Telor(g_u()[0].subs(x, x1).subs(y, y1).subs(z, z1))])


def g_u():
    return Matrix([[Telor(exp((6 * x + y) ** 2) - 1 + x ** 2 * cos(z))]])


psi = Matrix([[a * x ** 2 + b * x * y + c * y ** 2 + d * x ** 3 + f * x ** 2 * y + g * x * y ** 2 + h * y ** 3]])
psi_sh = Matrix([[diff(psi[0], x), diff(psi[0], y)]])


# def u_u():
#     return A_minus + ff_u(x, y, z)


def v_v():
    return A_plus + gg_u(x, y, z)


# print(expand((x+y)**3))
# print(ff_u(x,y,psi[0]))
# print(A_plus @ psi + gg_u(x,y,psi[0])-psi_sh@(A_minus @ Matrix([[x], [y]]) ))
# print(expand(    Telor((A_plus @ psi + gg_u(x, y, psi[0]) - psi_sh @ (A_minus @ Matrix([[x], [y]]) + ff_u(x, y, psi[0])))[0])))
expresion1 = expand(
    Telor((A_plus @ psi + gg_u(x, y, psi[0]) - psi_sh @ (A_minus @ Matrix([[x], [y]]) + ff_u(x, y, psi[0])))[0]))
all_combin = ((x, x, x, 6), (x, x, y, 2), (x, y, y, 2), (y, y, y, 6), (x, x, 2), (x, y, 1), (y, y, 2))
all_variables_abc = (a, b, c, d, f, g, h)
expr_buf = expresion1.copy()
matrix_gaus_psi = []
b_psi = []
for i in all_combin:
    rr = expr_buf
    def req(pos):
        global rr, i
        if pos != len(i) - 1:
            rr = diff(rr, i[pos])
            req(pos + 1)
        else:
            return rr


    req(0)
    rr_buf = rr.copy()


    def req1(pos):
        global rr_buf, i
        if pos != len(i) - 1:
            rr_buf = integrate(rr_buf, i[pos])
            req1(pos + 1)
        else:
            return rr_buf


    req1(0)
    rr_buf = expand(rr_buf)
    expr_buf = expr_buf - rr_buf
    rr = rr / i[-1]
    rr = expand(rr)
    cur = []
    for j in all_variables_abc:
        cur_buf = float(diff(rr, j))
        cur.append(cur_buf)
        rr = rr - integrate(cur_buf, j)
    matrix_gaus_psi.append(cur)
    b_psi.append(float(rr))
res = method_gauss(matrix_gaus_psi, b_psi)
print('psi(x,y)=', end=' ')
for i in range(len(all_combin)):
    bb = res[i]
    if bb >= 0:
        print('+', end='')
    for j in range(len(all_combin[i]) - 1):
        bb *= all_combin[i][j]
    print(bb, end=' ')
print()

phi1 = s * z ** 2 + q * z ** 3
phi2 = r * z ** 2 + u * z ** 3
phi1_sh = 2 * s * z + 3 * q * z ** 2
phi2_sh = 2 * r * z + 3 * u * z ** 2

# print(A_minus @ Matrix([phi1, phi2]) + ff_u(phi1, phi2, z)-Matrix([phi1_sh,phi2_sh])@(A_plus*z+gg_u(phi1,phi2,z)))
expresion2 = A_minus @ Matrix([phi1, phi2]) + ff_u(phi1, phi2, z)-Matrix([phi1_sh,phi2_sh])@(A_plus*z+gg_u(phi1,phi2,z))
expresion2 = Matrix([Telor(expresion2[0]),Telor(expresion2[1])])
all_combin_phi = ((z,z,z,6),(z,z,2))
all_variables_abc1 = (s,q,r,u)
expr_buf = expresion2.copy()
matrix_gaus_phi = []
b_phi = []
# print(expresion2)
for i in all_combin_phi:
    for h in range(len(expr_buf)):
        rr = expr_buf[h]
        def req(pos):
            global rr, i
            if pos != len(i) - 1:
                rr = diff(rr, i[pos])
                req(pos + 1)
            else:
                return rr
        req(0)
        rr_buf = rr.copy()
        def req1(pos):
            global rr_buf, i
            if pos != len(i) - 1:
                rr_buf = integrate(rr_buf, i[pos])
                req1(pos + 1)
            else:
                return rr_buf

        req1(0)
        rr_buf = expand(rr_buf)
        expr_buf[h] = expr_buf[h] - rr_buf

        rr = rr / i[-1]
        rr = expand(rr)
        cur = []
        for j in all_variables_abc1:
            cur_buf = float(diff(rr, j))
            cur.append(cur_buf)
            rr = rr - integrate(cur_buf, j)

        matrix_gaus_phi.append(cur)
        b_phi.append(float(rr))
# for i in matrix_gaus_phi:
#     print(i)
# print(b_phi)

res1 = method_gauss(matrix_gaus_phi, b_phi)
# print(res1)
print('phi1(z)=', end=' ')
for i in range(len(all_combin_phi)):
    bb = res1[i]
    if bb >= 0:
        print('+', end='')
    for j in range(len(all_combin_phi[i]) - 1):
        bb *= all_combin_phi[i][j]
    print(bb, end=' ')
print()
print('phi2(z)=', end=' ')
for i in range(len(all_combin_phi)):
    bb = res1[i+len(all_combin_phi)]
    if bb >= 0:
        print('+', end='')
    for j in range(len(all_combin_phi[i]) - 1):
        bb *= all_combin_phi[i][j]
    print(bb, end=' ')
print()
# print(A_minus @ Matrix([[x],[y]])+ff_u(x,y,psi[0]))
