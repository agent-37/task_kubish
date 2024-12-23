from cmath import sqrt
import matplotlib.pyplot as plt

from sympy import Matrix, sin, integrate, simplify, expand, solve, diff, pretty, div, collect
from sympy.abc import x, y, l, a

# from cmath import sin,cos

all_root, all_func = [],[]

gamma = Matrix([sin(y), sin(2*y), sin(3*y), y**2, y**4])
v2 = Matrix([sin(y), 2 * sin(2*y), sin(3*y), y**2, y**4])
for _ in range(5):
    mm = Matrix([[0.0, 0.0, 0.0, 0.0, 0.0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]])

    # rr = 0
    # for i in range(5):
    #     rr += v1[i]
    #
    # rr = expand(rr)
    # print(rr)
    # gl_int_sum = sqrt((integrate(expand(rr * rr), (y, 0, 1))).evalf()).real


    for i in range(5):
        for j in range(5):
            if i == j:
                mm[i, j] = (1 - l * (integrate(expand(gamma[j] * v2[i]), (y, 0, 1)))).evalf()
            else:
                mm[i, j] = - l * (integrate(expand(gamma[j] * v2[i]), (y, 0, 1))).evalf()
    print(mm)


    def det(mm):
        bb = 0
        for i1 in range(5):
            for i2 in range(5):
                for i3 in range(5):
                    for i4 in range(5):
                        for i5 in range(5):
                            qq = {i1, i2, i3, i4, i5}
                            if len(qq) == 5:
                                v = [i1, i2, i3, i4, i5]
                                count = 0
                                for i in range(5):
                                    for j in range(5):
                                        if i < j and v[i] > v[j]:
                                            count += 1
                                # print((-1)**count,i1,i2,i3,i4,i5)
                                bb += (-1) ** count * mm[0, i1] * mm[1, i2] * mm[2, i3] * mm[3, i4] * mm[4, i5]
        return bb


    bb = det(mm)
    # print(mm)
    print(expand(bb))
    print(solve(bb))
    res_l = solve(bb)
    for i in range(len(res_l)):
        for j in range(len(res_l)):
            if abs(res_l[i])<abs(res_l[j]):
                res_l[i],res_l[j]=res_l[j],res_l[i]
    root = res_l[0]

    new_mm = mm.copy()
    for i in range(len(new_mm)):
        new_mm[i] = new_mm[i].subs(l, root)
    n = 5
    bb = []
    for i1 in range(5):
        bbb = []
        for i2 in range(5):
            bbb.append(new_mm[i1, i2])
        bb.append(bbb)
    new_mm = bb.copy()
    # for i in new_mm:
    #     print(i)
    # print()

    for i in range(n):
        # if new_mm[i][i] == 0:
        #     for j in range(i + 1, n):
        #         if new_mm[j][i] != 0:
        #             new_mm[i], new_mm[j] = new_mm[j], new_mm[i]
        #             fl = 1
        #             break

        bb = new_mm[i][i]
        for j in range(i, n):
            new_mm[i][j] /= bb

        for j in range(i + 1, n):
            if new_mm[j][i] != 0:
                kk = new_mm[j][i]
                for h in range(i, n):
                    new_mm[j][h] -= kk * new_mm[i][h]
    # for i in new_mm:
    #     print(i)
    # print()
    b = [0, 0, 0, 0, 1]

    for i in range(n - 1, -1, -1):
        for j in range(i - 1, -1, -1):
            b[j] -= b[i] * new_mm[j][i]
            new_mm[j][i] = 0
    # print(b)

    u = 0
    for i in range(5):
        u += b[i] * gamma[i]
    # print(integrate(expand(u * u), (y, 0, 1)).evalf(),"!!!!!!!!!!!!!")
    u /= sqrt(integrate(expand(u * u), (y, 0, 1)).evalf()).real

    ww = 0
    for i in range(5):
        ww += gamma[i]* v2[i].subs(y, x)

    print()
    print(root)
    print(u)
    all_root.append(root)
    all_func.append(u)
    print(ww - u * u.subs(y, x) / root)
    print(expand(ww - u * u.subs(y, x) / root))
    www = expand(ww - u * u.subs(y, x) / root)
    k_buf = www.copy()
    for i in gamma:
        k_buf = collect(k_buf,i)
    for i in range(len(gamma)):
        v2[i]= k_buf.coeff(gamma[i]).subs(x,y)
    print(k_buf)
    print(v2)
    print('-------------------------------------------------------------------------------')
for i in range(len(all_root)):
    print(f'lambda_{i+1} = {all_root[i]}; func_{i+1} = {all_func[i].subs(y,x)}')
f = sin(3*x) + a*x**2
for i in range(len(all_root)):
    # print())
    bb  = integrate(expand(f * all_func[i].subs(y,x)), (x, 0, 1)).evalf()
    print( bb, '///', (bb- a* bb.coeff(a)) / bb.coeff(a) )
    # print(f'f_{i+1}={ sqrt(integrate(expand(f * all_func[i].subs(y,x)), (x, 0, 1)).evalf()).real}')
dh=0.001
points=[]
pos =0
res = [[],[],[],[],[]]
while pos<1:
    for i in range(len(all_func)):
        # print(all_func[i].subs(y,pos).evalf(),pos)
        res[i].append(all_func[i].subs(y,pos).evalf())
    points.append(pos)
    pos+= dh
for i in range(len(res)):
    plt.plot(points,res[i],label=f'u_{i+1}')
plt.legend(loc='upper left')
plt.show()
