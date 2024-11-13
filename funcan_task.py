from PIL.ImageOps import expand
from mpmath import zeros
from sympy import Matrix, sin, integrate, simplify, expand, solve, diff
from sympy.abc import x, y, l
# from cmath import sin,cos
v1 = Matrix([2*sin(y), sin(3*y),y**2,y**3,y**4])
mm = Matrix([[0.0,0.0,0.0,0.0,0.0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0]])

bonus = 1

for i in range(5):
    for j in range(5):
        if i==j:
            mm[i,j]=(1-l*(integrate(v1[i]*v1[j],(y,0,1))))*bonus
        else:
            mm[i, j] = - l * (integrate(v1[i]*v1[j],(y,0,1)))*bonus
print(mm)
bb=0
for i1 in range(5):
    for i2 in range(5):
        for i3 in range(5):
            for i4 in range(5):
                for i5 in range(5):
                    qq={i1,i2,i3,i4,i5}
                    if len(qq)==5:
                        v=[i1,i2,i3,i4,i5]
                        count=0
                        for i in range(5):
                            for j in range(5):
                                if i<j and v[i]>v[j]:
                                    count+=1
                        # print((-1)**count,i1,i2,i3,i4,i5)
                        bb+= (-1)**count*mm[0,i1]*mm[1,i2]*mm[2,i3]*mm[3,i4]*mm[4,i5]
# print(mm)
print(expand(bb))
print(diff(expand(bb),l,l,l,l,l)/120)
bb= expand(bb).evalf()
# bb = bb- l**5*diff(expand(bb),l,l,l,l,l)/120
print(bb)
print(solve(bb))
res_l = solve(bb)

for i in range(5):
    for j in range(5):
        mm[i,j]=mm[i,j].evalf()/bonus
for root in res_l:
    new_mm = mm.copy()
    for i in range(len(new_mm)):
        new_mm[i]= new_mm[i].subs(l,root)
    print(new_mm.det())
    n = 5
    bb=[]
    for i1 in range(5):
        bbb=[]
        for i2 in range(5):
            bbb.append(new_mm[i1,i2])
        bb.append(bbb)
    new_mm=bb
    # for i in new_mm:
    #     print(i)
    # print()
    # for i in range(n):
    #     if new_mm[i][i] == 0:
    #         for j in range(i + 1, n):
    #             if new_mm[j][i] != 0:
    #                 new_mm[i], new_mm[j] = new_mm[j], new_mm[i]
    #                 fl = 1
    #                 break
    #
    #
    #     bb = new_mm[i][i]
    #     for j in range(i, n):
    #         new_mm[i][j] /= bb
    #
    #
    #     for j in range(i + 1, n):
    #         if new_mm[j][i] != 0:
    #             kk = new_mm[j][i]
    #             for h in range(i, n):
    #                 new_mm[j][h] -= kk * new_mm[i][h]
    # for i in new_mm:
    #     print(i)
    # print()
    # for i in range(n - 1, -1, -1):
    #     for j in range(i - 1, -1, -1):
    #         kk = new_mm[j][i]
    #         new_mm[j][i] = 0


