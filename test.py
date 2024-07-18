from sympy import *
from sympy.abc import *

x = symbols('x0:11')
ode = [x[1],-x[0]] #x0'(t) = ode[0], x1'(t) = ode[1], ...
init = [1,0]   #x0(0) = init[0], x1(0) = init[1], ...
n = 4

def picard(ode,init,deg): 
    tmp = [0 for _ in range(len(init))]
    n = deg
    values = dict()
    for i in range(len(init)):
        values[x[i]] = init[i]


    for i in range(n):
        for j in range(len(init)):
            tmp[j] = init[j] + integrate(ode[j].subs(values),t)
        for j in range(len(init)):
            while(degree(tmp[j],t)>i+1):
                tmp[j] = tmp[j]-LT(tmp[j])
            values[x[j]] = tmp[j]

    return values

def f(n):
    print(n)
    if n > 100:
        return n-10
    return f(f(n+11))

print(f(1))