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

def continued_fraction(a, b):
    """
    連分数展開のアルゴリズム
    """
    terms = []
    while b != 0:
        q = floor(a / b)
        terms.append(q)
        a, b = b, a - q * b
    return terms

def fraction_from_cf(terms):
    if not terms:
        return 0
    frac = terms[-1]
    for term in reversed(terms[:-1]):
        frac = term + 1 / frac
    return frac

a = 31415926535897932384626433832795028841971
b = 10000000000000000000000000000000000000000
terms = continued_fraction(a, b)
print(str(len(str(a)))+'桁','a =',a)
print(str(len(str(b)))+'桁','b =',b)
print('連分数表示',' 長さ',len(terms))
print(terms)
print()
for i in range(len(terms)):
    ans = fraction_from_cf(terms[:i+1])
    print('打切',i+1)
    print('近似',ans)
    print('誤差',(ans-Rational(a,b)).evalf(30))
    if (ans-Rational(a,b)) >= 0:
        tmp = '+'
    else:
        tmp = '-'
    print('正負',tmp)
    print()

# 0.499996184243842 - 4.48922186945969*e
# 4.47074439125908*e + 0.50000315299614
a = 12794935963538478830399/26895352859468526256128
b = 1138012408557978581962263057189318521782674817/2170080016305965495537271176977755960772657152
print(a,b)