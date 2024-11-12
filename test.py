import sympy
from sympy import *

def rump():
    a = 77617.0
    b = 33096.0
    a = int(a)
    b = int(b)
    c = 333.75*b**6 + a**2*(11*a**2*b**2-b**6-121*b**4-2)+5.5*b**8+a*1/(2*b)
    print(c)
    a = float(a)
    b = float(b)
    c = 333.75*b**6 + a**2*(11*a**2*b**2-b**6-121*b**4-2)+5.5*b**8+a*1/(2*b)
    print(c)
    a = Rational(a)
    b = Rational(b)
    c = Rational(33375,100)*b**6 + a**2*(11*a**2*b**2-b**6-121*b**4-2)+ Rational(55,10)*b**8+a*1/(2*b)
    print(c.evalf())

def seven_square_sum(n):
    n = Integer(n)
    ans = 1/24*n**2*(n+1)**2*(3*n**4+6*n**3-n**2-4*n+2)
    print(N(ans,100))
    n = Float(n)
    ans = 1/8*n**2*(n+1)**2*(n**4+2*n**3-n**2/3-4*n/3+2/3)
    print(N(ans,100))
    n = Rational(n,1)
    ans = Rational(1,8)*n**2*(n+1)**2*(n**4+2*n**3-Rational(n**2,3)-Rational(4*n,3)+Rational(2,3))
    print(ans)
    ans = n**2*(n+1)**2*(3*n**4+6*n**3-n**2-4*n+2)*1/24
    print(ans)
    
def continued_fraction(a, b):
    """
    連分数展開のアルゴリズム
    """
    terms = []
    i = 0
    while b != 0:
        q = floor(a / b)
        terms.append(q)
        a, b = b, a - q * b
        i += 1
        if i == 101:
            break
    return terms

def fraction_from_cf(terms,n):
    if not terms:
        return 0
    frac = terms[n]
    for term in reversed(terms[:n]):
        frac = term + 1 / frac
    return frac

def main():
    a = 2**(Rational(1,3))
    b = 1
    terms = continued_fraction(a,b)
    print(terms)
    for n in range(len(terms)):
        ans = fraction_from_cf(terms,n)
        diff = (ans - a/b)/(a/b)
        print(ans)
        print(diff.evalf())

if __name__ == '__main__':
    main()