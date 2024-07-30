from sympy import *

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

def polynomial_continued_fraction(coeffs,n):
    """
    多項式の係数を連分数展開で有理数近似
    coeffs: 多項式の係数
    n: 打ち切り回数 
    """
    for i,coeff in enumerate(coeffs):
        terms = continued_fraction(numer(coeff),denom(coeff))
        coeff = fraction_from_cf(terms[:n])
        coeffs[i] = coeff
    return coeffs

if __name__ == '__main__':
    a = 90909
    b = 100000
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