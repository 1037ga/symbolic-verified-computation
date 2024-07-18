from sympy import *

def picard(ode,init,deg):
    P = Function('P')
    P = ode 
    Q = Function('Q')
    Q = init 
    n = deg

    for i in range(n):
        if degree(P,x) != 0:
            Q = init + integrate(P.subs(x,Q),t)
            print(Q)
        else:
            Q = init + integrate(P,t)
        # while((degree(Q,t)>i+1)):
        #     Q = Q-LT(Q)
        Q = poly(Q,t)
        Q = expand(extract_polynomial(Q,i+1))

    return Q

def  picard_second_derivative(ode,init1,init2,deg):
    P = Function('P')
    P = ode 
    Q1 = Function('Q1')
    Q1 = init1 
    Q2 = Function('Q2')
    Q2 = init2
    n = deg 

    for i in range(n):
        Q1, Q2 = init1 + integrate(Q2,t), init2 + integrate(ode.subs(x,Q1),t)
        print(Q1)
        print(Q2)
    return Q1

def extract_polynomial(poly, degree):
    terms = [term for term in poly.terms() if term[0][0] <= degree]
    return sum(coeff * t**exp for (exp,), coeff in terms)

if __name__ == "__main__":
    a,x,t,e = symbols('a x t e')
    ode = -x
    init1 = 1
    init2 = 0
    n = 4
    # ans = picard(ode,a,n)
    ans = picard_second_derivative(ode,init1,init2,n)
    print('ans =',ans)