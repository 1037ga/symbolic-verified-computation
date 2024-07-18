import sympy  
import sympy.sets
from sympy import *
from sympy.solvers.solveset import solvify
from tabulate import tabulate
import logging
import time

class Interval:
    def __init__(self, lower, upper):
        self.lower = lower
        self.upper = upper
        if compare_polynomials(lower <= upper):
            return
        else:
            raise ValueError(str(self)+': '+str(lower)+' > '+str(upper))

    def __repr__(self):
        return f'[{self.lower}, {self.upper}]'

    def __add__(self, other):
        if isinstance(other, Interval):
            return Interval(self.lower + other.lower, self.upper + other.upper)
        else:
            other = Interval(other,other)
            return Interval(self.lower + other.lower, self.upper + other.upper)
        
    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, Interval):
            return Interval(self.lower - other.upper, self.upper - other.lower)
        else:
            other = Interval(other,other)
            return Interval(self.lower - other.upper, self.upper - other.lower)
    
    def __rsub__(self, other):
        return self.__sub__(other)
    
    def __mul__(self, other):
        if isinstance(other, Interval):
            products = [self.lower * other.lower, self.lower * other.upper,
                        self.upper * other.lower, self.upper * other.upper]
            if (compare_polynomials(self.upper < 0)):
                if (compare_polynomials(other.upper < 0)):
                    return Interval(products[3],products[0])
                elif (compare_polynomials(other.lower > 0)):
                    return Interval(products[1],products[2])
                elif (compare_polynomials(other.lower <= 0)) and (compare_polynomials(other.upper >= 0)):
                    return Interval(products[1],products[0])
                else:
                    raise ValueError('__mul__error1'+'\n'+str(self)+'\n'+str(other))
            elif (compare_polynomials(self.lower > 0)):
                if (compare_polynomials(other.upper < 0)):
                    return Interval(products[2],products[1])
                elif (compare_polynomials(other.lower > 0)):
                    return Interval(products[0],products[3])
                elif (compare_polynomials(other.lower <= 0)) and (compare_polynomials(other.upper >= 0)):
                    return Interval(products[2],products[3])
                else:
                    raise ValueError('__mul__error2'+'\n'+str(self)+'\n'+str(other))
            elif (compare_polynomials(self.lower <= 0)) and (compare_polynomials(self.upper >= 0)):
                if (compare_polynomials(other.upper < 0)):
                    return Interval(products[2],products[0])
                elif (compare_polynomials(other.lower > 0)):
                    return Interval(products[1],products[3])
                elif (compare_polynomials(other.lower <= 0)) and (compare_polynomials(other.upper >= 0)):
                    if compare_polynomials(self.lower * other.upper <= self.upper * other.lower):
                        min = products[1]
                    else:
                        min = products[2]
                    if compare_polynomials(self.lower * other.lower <= self.upper * other.upper):
                        max = products[3]
                    else:
                        max = products[0]
                    return Interval(min,max)
                else:
                    raise ValueError('__mul__error3'+'\n'+str(self)+'\n'+str(other))
            else:
                raise ValueError('__mul__error4'+'\n'+str(self)+'\n'+str(other))
        else:
            if compare_polynomials(other >= 0):
                return Interval(self.lower * other, self.upper * other)
            elif compare_polynomials(other < 0):
                return Interval(self.upper * other, self.lower * other)
            else:
                raise ValueError('__mul__error5'+'\n'+str(self)+'\n'+str(other))

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
            if isinstance(other,Interval):
                # if other.lower <= 0 <= other.upper:
                #     raise ValueError("Interval division by an interval containing zero is not defined.")
                # quotients = [self.lower / other.lower, self.lower / other.upper,
                #             self.upper / other.lower, self.upper / other.upper]
                # return Interval(min(quotients), max(quotients))
                raise ValueError('未実装')
            elif isinstance(other,(int,float)):
                other = Interval(other,other)
                if other == Interval(0,0):
                    raise ValueError("Interval division by an interval containing zero is not defined.") 
                if compare_polynomials(self.upper/other.lower >= self.lower/other.lower):
                    return Interval(self.lower/other.lower,self.upper/other.lower)
                elif compare_polynomials(self.upper/other.lower < self.lower/other.lower):
                    return Interval(self.upper/other.lower,self.lower/other.lower)
                else:
                    raise ValueError('__truediv__error1'+'\n'+str(self)+'\n'+str(other))
            else:
                raise ValueError('__truediv__error2'+'\n'+str(self)+'\n'+str(other))

    def __neg__(self):
        return -1*Interval(self.lower, self.upper)
    
    def __pow__(self, power):
        if compare_polynomials(self.lower < 0) and not isinstance(power, int):
            raise ValueError("Negative base with non-integer exponent is not supported in interval arithmetic.")
        if isinstance(power, int):
            if power % 2 == 0:
                if compare_polynomials(self.lower < 0) and compare_polynomials(self.upper > 0):
                    if compare_polynomials(abs(self.lower) < abs(self.upper)):
                        return Interval(0, self.upper**power)
                    elif compare_polynomials(abs(self.lower) >= abs(self.upper)):
                        return Interval(0, self.lower**power)
                    else:
                        raise ValueError('__pow__error1'+'\n'+str(self)+'\n'+str(power))
                elif compare_polynomials(self.lower >= 0) or compare_polynomials(self.upper <= 0):
                    if compare_polynomials(abs(self.lower) < abs(self.upper)):
                        return Interval(self.lower**power, self.upper**power)
                    elif compare_polynomials(abs(self.lower) >= abs(self.upper)):
                        return Interval(self.upper**power, self.lower**power)
                    else:
                        raise ValueError('__pow__error2'+'\n'+str(self)+'\n'+str(power))
            else:
                return Interval(self.lower**power, self.upper**power)
        else:
            return Interval(self.lower**power, self.upper**power)

    def __eq__(self, other):
        if isinstance(other,Interval):
            if (self.lower == other.lower) and (self.upper == other.upper):
                return True
            else:
                return False
        else:
            raise ValueError('Error_ep')

    def norm(self):
        """
        intervalクラスの絶対値の最大値
        """
        if compare_polynomials(self.lower > 0):
            return self.upper
        elif compare_polynomials(self.upper < 0):
            return -self.lower
        elif (compare_polynomials(self.lower <= 0)) and (compare_polynomials(self.upper >= 0)):
            if compare_polynomials(-self.lower > self.upper):
                return  -self.lower
            elif compare_polynomials(-self.lower <= self.upper):
                return self.upper
            else:
                raise ValueError('norm_error1'+'\n'+str(self))
        else:
            raise ValueError('norm_error2'+'\n'+str(self))

    
    def expand(self):
        """
        intervalクラスの上限,下限の展開
        """
        return Interval(expand(self.lower),expand(self.upper))
    
    def subset(self,other):
        '''
        selfがotherの部分集合であるかをbool値で返す関数
        '''
        if isinstance(other,Interval):
            tmp1 = self.lower - other.lower
            tmp2 = self.upper - other.upper
            # 要改善 tmp3=True,tmp4=False,の時にtmp3がtrueの時のeの条件を保持してしまっている. 
            tmp3 = compare_polynomials(tmp1 >= 0)
            tmp4 = compare_polynomials(tmp2 <= 0)
            return tmp3 and tmp4
        else:
            raise ValueError("error_subset")

    def width(self):
        return self.upper - self.lower

def abs(expr):
    if compare_polynomials(expr >= 0):
        return expr
    elif compare_polynomials(expr < 0):
        return -expr
    else:
        raise ValueError('abs_error'+'\n'+str(expr))

def sym():
    """
    symbolic変数の更新
    """
    global e,t,a,b,c,x,syms
    e,t,a,b,c,x = symbols('e t a b c x')
    e = Symbol('e', real = true)
    syms = [e,t,a,b,c,x]

def horner(coefficients, t):
    """
    horner法を用いた多項式の評価関数
    coefficients:[t^0の係数,t^1の係数,t^2の係数,...]
    t:多項式の変数
    """
    coefficients.reverse()
    result = coefficients[0]
    for coeff in coefficients[1:]:
        result = result * t + coeff
    return result

def picard(ode,init,deg):
    P = sympy.Function('P')
    P = ode #常微分方程式の右辺
    Q = sympy.Function('Q')
    Q = init #関数x(t)の初期値
    n = deg #最高次数

    for i in range(n):
        if degree(P,x) != 0:
            Q = init + integrate(P.subs(x,Q),t)
        else:
            Q = init + integrate(P,t)
        while((degree(Q,t)>i+1)):
            Q = Q-LT(Q)

    return Q
    
        
def Guaranteed_accuracy(v_c,coefficients):
    itv()
    ans = [0 for _ in range(len(coefficients))]
    print(coefficients)
    for i in range(len(coefficients)):
        ans[len(coefficients)-i-1] = eval(str(coefficients[i]))
    v = Interval.expand(degree_n+(Rational(1,n+1)*horner(ans[n:],t))*t)
    ans = Interval.subset(v,v_c)
    tmp = [v.lower,v.upper,v_c.lower,v_c.upper,ans]
    output_table.append(tmp)
    if ans == False:   
        return False
    else:
        return v
           
def polynomial_coefficients(expr,t):
    expr = poly(expr,t)
    coefficients = expr.all_coeffs()
    return coefficients

def output(output_table,option=0):
    if option == 0:
        header = ['num','V.lower','V.upper','V_c.lower','V_c.upper','V ⊆ V_c']
        print(tabulate(output_table, headers=header),'\n')
    else:
        for o in output_table:
            print(o[0],'回目')
            print('V = '+'['+str(o[1])+', '+str(o[2])+']')
            print('V_c = '+'['+str(o[3])+', '+str(o[4])+']')
            print('V ⊆ V_c =',o[5],'\n')

def sympy_interval_to_inequality(itv,var):
    """
    sympy.Interval型の変数 itv を var の不等式(str型)にして返す関数
    """
    if isinstance(itv,(sympy.Interval)):
        if itv.is_open:
            return str(itv.inf)+' < '+ var +' < '+str((itv.sup))
        elif itv.left_open:
            return str(itv.inf)+' < '+ var +' ≦ '+str((itv.sup))
        elif itv.right_open:
            return str(itv.inf)+' ≦ '+ var +' < '+str((itv.sup))
        else:
            return str(itv.inf)+' ≦ '+ var +' ≦ '+str((itv.sup))
    elif isinstance(itv,(Union)):
        return str(itv)
    elif isinstance(itv,FiniteSet):
        return str(var)+' = '+str((itv))
    else:
        raise ValueError('type('+str(itv)+') is '+ str(type(itv)) +'. It is not sympy.Interval')

def has_CRootOf(poly):
    """
    多項式の係数にCRootOfが含まれているかどうかを判別する関数
    """
    # 多項式の係数を取得
    coeffs = poly.coeffs()
    
    # 係数にCRootOfが含まれているかチェック
    for coeff in coeffs:
        if any(isinstance(sympy.Abs(arg), CRootOf) for arg in coeff.args):
            return True
    return False

def compare_polynomials(expr):
    """
    expr を満たすならば True, 満たさないならば False を返す関数
    expr: 不等式
    """
    global e_range

    if isinstance(expr,bool):
        return expr
    if isinstance(expr,sympy.logic.boolalg.BooleanTrue):
        return True
    elif isinstance(expr,sympy.logic.boolalg.BooleanFalse):
        return False
    logging.debug('e_range = '+str(e_range))
    logging.debug('expr = '+str(expr))
    comparison = str(expr.rel_op)
    expr = expand(expr.lhs - expr.rhs)
    start = time.time()
    tmp1 = solve_poly_inequality(poly(expr,e),comparison)
    end = time.time()
    if tmp1 == []:
        logging.debug(str(False)+'\n')
        return False
    ans =  Intersection(tmp1[0],e_range)
    for i in range(len(tmp1)):
        if i > 0:
            tmp = Intersection(tmp1[i],e_range)
            if ans == EmptySet:
                ans = tmp
            else:
                if tmp != EmptySet:
                    ans = Union(ans,tmp)
    logging.debug('time: '+str(end-start))
    logging.debug('ans = '+str(ans))
    if ans == EmptySet:
        logging.debug(str(False)+'\n')
        return False
    elif approach_interval.is_subset(ans):
        if ans == e_range:
            logging.debug(str(True)+'\n')
            return True
        else:
            e_range = ans
            logging.debug(str(True))
            logging.debug(sympy_interval_to_inequality(e_range,'e') +'\n       の範囲で成立'+'\n')
            return True
    else:
        logging.debug(str(False)+'\n')
        return False

def itv():
    """
    intervalクラスの更新
    具体値はRationalで与える
    """
    global a,b,t
    # a = Interval(Rational(1,10),Rational(1,10))
    # t = Interval(0,50)
    a = Interval(Rational(1,1)-e,Rational(1,1)+e)
    t = Interval(0,Rational(1,10))

def main():
    global t, degree_n, n, output_table, e_range, approach_interval

    sym()
    e_range = sympy.Interval(-oo,oo)
    approach_value = 0
    approach_interval = sympy.Interval.open(approach_value,approach_value+Rational(1,10**10))
    # approach_interval = sympy.Interval.open(approach_value-Rational(1,10**10),approach_value)
    # approach_interval = EmptySet
    output_table = []
    # ode = k*(1-x/L)*x
    # ode = a*(1 - x/(500))*x
    ode = -x**2
    init = a
    n = 10
    print('dx/dt =',ode)

    # 近似解の生成 
    expr1 = picard(ode,init,n)
    print('x(t) =',expand(expr1))
    itv()
    print('x(0) =',init,'\n')
    sym()

    # 解候補区間の生成
    expr2 = expand(ode.subs(x,expr1),t)
    coeffs = polynomial_coefficients(expr2,t)
    itv()
    ans = [0 for _ in range(len(coeffs))]
    for i in range(len(coeffs)):
        ans[len(coeffs)-i-1] = eval(str(coeffs[i]))
    degree_n = Rational(1,n)*ans[n-1]
    r = Interval.norm(degree_n+(Rational(1,n+1)*horner(ans[n:],t))*t-degree_n)
    v_c = Interval.expand(degree_n + 2*r*Interval(-1,1))
    tmp = v_c

    # 解の精度保証  
    sym()
    v_c = symbols('v_c')
    expr3 = expand(ode.subs(x,(v_c*t**n+expr1-LT(expr1))),t)
    coeffs = polynomial_coefficients(expr3,t)
    tmp = Guaranteed_accuracy(tmp,coeffs)
    if isinstance(tmp,Interval):
        v_c = tmp
        for i in range(0):
            v_c = Guaranteed_accuracy(v_c,coeffs)
            if isinstance(v_c,bool):
                break
    output_table = [[i + 1] + row for i, row in enumerate(output_table)]
    output(output_table,1)
    
    # X(t.upper)の計算
    if isinstance(tmp,Interval):
        t = Interval(t.upper,t.upper)
        X = (Interval.expand(eval(str(expr1-LT(expr1)))+v_c*t**n))
        print('x(' + str(t.upper) + ') =', X)
        # print('width(x) =', Interval.width(X))
        print('N(width(x)) =', N(limit(Interval.width(X),e,0)),'\n')
        print(sympy_interval_to_inequality(e_range,'e')+'の範囲で成立'+'\n')
        print(e_range.inf, ' = ', (e_range.inf).evalf(20))
        print(e_range.sup, ' = ', (e_range.sup).evalf(20))


if __name__ == "__main__":
    # logging.basicConfig(level=logging.DEBUG, format='%(levelname)s: %(message)s')
    main()
    sym()
    e_range=sympy.Interval(0,oo)
    approach_interval = 0