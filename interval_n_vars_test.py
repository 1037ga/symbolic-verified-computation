import sympy  
import sympy.sets
from sympy import *
from sympy.solvers.solveset import solvify
import logging
import time
import sys
from continued_fraction import *

FRACTION_LOWER_NUM = 31
FRACTION_UPPER_NUM = FRACTION_LOWER_NUM+1

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
    
# def pprint(expr):
    expr = str(expr)
    expr = expr.replace('**','^')
    expr = expr.replace('*','')
    expr = expr.replace('e','ε')
    print(expr)

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
    
        
def Guaranteed_accuracy(v_c,coefficients,degree_n):
    t,a,b,c = itv_list
    length = len(coefficients)
    v_c0,v_c1,v_c2,v_c3,v_c4,v_c5,v_c6,v_c7,v_c8,v_c9 = v_c
    Ans = [0 for _ in range(length)]
    for i in range(length):
        t,a,b,c = itv_list
        ans = [0 for _ in range(len(coefficients[i]))]
        for j in range(len(coefficients[i])):
            ans[j] = eval(str(coefficients[i][j]))
        ans.reverse()
        Ans[i] = ans
    is_subset = [0 for _ in range(length)]
    v = [0 for _ in range(10)]
    for i in range(length):
        if Ans[i][n:] == []:
            v[i] = degree_n[i]+0*t
        else:
            v[i] = Interval.expand(degree_n[i]+(Rational(1,n+1)*horner(Ans[i][n:],t))*t)
            coeffs_lower = poly_to_coeffs(v[i].lower,e)
            coeffs_lower.reverse()
            coeffs_lower = polynomial_continued_fraction(coeffs_lower,FRACTION_LOWER_NUM)
            coeffs_upper = poly_to_coeffs(v[i].upper,e)
            coeffs_upper.reverse()
            coeffs_upper = polynomial_continued_fraction(coeffs_upper,FRACTION_UPPER_NUM)
            v[i] = Interval(coeffs_to_poly(coeffs_lower,e),coeffs_to_poly(coeffs_upper,e))
        print('v_c'+str(i),'= ',v_c[i])
        print(' v'+str(i),' = ',v[i])
        is_subset[i] = Interval.subset(v[i],v_c[i])
        print('v'+str(i),'⊆','v_c'+str(i),'=',is_subset[i],'\n')
    if False in is_subset:
        return False
    else:
        return v
    
def poly_to_coeffs(expr,var):
    expr = poly(expr,var)
    coefficients = expr.all_coeffs()
    return coefficients

def sympy_interval_to_inequality(itv,var,decimal = false):
    """
    sympy.Interval型の変数 itv を var の不等式(str型)にして返す関数
    """
    if isinstance(itv,(sympy.Interval)):
        if decimal == False:
            if itv.is_open:
                return str(itv.inf)+' < '+ var +' < '+str((itv.sup))
            elif itv.left_open:
                return str(itv.inf)+' < '+ var +' ≦ '+str((itv.sup))
            elif itv.right_open:
                return str(itv.inf)+' ≦ '+ var +' < '+str((itv.sup))
            else:
                return str(itv.inf)+' ≦ '+ var +' ≦ '+str((itv.sup))
        else:
            inf = itv.inf.evalf(20)
            sup = itv.sup.evalf(20)
            if itv.is_open:
                return str(inf)+' < '+ var +' < '+str((sup))
            elif itv.left_open:
                return str(inf)+' < '+ var +' ≦ '+str((sup))
            elif itv.right_open:
                return str(inf)+' ≦ '+ var +' < '+str((sup))
            else:
                return str(inf)+' ≦ '+ var +' ≦ '+str((sup))
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
    tmp1 = solve_poly_inequality(poly(expr,e),comparison)
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
    logging.debug('ans = '+str(ans))
    if ans == EmptySet:
        logging.debug(str(False)+'\n')
        return False
    elif approach_interval.is_subset(ans):
        if ans == e_range:
            logging.debug(str(True)+'\n')
            return True
        else:
            if ans.sup.has(CRootOf):
                terms = continued_fraction(floor(ans.sup.evalf(30)*10**30),10**30)
                newans = sympy.Interval.open(ans.inf,fraction_from_cf(terms[:31]))
                e_range = newans
            else:
                e_range = ans
            logging.debug(str(True))
            logging.debug(sympy_interval_to_inequality(e_range,'e') +'\n の範囲で成立'+'\n')
            return True
    else:
        logging.debug(str(False)+'\n')
        return False

def coeffs_to_poly(coeffs,var):
    poly = sum(coeff * var **i for i, coeff in enumerate(coeffs))
    return poly

def psa_varified(ode,init):
    global e

    # 近似解の生成 
    init_length = len(init)
    expr1 = picard(ode,init,n)
    for i in range(init_length):
        print('d'+str(x[i])+'/dt = '+str(ode[i]))
        print(str(x[i])+'(t) =',expand(expr1[x[i]]))
        t,a,b,c = itv_list
        print(str(x[i])+'('+str(t.upper*(div_num))+') =',str(init[i]),'=',eval(str(init[i])),'\n')
        sym()
        t,a,b,c = symbols('t a b c')

    # 解候補区間の生成
    expr2 = [0 for _ in range(init_length)]
    coeffs = [0 for _ in range(init_length)]
    for i in range(init_length):
        expr2[i] = expand(ode[i].subs(expr1),t)
        coeffs[i] = poly_to_coeffs(expr2[i],t)
    Ans = [0 for _ in range(init_length)]
    t,a,b,c = itv_list
    for i in range(init_length):
        ans = [0 for _ in range(len(coeffs[i]))]
        for j in range(len(coeffs[i])):
                if isinstance(coeffs[i][j],Rational):
                    ans[j] = nsimplify(eval(str(coeffs[i][j])))
                else:
                    ans[j] = eval(str(coeffs[i][j]))
        ans.reverse()
        Ans[i] = ans
    degree_n = [0 for _ in range(init_length)]
    v_c = [0 for _ in range(10)]
    for i in range(init_length):
        if len(Ans[i]) < n:
            degree_n[i] = 0
        else:
            degree_n[i] = Rational(1,n)*Ans[i][n-1]
        if Ans[i][n:] == []:
            r = Rational(1,10**(n-2))
            v_c[i] = Interval.expand(degree_n[i] + 2*r*Interval(-1,1)) 
        else:
            r = Interval.norm((degree_n[i]+(Rational(1,n+1)*horner(Ans[i][n:],t))*t)-degree_n[i])
            # r = Rational(1,10**(n-2))
            if r.subs(e,0) == 0:
                r = Rational(1,10**(n-1))
            v_c[i] = Interval.expand(degree_n[i] + 2*r*Interval(-1,1)) 
            coeffs_lower = poly_to_coeffs(v_c[i].lower,e)
            coeffs_lower.reverse()
            coeffs_lower = polynomial_continued_fraction(coeffs_lower,FRACTION_LOWER_NUM)
            coeffs_upper = poly_to_coeffs(v_c[i].upper,e)
            coeffs_upper.reverse()
            coeffs_upper = polynomial_continued_fraction(coeffs_upper,FRACTION_UPPER_NUM)
            v_c[i] = Interval(coeffs_to_poly(coeffs_lower,e),coeffs_to_poly(coeffs_upper,e))
    tmp = v_c 

    # 解の精度保証  
    sym()
    t,a,b,c = symbols('t a b c')
    v_c = symbols('v_c0:10')
    expr3 = [0 for _ in range(init_length)]
    for i in range(init_length):
        if degree(LT(expr1[x[i]]),t) == n:
            expr1[x[i]] = v_c[i]*t**n+expr1[x[i]]-LT(expr1[x[i]])
        else:
            expr1[x[i]] = v_c[i]*t**n+expr1[x[i]]
    for i in range(init_length):
        expr3[i] = expand(ode[i].subs(expr1),t)
        coeffs[i] = poly_to_coeffs(expr3[i],t)
    print('----- 1 -----')
    tmp = Guaranteed_accuracy(tmp,coeffs,degree_n)
    if isinstance(tmp,bool):
        return 0
    else:
        for i in range(1,3):
            print('----- '+str(i+1)+' -----')
            tmp = Guaranteed_accuracy(tmp,coeffs,degree_n)
    print()

    value = [0 for _ in range(init_length)]
    for i in range(init_length):
        sym()
        t,a,b,c = symbols('t a b c')
        coeffs[i] = poly_to_coeffs(expr1[x[i]],t)
        t,a,b,c = itv_list
        for j in range(1,n+1):
            coeffs[i][j] = eval(str(coeffs[i][j]))
        coeffs[i][0] = tmp[i]
        coeffs[i].reverse()
        out = dict()
        for j,c in enumerate(coeffs[i]):
            out[j] = c
        print(str(x[i])+'(t) =',out)
        t = Interval(t.upper,t.upper)
        value[i] = Interval.expand(horner(coeffs[i],t))
        print(str(x[i])+'('+str(t.upper*(div_num+1))+') =',value[i])
        value_width = Interval.width(value[i])
        print('width('+str(x[i])+'('+str(t.upper)+'))','=',value_width,'=',value_width.evalf(20))
        print()

    print(sympy_interval_to_inequality(e_range,'e')+' の範囲で成立'+'\n')
    print(sympy_interval_to_inequality(e_range,'e',True)+' の範囲で成立'+'\n')
    
    for i in range(init_length):
        value_lower = poly_to_coeffs(value[i].lower,e)
        value_upper = poly_to_coeffs(value[i].upper,e)
        value_lower.reverse()
        value_upper.reverse()
        if e_range.sup.has(CRootOf):
            e = Interval(e_range.inf,floor(e_range.sup.evalf(30)*10**30)/10**30)
        else:
            e = Interval(e_range.inf,e_range.sup)
        if  len(value_lower) > 2:
            value_lower_lower = horner(value_lower[1:],e).lower
        elif len(value_lower) == 2:
            value_lower_lower = value_lower[1]
        else:
            value_lower_lower = 0
        if len(value_upper) > 2:
            value_upper_upper = horner(value_upper[1:],e).upper
        elif len(value_upper) == 2:
            value_upper_upper = value_upper[1]
        else:
            value_upper_upper = 0
        e = symbols('e',real=True)
        coeffs_lower = [value_lower[0],value_lower_lower]
        coeffs_upper = [value_upper[0],value_upper_upper]
        coeffs_lower = polynomial_continued_fraction(coeffs_lower,FRACTION_LOWER_NUM)
        coeffs_upper = polynomial_continued_fraction(coeffs_upper,FRACTION_UPPER_NUM)
        value[i] = Interval(coeffs_lower[0]+coeffs_lower[1]*e,coeffs_upper[0]+coeffs_upper[1]*e)
    return value

def sym():
    """
    symbolic変数の更新
    """
    global a,b,c,x,y,e,t
    t,a,b,c,y = symbols('t a b c y')
    e = Symbol('e', real = true)
    x = symbols('x0:11')

def main():
    global e_range,approach_interval,n,itv_list,div_num
    # logging.basicConfig(level=logging.DEBUG, format='%(levelname)s: %(message)s')
    start = time.time()
    sym()
    e_range = sympy.Interval(-oo,oo)
    approach_interval = sympy.Interval.open(0,0+Rational(1,10**10))

    t_start,t_end = 0,Rational(1,1)
    ode = [-x[0]**2]
    init = [a]
    n = int(sys.argv[1])
    div = int(sys.argv[2])
    itv_list = [0 for _ in range(4)] #t,a,b,c
    itv_list[0] = Interval(0,Rational(t_end-t_start,div))
    itv_list[1] = Interval(Rational(1,1)-e,Rational(1,1)+e)
    itv_list[2] = Interval(1,1)
    itv_list[3] = Interval(1,1)
    ans = []
    for div_num in range(div):
        # print(div_num)
        value = psa_varified(ode,init)
        if isinstance(value,int):
            return 0
        for i in range(len(value)):
            itv_list[i+1] = value[i]
        ans.append(value)
    f = open('FRACTION_LOWER_NUM=31_2.txt', 'a')
    f.write('n = '+str(n)+', div = '+str(div)+', ode = '+str(ode)+', t_end = '+str(t_end)+'\n')
    f.write(str(x[i])+'('+str(t_end)+') = '+str(ans[div-1][0])+'\n')
    tmp = len(str(x[i])+'('+str(t_end)+') ')
    tmp = ' '*tmp
    f.write(tmp + '= ['+str(ans[div-1][0].lower.evalf(30))+', '+str(ans[div-1][0].upper.evalf(30))+']'+'\n')
    tmp = Interval.width(ans[div-1][0])
    f.write('width = '+str(tmp)+'\n')
    f.write('      = '+str(tmp.evalf(30))+'\n')
    f.write('e_range_sup = '+str(e_range.sup)+'\n')
    f.write('            = '+str(e_range.sup.evalf(30))+'\n')
    end = time.time()
    time_diff = end - start
    f.write('time:' + str(time_diff)+'\n'+'\n')

main()
sym()
