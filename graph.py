from sympy import *
from sympy.abc import *
from spb import *

x, a, b, c = symbols("x, a, b, c")
plot(2*atan(tanh(t/2)),-0.000253131513548180*t**9 - 61*t**7/5040 + t**5/24 - t**3/6 + t,0.000253131513548180*t**9 - 61*t**7/5040 + t**5/24 - t**3/6 + t)
print(2*atan(tanh(1/10/2)),-0.000253131513548180*t**9 - 61*t**7/5040 + t**5/24 - t**3/6 + t,0.000253131513548180*t**9 - 61*t**7/5040 + t**5/24 - t**3/6 + t)