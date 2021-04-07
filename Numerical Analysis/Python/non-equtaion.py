# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 16:31:52 2021

@author: yzz
"""

def set_function(x):
    """
    定义函数
    """
    return x*x*x -x-1

def set_a_and_b(x):
    """
    求解区间[a,b]
    以x=0.6为基础，以精度为0.1向两边寻找区间[a,b]
    """
    a = x -0.1
    b= x + 0.1
    while set_function(a)*set_function(b)>0:
        a = a -0.1
        b = b +0.1
    return a,b

def dic_solve_non_liner(x):
    """
    二分法
    1、获取区间[a,b]
    2、设置误差项为0.0001---求解
    3、输出(x*,k)-----精确解 x*,迭代步数 k
    """
    a = set_a_and_b(x)[0]
    b = set_a_and_b(x)[1]
    k = 0
    
    while set_function((a+b)/2)>0.00001 or set_function((a+b)/2)<-0.00001:
        
        if set_function((a+b)/2)*set_function(b)<0:
            a = (a+b)/2
        else:
            b = (a + b) / 2
        k = k + 1   
    return (a+b)/2,k



def interation_solve_non_liner(x):
    """
    不动点迭代
    1、构造迭代公式
    2、输出精确解x*,迭代步数
    """
    import math
    k = 0
    
    while (x-math.pow(x+1, 1/3))>0.00001 or (x-math.pow(x+1, 1/3))<-0.00001:
        
        x = math.pow(x+1, 1/3)
        k = k + 1   
    return x,k


def Aiyken_interation_solve_non_liner(x):
    """
    Aiyken不动点迭代
    1、利用其迭代公式....
    2、输出精确解x*,迭代步数
    """
    import math
    k = 0
    y = math.pow(x+1, 1/3)
    z = math.pow(y+1, 1/3)
    while (math.pow(x+1, 1/3)-x+((y-x)*(y-x)/(z-2*y+x)))>0.00001 or(math.pow(x+1, 1/3)-x+((y-x)*(y-x)/(z-2*y+x)))<-0.00001:
        x=x-((y-x)*(y-x)/(z-2*y+x))
        y = math.pow(x+1, 1/3)
        z = math.pow(y+1, 1/3)
        k = k + 1   
    return x,k

def newton_interation_solve_non_liner(x):
    """
    牛顿下山法
    1、lamada=1，构造迭代公式....
    2、输出精确解x*,迭代步数
    """
    import math
    k = 0
    while ((math.pow(x, 3)-x-1)/(3*x*x-1))>0.00001 or ((math.pow(x, 3)-x-1)/(3*x*x-1))<-0.00001:
        x = x -(math.pow(x, 3)-x-1)/(3*x*x-1)
        k = k+1
    return x,k




def newton_l_interation_solve_non_liner(x):
    """
    牛顿下山法
    1、lamada=1，构造迭代公式....
    2、输出精确解x*,迭代步数
    """
    import math
    lamada = 1
    k = 0
    while (lamada*(math.pow(x, 3)-x-1)/(3*x*x-1))>0.00001 or (lamada*(math.pow(x, 3)-x-1)/(3*x*x-1))<-0.00001:
        y = x -lamada*(math.pow(x, 3)-x-1)/(3*x*x-1)
     
        if abs(set_function(y))-abs(set_function(x))>0:
            lamada = lamada/2
        else:
             break
        x = y
        k = k+1
    return x,k




def newton_x_interation_solve_non_liner(x,x1):
    """
    牛顿---弦截法
    1、设置新的初值x1，构造迭代公式....
    2、输出精确解x*,迭代步数
    """
    import math
    k=0
    while ((x1-x)/(1-set_function(x)/set_function(x1)))>0.00001 or ((x1-x)/(1-set_function(x)/set_function(x1)))<-0.00001:
        #注意为啥直接弄了一个y=x出来
        y = x
        x = x1
        x1 = x1-((x1-y)/(1-set_function(y)/set_function(x1)))
        
        k = k+1
    return x1,k

