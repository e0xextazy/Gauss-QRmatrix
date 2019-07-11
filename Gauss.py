import numpy as np
import math
import pylab
from matplotlib import mlab
import random

def func(x):
    return math.sin(x)

def ScalMul(a, b): # Скалярное произведение
    n = len(a)
    y = 0
    for i in range(n):
        y += a[i]*b[i]
    return y

def Coefs(x, k):
    coefs = [0] * (k) # Заполнение матрицы с коэф перед С_i
    for i in range(k):
        coefs[i] = [0] * (k)
    
    for i in range(k):
        for j in range(k):
            coefs[i][j] = ScalMul([pow(p,i) for p in x], [pow(p,j) for p in x])
    return coefs

def RightVector(x, y, k):
    z = [0] * (k)
    for i in range (k):
        z[i] = ScalMul([pow(p,i) for p in x], y)
    return z

def QR(A, z): # QR Разложение матрицы
    n = len(z)
    c = [0] * (n) 
    Result = [0] * (n) # Вектор содержащий C_i
    for i in range(n):
        c[i] = [0] * (n)
    s = [0] * (n) 
    for i in range(n):
        s[i] = [0] * (n)
    # Прямой ход
    for k in range(0, n-1, 1): # Исключение переменных
        for l in range(k+1, n, 1):
            c[k][l] = A[k][k] / (sqrt( A[k][k]*A[k][k] + A[l][k]*A[l][k] ))
            s[k][l] = A[l][k] / (sqrt( A[k][k]*A[k][k] + A[l][k]*A[l][k] ))
            # Умножение матрицы A[][] на G[k][l]
            akk = A[k][k]
            alk = A[l][k]
            akl = A[k][l]
            all = A[l][l]
            A[k][k] =  akk*c[k][l] + alk*s[k][l]
            A[k][l] =  akl*c[k][l] + all*s[k][l]
            A[l][k] = -akk*s[k][l] + alk*c[k][l]
            A[l][l] = -akl*s[k][l] + all*c[k][l]
            # Вектор свободных членов умножается на G[k][l]
            zk = z[k]
            zl = z[l]
            z[k] =  zk*c[k][l] + zl*s[k][l]
            z[l] = -zk*s[k][l] + zl*c[k][l]
    # Теперь матрица A[][] вверх-диагональная 
    # Обратный ход
    Result[n-1] = z[n-1] / A[n-1][n-1]
    for l in range(n-1, 0, -1):
        h = z[l-1]
        for k in range(l+1, n+1, 1):
            h = h - Result[k-1] * A[l-1][k-1]
        Result[l-1] = h / A[l-1][l-1]
    return Result

def Approx(ABC, k, x): # Строим функцию
    Res = 0
    for i in range(k):
        Res += ABC[i] * pow(x, i)
    return Res
    
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    x = np.array([-2,-1,0,1,2])
    #y = np.array([5.3, 6.3, 4.8, 3.8, 3.3]) # Задача y с помощью набора точек
    size = len(x)-1
    y = [func(l)+random.uniform(0, 0.5)*math.pow(-1,l) for l in x] # Задача y с помощью функции
    n = len(x) # Число точек
    k = 4 # Размерность пр-ва базисных фун-й
    C = Coefs(x, k) # Матрица с коэф перед C_i
    z = RightVector(x, y, k)
    
    xlist = mlab.frange( x[0]-2, x[size]+2, 0.01 )
    ylist = [ Approx(QR(C, z), k, v) for v in xlist ] # Значения приближенной функции
    yFlist = [ func(p) for p in xlist ] # Значения данной функции
    
    print(Approx(QR(C, z), k, 1))
    plt.scatter(x, y, color="black") # Узлы функции
    pylab.plot(xlist, ylist, color="blue") # Построение приближенной функции
    pylab.plot(xlist, yFlist,"--", color="red") # Построение функции
    pylab.show()