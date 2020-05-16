# Codigo de metodo variacional a un oscilador armonico unidimensional
import numpy as np
from matplotlib import pyplot as plt

print('---- Metodo Variacional aplicado en Oscilador Unidimensional----')
#Definir constantes
hbar = 1
k = 1
mu = 1
omega = np.sqrt(k/mu)

#vector x
x = np.linspace(-20,20,500)

#funcion de onda de prueba
def FOPrueba(alpha,xt):
    return np.exp(-alpha*xt*xt)

#funcion que toma las derivadas
def dfdx(ft,xt):
    dx=x[1]-x[0]
    ftp = np.zeros_like(ft)
    for i in range(0,len(ft)):
        if(i<(len(ft)-1)):
            dif = ft[i+1]-ft[i]
            ftp[i] = dif/dx
        else:
            dif = ft[i]-ft[i-1]
            ftp[i] = dif/dx

    return ftp

#crear una funcion para ver como el operador T
# actua sobre Phi (que son ft en la funcion de python)
def TPhi(ft,xt):
    ftp = dfdx(ft,xt)
    ftpp = dfdx(ftp,xt)
    return -0.5*ftpp

#crear funcion para ver como V actua sobre phi
#(es ft en esta funcion de python)

def VPhi(ft,xt):
    ftv = 0.5*k*(xt**2)*ft
    return ftv

#ahora vamos a usar esto para encontrar las funcionales T y V
def TFuncional(ft,xt):
    tphi = TPhi(ft,xt)
    dx = x[1]-x[0]
    num = 0
    denom = 0
    #hacer la integral
    for i in range(0,len(ft)):
        num = num + ft[i]*tphi[i]*dx
        denom = denom + ft[i]*ft[i]*dx

    return num/denom

def VFuncional(ft,xt):
    vphi = VPhi(ft,xt)
    dx = x[1] - x[0]
    num = 0
    denom = 0
    # hacer la integral
    for i in range(0, len(ft)):
        num = num + ft[i] * vphi[i] * dx
        denom = denom + ft[i] * ft[i] * dx

    return num / denom

#crear array diferentes valores de alpha
#para hacer un grafico
alpha = np.zeros(20)
Evals = np.zeros(20)
Tvals = np.zeros(20)
Vvals = np.zeros(20)
Eg = np.zeros(20)


#Loop para ver las curvas
for j in range(0,20):
    alpha[j] = 0.05*(j+3)
    Phi = FOPrueba(alpha[j],x)
    Tvals[j] = TFuncional(Phi,x)
    Vvals[j] = VFuncional(Phi,x)
    Evals[j] = Tvals[j] + Vvals[j]
    Eg[j] = np.sqrt(k/mu)*0.5

plt.xlabel(r'$\alpha$')
plt.ylabel(r'$Energy$')
plt.title('Variational Oscillator´s Energy', fontsize=12)
plt.plot(alpha,Evals,'red', label = 'Total Energy')
plt.plot(alpha,Tvals,'blue',label = 'Kinetic Energy')
plt.plot(alpha,Vvals,'purple',label = 'Potential Energy')
plt.plot(alpha,Eg,'black',label = 'Ground State Energy')
plt.legend(loc='lower right', fontsize= 7)

plt.show()

#Minimizacion E con respecto a alpha.
j = 0
alpha_min = 0
Emin = 0
while (Evals[j] - Eg[1]) > 0.0001:
    alpha_min = alpha[j]
    Emin = Evals[j]
    j+=1

print('Resultados:')
print('Energia del estado base:  ', Eg[1])
print('Energia total minima:  ', Emin)
print('Energia total - Energia estado base:  ', Emin - Eg[1])
print('Error de la minimizacion (%):  ', ((Emin - Eg[1])/Eg[1])*100)
print('El valor alpha de minimizacion:  ', alpha_min)

#Grafico de la funcion de onda
plt.xlabel('x')
plt.ylabel(r'$\psi$')
plt.title('Oscillator´s Wave Function', fontsize=12)
Phi = FOPrueba(alpha_min,x)
Phip =dfdx(Phi,x)
Phipp =dfdx(Phip,x)
plt.plot(x,Phi,'red',label = r'$\psi$')
plt.plot(x,Phip,'blue', label =r'$\frac{d \psi}{dx}$', linewidth=0.8)
plt.plot(x,Phipp,'black', label =r'$\frac{d^2 \psi}{dx^2}$', linewidth=0.8)
plt.legend(loc='lower right', fontsize= 10)
plt.text(13.0, 0.70, r'$\psi = e^{-\alpha x^2}$ ', fontsize=12)
plt.text(13.0, 0.50, r'$\alpha$ = %.1e' %(alpha_min), fontsize=10)
plt.show()