# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 11:24:44 2020

@author: Titi
"""

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import sympy as sym

np.set_printoptions(suppress=True)

#Fizikai paraméterek
m=0.1
l=1.
k=5.
F=11.

#Kezdeti értékek
Fi_10=0.0
Fi_20=-0.1
dFi_10=0.
dFi_20=0.

#Szimulációs paraméterek
dt=0.02 #ha túl nagy akkor a program nem tud futni (talán az egyenletrendszer megoldó talál több lehetséges megoldást)
steps=1000 #lépésszám 

#A rendszer pillanatbeli pozícióját, sebességét, gyorsulását
# egy 2*3-as mátrixba fogom tárolni
#Ezeket a 2*3-as mátrixokat rendezem egy 3D mátrixba lépésenként

ddFi_1 = sym.symbols('ddFi_1')
ddFi_2 = sym.symbols('ddFi_2')

#Fi hipermátrixba beírja az utolsó pillanatbeli gyorsulás vektort
def ddFI(Fi_T):
    eq1=m*l**2*(4/3*ddFi_1+1/2*ddFi_2*np.cos(Fi_T[-1,0,0]-Fi_T[-1,0,1])+1/2*Fi_T[-1,0,1]**2*np.sin(Fi_T[-1,0,0]-Fi_T[-1,0,1]))+2*k*Fi_T[-1,0,0]-k*Fi_T[-1,0,1]-F*l*np.sin(Fi_T[-1,0,0]-Fi_T[-1,0,1])
    eq2=m*l**2*(1/3*ddFi_2+1/2*ddFi_1*np.cos(Fi_T[-1,0,0]-Fi_T[-1,0,1])-1/2*Fi_T[-1,0,0]**2*np.sin(Fi_T[-1,0,0]-Fi_T[-1,0,1]))+k*Fi_T[-1,0,1]-k*Fi_T[-1,0,0]
    ddFi=sym.solve((eq1,eq2), (ddFi_1, ddFi_2))
    Fi_T[-1,2,:]=[ddFi[ddFi_1],ddFi[ddFi_2]]

#Fi hipermátrixból visszaadja a XY koordinátákat 
def conv(Fi):
    x1=l*np.cos(Fi[0,0,0])
    y1=l*np.sin(Fi[0,0,0])
    x2=x1+l*np.cos(Fi[0,0,1])
    y2=y1+l*np.sin(Fi[0,0,1])
    XY=np.array([[x1, y1, x2, y2]])
    for i in range(steps-1):
        x1=l*np.cos(Fi[i+1,0,0])
        y1=l*np.sin(Fi[i+1,0,0])
        x2=x1+l*np.cos(Fi[i+1,0,1])
        y2=y1+l*np.sin(Fi[i+1,0,1])
        XY_next=np.array([[x1, y1, x2, y2]])
        XY=np.append(XY, XY_next, axis=0)
    return(XY)

#Fi null
Fi = np.array([[[Fi_10, Fi_20], [dFi_10, dFi_20],[0., 0.]]])
ddFI(Fi)


#Fi 1
Fi_next=np.array([[[10., 10.], [0., 0.],[0., 0.]]])
Fi_next[0,1,:]=Fi[-1,1,:]+dt*Fi[-1,2,:]
Fi_next[0,0,:]=Fi[-1,0,:]+dt*(Fi[-1,1,:]+Fi_next[0,1,:])/2
Fi=np.append(Fi, Fi_next, axis=0)
ddFI(Fi)


#Fi 2-steps -> (Runge Kutta (dd-1, -2 -> d | d 0, -1 ->Fi | Fi 0, d 0 -> dd ))
for i in range(steps-2):
    Fi_next[0,1,:]=Fi[-1,1,:]+dt*(3*Fi[-1,2,:]-Fi[-2,2,:])/2
    Fi_next[0,0,:]=Fi[-1,0,:]+dt*(Fi[-1,1,:]+Fi_next[0,1,:])/2
    Fi=np.append(Fi, Fi_next, axis=0)
    ddFI(Fi)

#Fi vektor kirajzoltatása (ez az igazán hasznos)
plt.figure(num=None, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
plt.title('m= %1.1f'  %m+ ' kg, l= %1.1f' %l+' m , k= %1.1f'%k+ 'N/rad, F= %1.1f' %F+ 'N')
plt.plot(range(steps),Fi[:,0,0], label="Fi1")
plt.plot(range(steps),Fi[:,0,1],label="Fi2")
plt.show()


#A mozgás kirajzolása egy GIF-be (Ez látványos) 
#koordináták konvertálása
XY=conv(Fi)
#egyes képek elkészítésa
fig = plt.figure(steps, figsize=(6, 6))
ims = []
for i in range(steps):
    im = plt.plot([0, XY[i,0],XY[i,2]], [0, XY[i,1],XY[i,3]],color='red')
    plt.axis([-2, 2, -2, 2])
    plt.title('m= %1.1f'  %m+ ' kg, l= %1.1f' %l+' m , k= %1.1f'%k+ 'N/rad, F= %1.1f' %F+ 'N')
    #plt.title('m= %1.2f, ' %m + 'm= %1.2f, ' %l )
    ims.append(im)
#képek összefűzése , és mentése 2dof_koveto.gif-ként.    
ani = animation.ArtistAnimation(fig, ims, interval=dt)
ani.save("2dof_koveto.gif",writer='pillow')


