import matplotlib.pyplot as plt
import math
import os

def zlisting(tpr_min,ppr_step,tpr_step,tpr_max,ppr_max,y0,tolerence):
    z_list=[]
    load = int(100/int((tpr_max-tpr_min)/tpr_step))
    i = 1
    print('LOADING...')
    while tpr_min <= tpr_max:
        z_list.append(zmatrix(tpr_min,ppr_step,ppr_max,y0,tolerence))
        tpr_min = tpr_min + tpr_step
        os.system('cls')
        print('LOADING', load*i, "%")
        i = i + 1
    return z_list

def zmatrix(tpr,ppr_step,ppr_max,y0,tolerence):
    z_array = []
    a = 0.0
    t = 1/tpr
    yerror = 5
    while a <= ppr_max:
        y0 = 0.06125 * ( a) * t * (math.exp((-1.2 * ((1 - t)**2)))) + 0.001
        yerror = yfinder(y0,a,t,tolerence,yerror)
        z_array.append(zfinder(yerror,a,t))
        a = a + ppr_step
    a = ppr_step
    y0 = 0.06125 * ( a) * t * (math.exp((-1.2 * ((1 - t)**2)))) + 0.001
    z_array[0] = zfinder(NotCorrected_yfinder(y0,a,t,tolerence),a,t)
    return z_array

def zfinder(y,ppr,t):
    z = (0.06125*ppr*t*math.exp(-1.2*(1-t)*(1-t)))/y
    return z

def absolute(abs):
    if abs < 0:
        return -abs
    else:
        return abs

def yfinder(y,ppr,t,tolerence,yerror):
    ystock = y
    i = 1
    y1 = y - (F(y,ppr,t)/Fprime(y,ppr,t))
    while abs(y1-y) > tolerence:
        if y1 < 0 or (abs(yerror-y1)>0.1 and ppr != 0):
            y1 = ystock/(1.001**i)
            i = i+1
        y = y1
        y1 = y1 - (F(y1,ppr,t)/Fprime(y1,ppr,t))
    return y1
def NotCorrected_yfinder(y,ppr,t,tolerence):
    ystock = y
    i = 1
    y1 = y - (F(y,ppr,t)/Fprime(y,ppr,t))
    while abs(y1-y) > tolerence:
        if y1 < 0:
            y1 = ystock/(1.001**i)
            i = i+1
        y = y1
        y1 = y1 - (F(y1,ppr,t)/Fprime(y1,ppr,t))
    return y1

def F(y,ppr,t):
    F1 = -0.06125*ppr*t*math.exp(-1.2*(1-t)*(1-t))
    F2 = (y+y**2+y**3-y**4)/((1-y)**3)
    F31 = 14.76*t
    F32 = -9.76*t*t
    F33 = 4.58*t*t*t
    F34 = y*y
    F3 = -(F31+F32+F33)*F34
    F41 = (90.7*t-242.2*t*t+42.4*t*t*t)
    F42 = (y**(2.18+2.82*t))
    F4 = F42*F41
    return F1+F2+F3+F4

def Fprime(y,ppr,t):
    s = (1+4*y+4*y*y-4*y*y*y+y*y*y*y)/((1-y)**4) - (29.52*t-19.52*t*t+9.16*t*t*t)*y + (2.18+2.82*t)*(90.7*t-242.2*t*t+42.4*t*t*t)*(y**(1.18+2.82*t))
    return s

def ppr_axis_gene(ppr_step,ppr_max):
    ppr = 0
    ppr_axis = []
    while ppr <= ppr_max:
        ppr_axis.append(ppr)
        ppr = ppr + ppr_step
    return ppr_axis

def tpr_axis_gene(tpr_min,tpr_step,tpr_max):
    tpr_axis = []
    while tpr_min <= tpr_max:
        tpr_axis.append(tpr_min)
        tpr_min = tpr_min + tpr_step
    return tpr_axis

y0 =  584
tolerence = 10**-3
ppr_step = 10**-2
ppr_max = 20
tpr_step = 0.1
tpr_min = 1.05
tpr_max = 3.05
z_list = []
ppr_axis = []
tpr_axis = []


ppr_axis = ppr_axis_gene(ppr_step,ppr_max)
tpr_axis = tpr_axis_gene(tpr_min,tpr_step,tpr_max)
z_list = zlisting(tpr_min,ppr_step,tpr_step,tpr_max,ppr_max,y0,tolerence)
# print(z_list[0])
fig, ax = plt.subplots()
i = 0
while i < ((tpr_max-tpr_min)/tpr_step):
    ax.plot(ppr_axis,z_list[i], label=('Tpr = '+str(tpr_axis[i])[:4]))
    i = i + 1
ax.set_xlabel('Pseudo reduce pressure Ppr')
ax.set_ylabel('z factor')
ax.legend()
plt.show()
