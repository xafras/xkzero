# -*- coding: utf-8 -*-
#   Xafras xyZw Krazos
#   Ouvert le : Sat Jan 14 03:42:48 2023
#   Fermé le :

#   Init

from xkzero import *
import matplotlib.pyplot as plt
import numpy as np

#   Func

get_f = lambda x : np.log(x)*x
get_fp = lambda x : np.log(x)+1
get_foverfp = lambda x : x*np.log(x)/(1+np.log(x))

#   Exec

EPS = 10e-6

L1,delta,_,N1 = Newton(1.6, get_f, get_fp, get_foverfp, return_list=(True), eps=EPS)
L2,delta,_,N2 = dichotomie(0.5, 1.6, get_f, return_list=(True), eps=EPS)

plt.figure()
plt.plot([n for n in range(0,N1+1)], L1, label='Newton')
plt.plot([n for n in range(0,N2+1)], L2, label='dichotomie')
plt.yscale('log')
plt.minorticks_on()
plt.grid('both')
plt.legend()
plt.show()

#   Fin