# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 22:03:53 2018

@author: klt11
"""

import numpy as np
import scipy.special as spe
from math import sqrt
import matplotlib.pyplot as plt


def Q(x):
    return 0.5*spe.erfc(x/np.sqrt(2))


def BER_SNR(M, x):
    r = np.power(10, x/10.0)
    BER_theory = 4.0*(1.0-1.0/np.sqrt(M))*Q(np.sqrt(3.0*r/(2.0*(M-1))))/np.log2(M)
    return BER_theory

x = np.linspace(-5, 25, 300)
xx = np.array(range(5, 20))

simu_data = np.array([[0.1035, 0.0787, 0.0606, 0.0391, 0.0225, 0.0136, 0.0061, 0.0024, 0.0008, 0.0002, 0.0001, 0.0000, 0.0000, 0.0000, 0.0000], 
                      [0.2179, 0.1953, 0.1634, 0.1357, 0.1131, 0.0815, 0.0600, 0.0426, 0.0233, 0.0124, 0.0067, 0.0018, 0.0008, 0.0003, 0.0000], 
                      [0.2376, 0.2124, 0.1878, 0.1653, 0.1418, 0.1191, 0.0975, 0.0775, 0.0600, 0.0425, 0.0296, 0.0174, 0.0096, 0.0045, 0.0020]])

    
QAM8_fit_func = np.polyfit(xx, simu_data[1], 3)
QAM8_fit_func = np.poly1d(QAM8_fit_func)


fig, ax = plt.subplots(figsize=(14,14))
ax.plot(x, BER_SNR(4, x), 'r-', lw=2, label='4QAM-Theory')
#ax.plot(x, BER_SNR(8, x), 'g-', lw=2, label='8QAM-Theory')
ax.plot(x, QAM8_fit_func(x), 'g-', lw=2, label='8QAM-Fit')
ax.plot(x, BER_SNR(16, x), 'b-', lw=2, label='16QAM-Theory')
ax.plot(xx, simu_data[0, :], 'ro', markersize=6, label='4QAM-Simulation')
ax.plot(xx, simu_data[1, :], 'go', markersize=6, label='8QAM-Simulation')
ax.plot(xx, simu_data[2, :], 'bo', markersize=6, label='16QAM-Simulation')
ax.set_xlim(-6, 26)
ax.set_ylim(-0.01, 0.36)
#plt.legend(['4QAM-Theory', '8QAM-Theory', '16QAM-Theory', 
#            '4QAM-Simulation', '8QAM-Simulation', '16QAM-Simulation'], loc = 'upper right')
plt.legend(['4QAM-Theory', '8QAM-Fit', '16QAM-Theory', 
            '4QAM-Simulation', '8QAM-Simulation', '16QAM-Simulation'], loc = 'upper right')
plt.xlabel('SNR (dB)', size=20)
plt.ylabel('BER ', size=20)
plt.show()