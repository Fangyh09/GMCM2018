# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 19:30:38 2018

@author: klt11
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 10:39:18 2018

@author: klt11
"""

import numpy as np
import scipy.special as spe
from math import sqrt
import matplotlib.pyplot as plt


# M = 4 or 16
def generate_square_points(M):
    s = int(np.sqrt(M))
    assert np.log2(M) % 2 == 0
    points = np.zeros((M, 2))
    for i in range(s):
        for j in range(s):
            points[i*s+j, 1] = -(s-1)+2*i
            points[i*s+j, 0] = -(s-1)+2*j
    return points
            

def generate_noise(points, sigma2):
    noise = np.random.normal(0.0, np.sqrt(sigma2), size=points.shape)
    return points + noise


def show_points(points):
    lim = np.abs(points).max() + 1.5
    fig, ax = plt.subplots(figsize=(10,10))
    ax.plot(points[:,0], points[:,1], 'o')
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    plt.show()


def mapping(points, points_with_noise):
    map_list = []
    num_points = points.shape[0]
    for i in range(num_points):
        dist = np.linalg.norm(points_with_noise[i, :] - points, axis=1)
        map_list.append(np.argmin(dist))
    return np.array(map_list)

def cal_signal_error(map_result, signal_prob):
    truth = np.array(range(len(map_result)))
    avg_err = ((map_result != truth) * signal_prob).sum()
    return avg_err

def cal_bit_error(map_result, signal_prob, coding):
    truth = np.array(range(len(map_result)))
    avg_err = ((np.abs(coding[map_result] - coding[truth])).sum(1) * signal_prob).sum()
    return avg_err

def get_signal_avg_power(points, signal_prob):
    return ((np.linalg.norm(points, axis=1)**2) * signal_prob).sum()
    
def Q(x):
    return 0.5*spe.erfc(x/np.sqrt(2))


def get_points(M):
    if M == 4 or M == 16:
        points = generate_square_points(M)
#        points = get_points1()
    elif M == 8:
#        points = np.array([[0.0, 1.0+sqrt(3)], [-1.0, 1.0], [1.0, 1.0], 
#                            [-1.0-sqrt(3), 0.0], [1.0+sqrt(3), 0], 
#                            [-1.0, -1.0], [1.0, -1.0], [0.0, -1.0-sqrt(3)]])
        # 1-2-3-2
        points = np.array([[0, 2*sqrt(3)], [-1, sqrt(3)], [1, sqrt(3)], [-2, 0], 
                            [0, 0], [2, 0], [-1, -sqrt(3)], [1, -sqrt(3)]])
        points[:, 1] = points[:, 1] - np.array(sqrt(3)/4)
#        # 3-2-3
#        points = np.array([[-2, sqrt(3)], [0, sqrt(3)], [2, sqrt(3)], [-1, 0], 
#                            [1, 0], [-2, -sqrt(3)], [0, -sqrt(3)], [2, -sqrt(3)]])
    return points

def get_coding(M):
    if M == 4:
        coding = np.array([[0,0], [0,1], [1,0], [1,1]])
    elif M == 16:
        coding = np.array([[1,0,1,1], [0,0,1,1], [0,1,1,1], [1,1,1,1], 
                           [1,0,0,1], [0,0,0,1], [0,1,0,1], [1,1,0,1], 
                           [1,0,0,0], [0,0,0,0], [0,1,0,0], [1,1,0,0], 
                           [1,0,1,0], [0,0,1,0], [0,1,1,0], [1,1,1,0]])
    elif M == 8:
#        coding = np.array([[1,0,1], [0,0,0], [0,0,1], [1,0,0], 
#                           [1,1,1], [0,1,0], [0,1,1], [1,1,0]])
        # 1-2-3-2
        coding = np.array([[1,1,1], [1,1,0], [0,1,1], [1,0,0], 
                           [0,1,0], [1,0,1], [0,0,0], [0,0,1]])
#        # 3-2-3
#        coding = np.array([[0,0,0], [1,0,0], [1,1,0], [0,0,1], 
#                           [1,0,1], [0,1,0], [0,1,1], [1,1,1]])
    return coding


def cos(x):
    return np.cos(x*np.pi/180)

def sin(x):
    return np.sin(x*np.pi/180)

def test_func1(p):
    p1 = p
    p2 = (1/4-p1)/3
    return -(4*p1*np.log2(p1)+12*p2*np.log2(p2))


def get_points1():
    points = list()
    points.append(np.array([-1, 1]))
    points.append(np.array([1, 1]))
    points.append(np.array([-1, -1]))
    points.append(np.array([1, -1]))
    for x in np.linspace(0, 330, 12):
        points.append(((sqrt(3)+sqrt(35))/4)*np.array([cos(x), sin(x)]))
    return np.asarray(points)


# Start
M = 16
num_iters = 10000

points = get_points(M)
coding = get_coding(M)
# uniform
#signal_prob = np.ones(M) / M
#signal_prob = np.array([0.2, 0.2, 0.2, 0.2, 1/60, 1/60, 1/60, 1/60, 
#                        1/60, 1/60, 1/60, 1/60, 1/60, 1/60, 1/60, 1/60])
#signal_prob = np.array([1/60, 1/60, 1/60, 1/60, 
#                        1/60, 1/5, 1/5, 1/60, 
#                        1/60, 1/5, 1/5, 1/60, 
#                        1/60, 1/60, 1/60, 1/60])
signal_prob = np.array([0.0158, 0.0158, 0.0158, 0.0158, 
                        0.0158, 0.2027, 0.2027, 0.0158, 
                        0.0158, 0.2027, 0.2027, 0.0158, 
                        0.0158, 0.0158, 0.0158, 0.0158])

Ps = get_signal_avg_power(points, signal_prob)

SNR = 13.3
Pn = Ps / (10.0 ** (SNR / 10.0))
#Pn = 0.01
#SNR = 10 * np.log10(Ps/Pn)
EbN0 = 10 * np.log10((Ps/np.log2(M))/(2.0*Pn))
print('Ps: {:.4f}, Pn: {:.4f}, SNR: {:.1f}, EbN0: {:.1f}'.format(Ps, Pn, SNR, EbN0))
    
show_points(points)


signal_err = 0
bit_err = 0
fig, ax = plt.subplots(figsize=(10,10))
for i in range(num_iters):
    points_with_noise = generate_noise(points, Pn)
    ax.plot(points_with_noise[:,0], points_with_noise[:,1], 'ro', markersize=1)
    map_result = mapping(points, points_with_noise)
    signal_err += cal_signal_error(map_result, signal_prob)
    bit_err += cal_bit_error(map_result, signal_prob, coding)

ax.plot(points[:,0], points[:,1], 'o')
lim = np.abs(points).max() + 1.5
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)
ax.set_title('{} QAM'.format(M))
plt.show()
print('SER_simu= {:.6f}, BER_simu = {:.6f}'.format(signal_err/num_iters, bit_err/(np.log2(M)*num_iters)))
#print('SER_simu= {:.6f}, BER_simu = {:.6f}'.format(signal_err/num_iters, bit_err/(3*num_iters)))

if M == 4 or M == 16:
    SER_theory = 4.0*(1.0-1.0/np.sqrt(M))*Q(np.sqrt(3.0*Ps/(2.0*(M-1)*Pn)))
    # Gray Coding
    BER_theory = SER_theory/np.log2(M)
    print('SER_theory = {:.6f}, BER_theory = {:.6f}'.format(SER_theory, BER_theory))