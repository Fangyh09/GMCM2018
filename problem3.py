# -*- coding: utf-8 -*-
"""
Created on Sun Sep 16 22:22:38 2018

@author: klt11
"""

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.optim import lr_scheduler
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


M = 8
num_out = 2*M
num_in = M
d_min = torch.tensor([1.9999])

num_iters = 1000
step_size = num_iters // 2
hidden = 200

class Net(nn.Module):
    def __init__(self):
        super(Net, self).__init__()
        self.fc1 = nn.Linear(num_in, hidden)
        self.fc2 = nn.Linear(hidden, hidden)
        self.fc3 = nn.Linear(hidden, num_out)
    
    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        x = self.fc3(x)
        return x.squeeze().view(M, 2)
    
def Loss1(x):
    return torch.norm(x)

def Loss2(x):
    loss = 0.0
    for i in range(M):
        for j in range(M):
            if j == i:
                break
            p1 = x[i, :]
            p2 = x[j, :]
            loss += 10 * F.relu(d_min - torch.norm(p1- p2))
    return loss


def showpoint(x):
    fig, ax = plt.subplots(figsize=(6,6))
    a = x.data.cpu().numpy()
    ax.plot(a[:,0], a[:,1], 'o')
    for i in range(x.size(0)):
        circle = mpatches.Circle(a[i,:], 1.0, color = "r")
        ax.add_patch(circle)
    ax.set_xlim(-10, 10)
    ax.set_ylim(-10, 10)
    plt.show()


x = torch.randn(1, num_in)

net = Net()

optimizer = optim.Adam(net.parameters(), lr=0.001)
scheduler = lr_scheduler.StepLR(optimizer, step_size=step_size, gamma=0.1)




for i in range(num_iters):
    if (i+1) % step_size == 0:
        scheduler.step()
        
    optimizer.zero_grad()
#    x = torch.randn(1, num_in)
    out = net(x)
    loss = Loss1(out) + Loss2(out)
    loss.backward()
    optimizer.step()
    
    if (i+1) % 10 == 0:
        print('Iter {}/{}, Loss{:.6f}'.format(i+1, num_iters, loss.item()))
        showpoint(out)