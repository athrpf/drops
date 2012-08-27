"""
This script generates a bunch of files for testing the ALE code

It uses
h(x,t) = y0 + a0*sin(k*x + omega*t)
as the height function. The relationship
q_t + h_x = 0
yields
q_t = -a0*kcos(k*x + omega*t)
and therefore
q(x,t) = q0 - a0*k/omega*sin(k*x+omega*t)
"""

import numpy as np

def genfunh(x, t, y0=1.0, a0=0.1, k=1.0, omega=1.0):
    retvalh = y0 + a0*np.sin(k*x + omega*t)
    return retvalh

def genfunq(x, t, q0=1.0, a0=0.1, k=1.0, omega=1.0):
    retvalq = q0 + a0*k/omega*np.sin(k*x+omega*t)
    return retvalq

def generate(startn=161, stopn=200):
    t0 = 0.0
    tinc = 1.0
    nx = 100
    for k in range(startn, stopn):
        hf = open('longh'+str(k) +'.dat', 'w')
        qf = open('longq'+str(k) +'.dat', 'w')
        for j in range(nx):
            hf.write('{}'.format(genfunh(float(j), float(k))))
            qf.write('{}'.format(genfunq(float(j), float(k))))
            if j is not nx-1:
                hf.write('\n')
                qf.write('\n')
        hf.close()
        qf.close()

if __name__=='__main__':
    generate()
