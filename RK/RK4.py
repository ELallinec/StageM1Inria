# Imports
import numpy as np
import matplotlib.pyplot as plt

# Functions


def k1(f, tn, yn):
    return f(tn, yn)


def k2(f, tn, yn, dt):
    return f(tn+dt/2, yn+dt*k1(f, tn, yn)/2)


def k3(f, tn, yn, dt):
    return f(tn+dt/2, yn+dt*k2(f, tn, yn, dt)/2)


def k4(f, tn, yn, dt):
    return f(tn+dt, yn+dt*k3(f, tn, yn, dt))

# f1 = J^-1 Grad(H(p,q)) avec H(p,q) = 1/2(p^2+q^2)


def small_RK4(f, t0, y0, dt):
    return y0 + 1/6 * dt * (k1(f, t0+dt, y0) + 2 * k2(f, t0+dt, y0, dt) + 2 * k3(f, t0+dt, y0, dt) + k4(f, t0+dt, y0, dt))


def tab_RK4(f, t0, T, y0, dt):
    y = y0
    m = round((T-t0)/dt)
    tab_y = np.zeros((len(y0), m+1))
    tab_y[:, 0] = y0
    for k in range(m):
        y = small_RK4(f, t0+k*dt, y, dt)
        tab_y[:, k+1] = y
    return tab_y


def RK4(f, t0, T, y0, dt):
    y = y0
    m = round((T-t0)/dt)
    for k in range(m):
        y = small_RK4(f, t0 + k*dt, y, dt)
    return y
