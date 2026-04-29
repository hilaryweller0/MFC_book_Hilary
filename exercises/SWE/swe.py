# Numerical Solution of the Linearised Shallow Water Equations

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass, field

def cosBell(x, xmid, width):
    halfWidth = 0.5*width
    a = xmid - halfWidth
    b = xmid + halfWidth
    return np.where((x >= a) & (x <= b),
                    0.5*(1-np.cos(2*np.pi*(x-a)/width)), 0)

def FB_Cgrid(h, u, H, g, dt, dx):
    h -= H*dt/dx*(u - np.roll(u,1))
    u -= g*dt/dx*(np.roll(h,-1) - h)

def FB_Agrid(h, u, H, g, dt, dx):
    h -= 0.5*H*dt/dx*(np.roll(u,-1) - np.roll(u,1))
    u -= 0.5*g*dt/dx*(np.roll(h,-1) - np.roll(h,1))

def ell2(h, h0):
    return np.sqrt(np.sum((h-h0)**2)/np.sum(h0**2))

@dataclass
class Params:
    L: float = 8e6
    nx: int = 40
    g: float = 10
    H: float = 1e4
    xmid: float = 0.5
    width: float = 1/3
    dx:    float = field(init=False)
    
    def __post_init__(self):
        self.dx = self.L/self.nx
    
    def x(self) -> np.ndarray:
        return np.arange(0,self.L,self.dx)
    
    def xh(self) -> np.ndarray:
        return np.arange(0.5*self.dx,self.L,self.dx)
    
    def init(self, x: np.ndarray) -> np.ndarray:
        return cosBell(x, self.L*self.xmid, self.L*self.width), \
                       np.zeros_like(x)
    
    def Tround(self) -> float:
        "Time for one complete revolution at the gravity wave speed"
        return self.L/np.sqrt(self.g*self.H)
    
    def stable_dt_nt(self) -> tuple[float,int]:
        "Find the largest stable time step that exactly divides Tround"
        dt = self.dx/np.sqrt(self.g*self.H)
        nt = int(np.ceil(self.Tround()/dt))
        dt = self.Tround()/nt
        return dt, nt

p = Params(nx=80)

# Grid and initial conditions
x = p.x()
xu = p.xh()
h0, u0 = p.init(x)
T = p.Tround()
dt,nt = p.stable_dt_nt()
dt *= 0.5; nt *= 2

# Simulate one complete revolution
h = h0.copy(); u = u0.copy()
plt.plot(x, h)
for n in range(nt):
    FB_Cgrid(h, u, p.H, p.g, dt, p.dx)
    if (n+1)%5 == 0:
        plt.plot(x, h)

plt.show()

l2C = ell2(h, h0)

# A-grid version
h = h0.copy(); u = u0.copy()
plt.plot(x, h)
for n in range(nt):
    FB_Agrid(h, u, p.H, p.g, dt, p.dx)
    if (n+1)%5 == 0:
        plt.plot(x, h)

plt.show()

l2A = ell2(h, h0)

