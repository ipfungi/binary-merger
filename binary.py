import numpy as np
from scipy import constants as cst

class Inspiral():

    # constructor
    def __init__(self, m1, m2, w, a): # self contient une instance de la classe
        self.m1 = m1
        self.m2 = m2
        self.w = w # initial orbital frequency
        self.a_i = a # distance between m1 and m1

        # protected attributes that cannot be modified manually since they depend on the orbital parameters above    
        self._beta_i = self.w*self.a_i/cst.c
        self._mu = self.m1*self.m2/(self.m1 + self.m2)
        self._mc = ((self.m1*self.m2)**(3/5))/(self.m1 + self.m2)**(1/5)    
        self._tc = (5/256)*((cst.G*(self.m1 + self.m2)**2)/(self._mu*cst.c**3))*(self._beta_i**(-8))   
        print("chirp time", self._tc) 

    # methods useful for the main instance method "inspiral"      
    def f_orbit(self, t):
        return 0.5*((768/(15*cst.c**5))**(-3/8))*np.pi**(-1)*((cst.G*self._mc)**(-5/8))/(self._tc - t)**(3/8)

    def a_orbit(self, f):
        return (cst.G * (self.m1 + self.m2) / (2*np.pi*f)**2)**(1/3)
    
    # main method
    def strain(self, radius1, steps = 10**4): # create OPTIONNAL ARGS !! 
        r1 = np.empty((steps,2))
        r2 = np.empty((steps,2))
        radius2 = np.abs(radius1 - self.a_i)

        Rs1 = 2*cst.G*self.m1/cst.c**2
        rmax = radius1*2
   
        h = np.ones((steps, steps, steps)) # contain the spherical waves for each radius (number == steps) and then for each step

        r, phi = np.linspace(Rs1*10**(-3), rmax, steps), np.linspace(0., 2*np.pi, steps )
        uno = np.ones((steps, steps))
        R, PHI = np.meshgrid(r, phi)

        delta_t = self._tc/steps
        t = 0

        for i in range(steps):   
            f = self.f_orbit(t) # always call a method with self
            a = self.a_orbit(f)
            phi_orbit = 2*np.pi*f*t

            r1[i][0] = radius1*np.cos(phi_orbit)
            r1[i][1] = radius1*np.sin(phi_orbit)
            r2[i][0] = -radius2*np.cos(phi_orbit) 
            r2[i][1] = -radius2*np.sin(phi_orbit)

            v = a*2*np.pi*f # angular velocity for a circular orbit
            f_gw = 2*f
            omega = 2*np.pi*f_gw
            k = omega/cst.c

            factor = -4*cst.G*self._mu*(v**2)/cst.c**4 # attention without r
           
            h[i] = factor*np.cos(omega*t*uno - k*R)
            print(t, "élement", i )

            t += delta_t

            radius1 = self.m2*a/(self.m1 + self.m2)
            radius2 = np.abs(radius1 - a)

        return r1, r2, h, R, PHI


class Ringdown:
    pass

class Binary: # contient la config du binary 
    pass 