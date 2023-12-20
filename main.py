import binary as bnr 
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm, ticker
import scipy.constants as cst

plt.style.use('dark_background')

# INPUT
M_dot = 2*10**30 # 1 solar mass
m1 = M_dot*10 
m2 = M_dot*10
fi = 20
w = 2*np.pi*fi # initial orbital frequency of 20 Hz

D = (2*cst.G*(m1+m2)/cst.c**2)*4 # distance between each extremity of the black hole if they are aligned and sticked each other
radius1=D*2
a = radius1*2
steps = 10**2

# CALL THE CLASS
blackholes = bnr.Inspiral(m1, m2, w, a) # instance de la classe "Binary"
r1, r2, h, R, PHI = blackholes.strain(radius1, steps)

#OUTPUT 

fig, ax = plt.subplots()

Rs1 = 2*cst.G*m1/cst.c**2
Rs2 = 2*cst.G*m1/cst.c**2
rmax = radius1*2 # distance from the source

X = R*np.cos(PHI)
Y = R*np.sin(PHI)
save_results_to = '/Users/johndoe/Desktop/Code/Inspiral/Drafts'

for i in range(steps):
      cs = ax.contourf(X, Y, h[i], levels=np.linspace(h.min(),h.max(),steps), cmap=cm.winter)
      #cs = ax.imshow(h[7], cmap='gray')
      #cbar = fig.colorbar(cs)

      ax.scatter(r1[i][0], r1[i][1], color = 'black', marker = ".", s = 5*10**2)
      ax.scatter(r2[i][0], r2[i][1], color = 'black', marker = ".", s = (5*10**2)*(Rs2/Rs1)**4)
      ax.set_xticklabels([])
      ax.set_yticklabels([])
      #ax.set_xlim(-1.2*radius1, 1.2*radius1)
      #ax.set_ylim(-1.2*radius1, 1.2*radius1)   
      ax.axis('off')
      plt.savefig("Drafts/" + str(i) + "_m1" + str(m1/M_dot) + "_m2" + str(m2/M_dot) + "_f" + str(fi) + "image.jpg")
      #plt.show()




# to polar coordinates
#r = np.linspace(Rs1, rmax, steps) #
#phi = np.linspace(10**-2, 2*np.pi, steps)


#x = np.linspace(-rmax, rmax, steps)
#y = np.linspace(-rmax, rmax, steps)

#X, Y = np.meshgrid(x, y)

#cs = ax.contourf(X, Y, h[8], levels=np.linspace(h.min(),h.max(),steps), locator=ticker.LogLocator(), cmap=cm.grey)


#ax.set_xlim(-1.2*radius1, 1.2*radius1)
#ax.set_ylim(-1.2*radius1, 1.2*radius1)


