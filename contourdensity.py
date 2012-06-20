
import matplotlib.pyplot as plt
from numpy import *
from scipy import *
from pylab import *

data = loadtxt("hist3d.out")
x,y,z = data[:,0], data[:,1], data[:,2]


rmax = 80
delta = y[1]-y[0]
steps = int(2*rmax/delta);
print y[0]
print y[steps]
print delta

 # Creating the grid of coordinates x,y 
x,y = mgrid[-rmax:rmax:delta,-rmax:rmax:delta]

f  = zeros((steps,steps))

for i in range( steps ):
   for j in range( steps):
      f[i][j] = z[i*steps+j]


# Make plot with vertical (default) colorbar
fig = plt.figure()
ax = fig.add_subplot(111)

cax = ax.imshow(f, extent=[-rmax,-rmax + steps*delta, -rmax, -rmax + steps*delta], interpolation='nearest')
ax.set_title(r'N = 12, $\omega$ = 0.01')
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')

#color_start = (int(min(energy)*100)*0.01)
#color_end = (int(max(energy)*100)*0.01)
#color_end = 0.015

#v = linspace(0, color_end, 5, endpoint=True)

# Add colorbar, make sure to specify tick locations to match desired ticklabels
cbar = fig.colorbar(cax)


plt.show()

savefig('test')
