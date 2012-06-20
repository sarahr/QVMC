import matplotlib.pyplot as plt
from numpy import *
from scipy import *
from pylab import *

# Input parameters
numfiles = 4
rmax = 80
numbins = 80

binwidth = rmax/numbins
counts = 0
density =  zeros((2*numbins,2*numbins))

for i in range(numfiles):
       data = loadtxt("position%i.dat"%(i))
       x,y = data[:,0], data[:,1]

       delta = y[1]-y[0]
       steps = 2*rmax/delta
    
       for j in range(len(x)):
                posx = x[j]
                posy = y[j]

                assignx = floor(posx/binwidth)+ numbins
                assigny = floor(posy/binwidth)+ numbins
                
                if((assignx < 2*numbins) & (assigny < 2*numbins) & (assignx>0) & (assigny>0)):
                	density[assignx][assigny] += 1
                	counts += 1;
 

density /= counts;

 # Creating the grid of coordinates x,y 
x,y = mgrid[-rmax:rmax:delta,-rmax:rmax:delta]

# Make plot with vertical (default) colorbar
fig = plt.figure()
ax = fig.add_subplot(111)

cax = ax.imshow(density, extent=[-rmax,-rmax + steps*delta, -rmax, -rmax + steps*delta], interpolation='nearest')
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



  


