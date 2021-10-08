import matplotlib.pylab as plt
import numpy as np

# Open the file and read everything as a string into 'f'
f = open('dist_finder_sph1_mass.dat','r')


lines = f.readlines()[0:]

lines = lines[0:]

#Close the file
f.close()

# Initialize lists to add the data to
x = []
y = []

for line in lines:						#goes over each 'line' in 'lines'
    # Split our line into columns
    data = line.split()
    x.append(float(data[0]))
    y.append(float(data[1]))


plt.title('SPH Particles mass (m)'+'\n'+'vs. Photon Packet Travel Distance (L) '+'\n'+'No. SPH particles = 1')#1st Order')#2nd Order + lnm1 System')
#Optical Distance (d)
#SPH Particles mass (m1=m2=m3=m)
#Impact parameter (c)					
plt.xlabel('m [code units]')
plt.ylabel('L [code units]')
#d [dimensionless]
#m [code units]
#c [code units]

ymax = max(y)
ymax = ymax + 0.1*ymax
ymin = min(y)
xmax = max(x)
xmax = xmax + 0.1*xmax
xmin = min(x)

plt.plot(x,y,color='black')

# Add a legend without the frame
plt.legend(frameon=False,loc='lower left')
plt.axis([xmin,xmax,ymin,ymax])

# save the plot as a PNG image
plt.savefig("test.png")

plt.show()
