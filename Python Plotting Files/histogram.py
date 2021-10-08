import matplotlib.pylab as plt
import numpy as np

# Open the file and read everything as a string into 'f'
f = open('scatters.dat','r')

#Separate the lines
lines = f.readlines()[0:]

#Close the file
f.close()

# Initialize lists to add the data to
scatters = []
for line in lines:						#goes over each 'line' in 'lines'
    # Split our line into columns
    data = line.split()
    scatters.append(float(data[0]))

print "minimum scatters=",min(scatters)
print "maximum scatters=",max(scatters)
print "len scatters=",len(scatters)

bins = range(0,16,1)

binsnew = []
i=0
for line in bins:
	binsnew.append(float(bins[i])/2.)
	i = i+1

plt.hist(scatters,binsnew,histtype='stepfilled')
plt.xlabel(r'Number of Scatterings per Photon',fontsize=16)
plt.ylabel(r'Number of Photons',fontsize=16)
plt.title(r'Histogram of Number of Scatterings per Photon'+'\n'+'vs. Number of Photons',fontsize=12)


plt.show()
