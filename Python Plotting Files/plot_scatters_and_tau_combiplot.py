import matplotlib.pylab as plt
import numpy as np
import math

datsize = 'huge'

print
print "***"
print "Data type = ", datsize, " combi.dat"
print "***"
print


plottitle = 'System Optical Depth (Tau) vs. Average Scatterings (<N_scattering>)' + '\n' + 'Quasi-Uniform Density Sphere' 
plottitle = plottitle + '\n' + 'Huge Data Set'

# Open the file and read everything as a string into 'f'
f = open('tau_scatters_136_' + datsize + '_combined.dat','r')
#f = open('tau_vs_scatters_large_data_set1.dat','r')

#Separate the lines
lines = f.readlines()[0:]

#Close the file
f.close()

# Initialize lists to add the data to
tau136 = []
scatter136 = []
tau_error136 = []
scatter_error136 = []

for line in lines:						#goes over each 'line' in 'lines'
    # Split our line into columns
    data = line.split()
    tau136.append(float(data[0]))
    scatter136.append(float(data[1]))
    tau_error136.append(float(data[2]))
    scatter_error136.append(float(data[3]))

# Open the file and read everything as a string into 'f'
f = open('tau_scatters_1002_' + datsize + '_combined.dat','r')
#f = open('tau_vs_scatters_large_data_set1.dat','r')

#Separate the lines
lines = f.readlines()[0:]

#Close the file
f.close()

# Initialize lists to add the data to
tau1002 = []
scatter1002 = []
tau_error1002 = []
scatter_error1002 = []

for line in lines:						#goes over each 'line' in 'lines'
    # Split our line into columns
    data = line.split()
    tau1002.append(float(data[0]))
    scatter1002.append(float(data[1]))
    tau_error1002.append(float(data[2]))
    scatter_error1002.append(float(data[3]))

# Open the file and read everything as a string into 'f'
f = open('tau_scatters_10659_' + datsize + '_combined.dat','r')
#f = open('tau_vs_scatters_large_data_set1.dat','r')

#Separate the lines
lines = f.readlines()[0:]

#Close the file
f.close()

# Initialize lists to add the data to
tau10659 = []
scatter10659 = []
tau_error10659 = []
scatter_error10659 = []

for line in lines:						#goes over each 'line' in 'lines'
    # Split our line into columns
    data = line.split()
    tau10659.append(float(data[0]))
    scatter10659.append(float(data[1]))
    tau_error10659.append(float(data[2]))
    scatter_error10659.append(float(data[3]))


scatter136_diff = []
i = 0
for line in tau136:
	s = (float(tau136[i]) + (float(tau136[i])**2)/2.0)
	diff = (1.0 - (float(scatter136[i])/s))*100. 
	scatter136_diff.append(float(diff))
	i = i + 1
	
i = 0
minval = 10.**(10)
maxval = -10.**(10)
for line in scatter136_diff:
	val = float(scatter136_diff[i])
	if (val < minval):
		minval = val
		minindex = i
	if (val > maxval):
		maxval = val
		maxindex = i
	i = i + 1

max136 = maxval
min136 = minval

scatter1002_diff = []
i = 0
for line in tau1002:
	s = (float(tau1002[i]) + (float(tau1002[i])**2)/2.0)
	diff = (1.0 - (float(scatter1002[i])/s))*100. 
	scatter1002_diff.append(float(diff))
	i = i + 1
	
i = 0
minval = 10.**(10)
maxval = -10.**(10)
for line in scatter1002_diff:
	val = float(scatter1002_diff[i])
	if (val < minval):
		minval = val
		minindex = i
	if (val > maxval):
		maxval = val
		maxindex = i
	i = i + 1

max1002 = maxval
min1002 = minval

scatter10659_diff = []
i = 0
for line in tau10659:
	s = (float(tau10659[i]) + (float(tau10659[i])**2)/2.0)
	diff = (1.0 - (float(scatter10659[i])/s))*100. 
	scatter10659_diff.append(float(diff))
	i = i + 1
	
i = 0
minval = 10.**(10)
maxval = -10.**(10)
for line in scatter10659_diff:
	val = float(scatter10659_diff[i])
	if (val < minval):
		minval = val
		minindex = i
	if (val > maxval):
		maxval = val
		maxindex = i
	i = i + 1

max10659 = maxval
min10659 = minval


print
print "Lowest Deviation from scatter theory for SPH136 = ", min136, " %"
print "Highest Deviation from scatter theory for SPH136 = ", max136, " %"
print
print "Lowest Deviation from scatter theory for SPH1002 = ", min1002, " %"
print "Highest Deviation from scatter theory for SPH1002 = ", max1002, " %"
print
print "Lowest Deviation from scatter theory for SPH10659 = ", min10659, " %"
print "Highest Deviation from scatter theory for SPH10659 = ", max10659, " %"
print

lst = [max(tau136),max(tau1002),max(tau10659)]
taumax = max(lst)

lst2 = [min(tau136),min(tau1002),min(tau10659)]
taumin = min(lst)

print "taumax =", taumax

taumax = int(math.ceil(taumax))

print "taumax new =", taumax

tau_r = [x/1000. for x in range(0,taumax*1000,1)]

scatter_theory = [(float(tau_r[x]) + (float(tau_r[x])**2)/2.0) for x in range(0,len(tau_r))]

xmin = 0.0#taumin
xmax = taumax
ymin = 0.0

ymax = max(scatter_theory)

#xmin = xmin - xmin*0.1
xmax = xmax + xmax*0.1
#ymin = ymin - ymin*0.1
ymax = ymax + ymax*0.1

# yerr = scatter_error
plt.errorbar(tau136,scatter136,xerr = tau_error136 ,color='cyan', fmt= '.', label='136 Particles - Data', elinewidth =0.5)
plt.plot(tau136,scatter136,color='cyan',label='136 Particles - Fit')
plt.xlabel(r'Tau',fontsize=16)
plt.ylabel(r'<N_scattering>',fontsize=16)
plt.title(plottitle,fontsize=12)
plt.axis([xmin,xmax,ymin,ymax])
plt.gca().set_aspect('auto', adjustable='box')

plt.errorbar(tau1002,scatter1002,xerr = tau_error1002 ,color='red', fmt= '.', label='1002 Particles - Data', elinewidth =0.5)
plt.plot(tau1002,scatter1002,color='red',label='1002 Particles - Fit')
plt.xlabel(r'Tau',fontsize=16)
plt.ylabel(r'<N_scattering>',fontsize=16)
plt.title(plottitle,fontsize=12)
plt.axis([xmin,xmax,ymin,ymax])
plt.gca().set_aspect('auto', adjustable='box')

plt.errorbar(tau10659,scatter10659,xerr = tau_error10659 ,color='blue', fmt= '.', label='10659 Particles - Data', elinewidth =0.5)
plt.plot(tau10659,scatter10659,color='blue',label='10659 Particles - Fit')
plt.xlabel(r'Tau',fontsize=16)
plt.ylabel(r'<N_scattering>',fontsize=16)
plt.title(plottitle,fontsize=12)
plt.axis([xmin,xmax,ymin,ymax])
plt.gca().set_aspect('auto', adjustable='box')


plt.plot(tau_r,scatter_theory,color='black',linewidth=2.5,label='Theory')
#plt.plot(tau_theory,scatter_theory,color='black',label='Theory')
plt.xlabel(r'Tau',fontsize=16)
plt.ylabel(r'<N_scattering>',fontsize=16)
plt.title(plottitle,fontsize=12)
plt.axis([xmin,xmax,ymin,ymax])
plt.gca().set_aspect('auto', adjustable='box')
#plt.axes().set_aspect('equal')


plt.legend(loc=2)

plt.show()
