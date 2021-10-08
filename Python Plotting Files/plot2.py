import matplotlib.pylab as plt
import numpy as np

# Open the file and read everything as a string into 'f'
f = open('output2.dat','r')


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

xyarray = np.column_stack((x,y))

def getKey(item):
    return item[0]

xyarray = sorted(xyarray, key=getKey)

xnew = []
ynew = []

i=0
for line in lines:						#goes over each 'line' in 'lines'
    # Split our line into columns
    data = line.split()
    xnew.append(float(xyarray[i][0]))
    ynew.append(float(xyarray[i][1]))
    i = i + 1 

plt.plot(xnew,ynew,color='black')

plt.title('tau vs l')#1st Order')#2nd Order + lnm1 System')							
plt.xlabel('tau')
plt.ylabel('l')

ymax = max(ynew)
ymin = min(ynew)
xmax = max(xnew)
xmin = min(xnew)

# Add a legend without the frame
plt.legend(frameon=False,loc='lower left')
plt.axis([xmin,xmax,ymin,ymax])

# save the plot as a PNG image
plt.savefig("test.png")

plt.show()
