# make histograms of particle locations
import numpy as np
import matplotlib.pyplot as plt

with open( 'prtclLocation.dat', 'r') as fo:
    data = np.genfromtxt(fo);
    nt = list(set(data[:,3]))

nt.sort()
x = [ [] for i in nt ]
y = [ [] for i in nt ]
z = [ [] for i in nt ]
# print nt
# print x, y, z


# each index of x,y,z corresponds to a specific time-step
# x[i] contains the x-coordinate of all particles at time-step i
for row in data:
    i = nt.index(row[3])
    # print i, row
    x[i].append(row[0])
    y[i].append(row[1])
    z[i].append(row[2])

plt.figure()

plt.subplot(311)
for xt in x:
    bins, edges = np.histogram(xt, 50, normed=0)
    left,right = edges[:-1],edges[1:]
    X = np.array([left,right]).T.flatten()
    Y = np.array([bins,bins]).T.flatten()
    plt.plot( X, Y, label=str(nt[x.index(xt)]))
plt.ylabel('x dist.')
plt.title('Unnormalized Distribution of Particle Locations w/ in-flux')

plt.subplot(312)
for yt in y:
    bins, edges = np.histogram(yt, 50, normed=0)
    left,right = edges[:-1],edges[1:]
    X = np.array([left,right]).T.flatten()
    Y = np.array([bins,bins]).T.flatten()
    plt.plot( X, Y)
plt.ylabel('y dist.')

plt.subplot(313)
for zt in z:
    bins, edges = np.histogram(zt, 50, normed=0)
    left,right = edges[:-1],edges[1:]
    X = np.array([left,right]).T.flatten()
    Y = np.array([bins,bins]).T.flatten()
    plt.plot( X, Y, label=str(nt[z.index(zt)]))
plt.ylabel('z dist.')
plt.legend( fontsize=10)

plt.savefig('prtcl_location_dist3.png', bbox_inches='tight')
plt.show()
