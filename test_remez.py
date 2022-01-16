import numpy
from pylab import *
import scipy.linalg
import scipy.optimize
import scipy.misc
import math

x=np.linspace(0,1,100)
y= sqrt(x)
m = 4
n = 2

c=np.arange(0,10,1)
print(c)
error =0.0
# The number of unknowns
num = n+m+2
print(np.power(3,c))
# Determine the length of x-array: num_pts
num_pts = len(x)
# Create x and y nodes
xyindex = np.linspace(0,num_pts-1,num)
xyindex = [int(i) for i in xyindex]
xnodes = np.zeros(num)
ynodes = np.zeros(num)
xnodes = x[xyindex]
ynodes = y[xyindex]
#print(ynodes)
for i in range(num):
	xnodes[i] = x[xyindex[i]]
	ynodes[i] = y[xyindex[i]]
#print(ynodes)
# Create the matrix and find the solution iteratively since the error term is nonlinear
P = numpy.zeros([num,num])
P1 = numpy.zeros([num,num])
# Fill the matrix, first the numerator terms
P1[0:num,0:m+1] = np.vander(xnodes, m+1,increasing = True)
#print(P1)
for i in range(num):
	for j in range(m+1):
			P[i][j]=(xnodes[i]**j)
# Fill the matrix, now the denominator terms
for i in range(num):
	for j in range(m+1,m+n+1):
		P[i][j]=( (-1**i)*error-ynodes[i])*(xnodes[i]**(j-m))
for i in range(num):
	for j in range(m+1,m+n+1):
		P1[i][j]=( (-1**i)*error-ynodes[i])*(xnodes[i]**(j-m))
# Now do the error term
P1[:,m+n+1] = np.power(-1,np.arange(0,num,1))
for i in range(num):
	P[i][m+n+1]=(-1)**(i)
for i in range(num):
	P1[i][m+n+1]=(-1)**(i)
print(P)
print(P1)
print(np.linalg.norm(P-P1))
# create the right vector using y
# solve for coefficients of polynomials and E
c=scipy.linalg.solve(P,ynodes)
c=scipy.linalg.solve(P1,ynodes)