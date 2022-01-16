# Abiy Tasissa
# Original version: 24 November 2013
# The Remez Exchange algorithm
import numpy
from pylab import *
import scipy.linalg
import scipy.optimize
import scipy.misc
import math

# Define the method
# m--degree of the numerator
# n- degree of the denominator
# x- arrays of x-data
# y- arrays of y-data 


def remez(m,n,x,y,tol):
 	# The number of unknowns
	num = n+m+2
    # Determine the length of x-array: num_pts
	num_pts = len(x)
	# Create x and y nodes
	xyindex = np.linspace(0,num_pts-1,num)
	xyindex = [int(i) for i in xyindex]
	xnodes = x[xyindex]
	ynodes = y[xyindex]
	# Create the matrix and find the solution iteratively since the error term is nonlinear
	iteration = 0.0
	max_iterations = 3
	#while True:
	for i in range(max_iterations):
		error = 0.0
		for i in range(max_iterations):
		#while True:
			P = numpy.zeros([num,num])
			# Fill the matrix, first the numerator terms
			P[0:num,0:m+1] = np.vander(xnodes, m+1,increasing = True)
			# Fill the matrix, now the denominator terms
			for i in range(num):
				for j in range(m+1,m+n+1):
					P[i][j]=( (-1**i)*error-ynodes[i])*(xnodes[i]**(j-m))
			# Now do the error term
			P[:,m+n+1] = np.power(-1,np.arange(0,num,1))
			# create the right vector using y
			# solve for coefficients of polynomials and E
			c=scipy.linalg.solve(P,ynodes)
			if math.fabs(c[n+m+1]-error)<=tol:
				break
			error=c[n+m+1]
		# Compute residuals
		res=np.zeros(num_pts)
		for i in range(num_pts):
			res[i]=rational(x[i],c,m,n)-y[i]
		iteration=iteration+1.0
		# check condition and create new xnodes
		if min(res)<=abs(c[m+n+1]) and max(res)<=abs(c[m+n+1]):
			break
		# Create new nodes for next iteration
		for i in range(num-1):
			a=min(res[xyindex[i]:xyindex[i+1]+1])
			b=max(res[xyindex[i]:xyindex[i+1]+1])
			i1=xyindex[i]+nonzero(res[xyindex[i]:xyindex[i+1]+1]==a)[0][0]
			i2=xyindex[i]+nonzero(res[xyindex[i]:xyindex[i+1]+1]==b)[0][0]
			if res[xyindex[i]]*a>0 :
				xyindex[i]=i1
				#print(xyindex[i])
			if res[xyindex[i+1]]*a>0:
				xyindex[i+1]=i1
				#print(xyindex[i+1])
			if res[xyindex[i]]*b>0 :
				xyindex[i]=i2
				#print(xyindex[i])
			if res[xyindex[i+1]]*b>0 :
				xyindex[i+1]=i2
				#print(xyindex[i+1])
			
		xnodes = np.zeros(num)
		ynodes = np.zeros(num)
		xnodes  = x[xyindex]
		ynodes =  y[xyindex] 
	return c


# Define the error function
def rational(x,c,m,n):
	p1=0.0
	for i in range(m+1):
		p1=p1+c[i]*(x**i)
	q1=1.0
	for i in range(m+1,n+m+1):
		q1=q1+c[i]*(x**(i-m))
	return p1/q1



# Test the algorithm
x = np.arange(0.1,1.8,0.1)
y = np.array([16.2750, 30.2720, 41.5290, 51.1910, 60.4920, 69.7630,79.1350, 89.4360, 102.5120, 118.6110, 
	138.3820, 168.2130, 244.4970, 523.4050, 1188.6980, 5022.4910, 11866.0000])
m = 4
n = 2
# Tolerance for the algorithm
tol = 1e-5
remez_sol =remez(m,n,x,y,tol)
print(remez_sol)
remez_sol[1]
# Remez Solution
remez_sol_num=0.0
for i in range(m+1):
	remez_sol_num=remez_sol_num+remez_sol[i]*(x**i)
remez_sol_dnm=1.0
for i in range(m+1,m+n+1):
	remez_sol_dnm=remez_sol_dnm+remez_sol[i]*(x**(i-4))
y_approx=remez_sol_num/remez_sol_dnm
# Visualize results by comparing true solution with numerical solution
plot(x,y,'r--',x,y_approx,'g*')
grid(True)
legend(['Ground solution','Numerical approximation'])
savefig('test_remez.png')
show()



