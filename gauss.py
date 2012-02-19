''' x,A = gaussNodes(m,tol=10e-9)
	Returns nodal abscissas {x} and weights {A} of
	Gauss-Legendre m-point quadrature.
'''
from math import cos,pi
from numpy import zeros

digits = 17
num = 10 + 1

def gaussNodes(m,tol=10e-9):

	def legendre(t,m):
		p0 = 1.0; p1 = t
		for k in range(1,m):
			p = ((2.0*k + 1.0)*t*p1 - k*p0)/(1.0 + k )
			p0 = p1; p1 = p
		dp = m*(p0 - t*p1)/(1.0 - t**2)
		return p,dp

	A = zeros(m)
	x = zeros(m)
	nRoots = (m + 1)/2
	for i in range(nRoots):
		t = cos(pi*(i + 0.75)/(m + 0.5))
		for j in range(30):
			p,dp = legendre(t,m)
			dt = -p/dp; t = t + dt
			if abs(dt) < tol:
				x[i] = t; x[m-i-1] = -t
				A[i] = 2.0/(1.0 - t**2)/(dp**2) # Eq.(6.25)
				A[m-i-1] = A[i]
				break
	return x,A

x, w = gaussNodes(num, tol=eval("10e-%d" % digits))

points = ', '.join(("%%.%df" % digits) % x for x in x)
weights = ', '.join(("%%.%df" % digits) % x for x in w)

print """package main

const N = %d

var (
	lob_pts = [...]float64{%s}
	lob_wts = [...]float64{%s}
)""" % (num, points, weights)