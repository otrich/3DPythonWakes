import scipy as sc
import numpy as np
import funcs

def CGEqns(x,y,u,dX,dY,N,M,n,F):
	''' This function sets up the system we can use to solve '''
	# using the ordering in the paper by McCue 
	eqns = np.zeros((1,2*N*M))
	[zetax, zeta1, phix, phi1] = funcs.reshapingUnknowns(u,M,N)
	phi = funcs.allVals(phi1,phix,dX,M,N)
	zeta = funcs.allVals(zeta1,zetax,dX,M,N)
	phiy = funcs.yDerivs(phi,y,M,N):
	zetay = funcs.yDerivs(zeta,y,M,N):

	phiH = funcs.halfMesh(phi,dX,N)
	zetaH = funcs.halfMesh(zeta,dX,N)
	phixH = funcs.halfMesh(phix,dX,N)
	zetaxH = funcs.halfMesh(zetax,dX,N)
	phiyH = funcs.halfMesh(phiy,dX,N)
	zetayH = funcs.halfMesh(zetay,dX,N)		

	# enforcing the condition at the surface (half mesh)
	for jj in range(M):
		for ii in range(N-1):
			eqns[ii*jj] = 1./2.*((1+etaxH**2.)*phiyH**2.+(1+etayH**2.)*phixH**2.-2.*zetaxH*zetayH*phixH*phiyH)/(1.+phixH**2.+phiyH**2.) + zetaH/F**2.-1./2.

	# enforcing the boundary integral (half mesh)
	for jj in range(M):
		for ii in range(N-1):
			eqns[ii*jj+M*(N-1)] = funcs.integral1 + funcs.integral2 + funcs.integral3 + funcs.integral4
			
	# enforcing the boundary condition
	eqns[(2*NM-2*M):4*M-1] = np.vstack((zeta1, phix-1, zetax[0,:], phix[0,:] - x[0,:]))
	''' want to impose the far field conditions
	1. phi goes like x + C
	2. zeta goes like 0
	3, phix goes like 1
	4. put something on zetax also


	return eqns
