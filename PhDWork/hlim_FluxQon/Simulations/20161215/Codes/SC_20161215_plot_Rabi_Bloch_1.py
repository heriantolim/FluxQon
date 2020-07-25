import numpy as np
from qutip import *

def psi0(theta=0,phi=0):
	theta/=2
	return np.array([np.cos(theta),np.exp(1j*phi)*np.sin(theta)])

def Ry(theta=0):
	theta/=2
	return np.array([[np.cos(theta),-np.sin(theta)],
		[np.sin(theta),np.cos(theta)]])
def Rz(theta=0):
	theta/=2
	return np.array([[np.exp(-1j*theta),0],
		[0,np.exp(1j*theta)]])

def Omega(Delta=0,n=1,g=1):
	return np.sqrt(4*g**2*n+Delta**2)

def alpha(Delta=0,n=1,g=1):
	om=Omega(Delta,n,g)
	return 2*np.arctan(np.sqrt((om-Delta)/(om+Delta)))

def U(t=0,Delta=0,n=1,g=1):
	return Ry(alpha(Delta,n,g)).dot(Rz(Omega(Delta,n,g)*t)).dot(
		Ry(-alpha(Delta,n,g)))

def Rabi(N,psi=psi0(0,0),Delta=0,n=1,g=1):
	t=np.linspace(0,2*np.pi/Omega(Delta,n,g),N)
	x=np.zeros(N)
	y=np.zeros(N)
	z=np.zeros(N)
	for i in np.arange(N):
		psi1=U(t[i]).dot(psi)
		r=np.absolute(psi1)
		a=np.angle(psi1)
		a=a[1]-a[0]
		x[i]=2*r[0]*r[1]*np.cos(a)
		y[i]=2*r[0]*r[1]*np.sin(a)
		z[i]=r[0]**2-r[1]**2
	return [x,y,z]

b=Bloch()
b.zlabel=['$\\vert g,n\\rangle$','$\\vert e,n-1\\rangle$']
b.zlpos=[1.25,-1.35]
b.add_points(Rabi(40,psi0(np.pi/4)))
b.show()
