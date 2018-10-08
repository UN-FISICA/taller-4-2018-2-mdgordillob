import math
import numpy as np
from scipy import optimize
class Zeros(object):
	def __init__(self,f,metodo,epsilon=1E-7,N=10000,store=False ):
		self.f=f
		self.epsilon=float(epsilon)
		self.metodo=metodo
		self.store=store
		self.N=N
	def __call__(self,x):
		metodo=self.metodo
		if(metodo=="Newton"):
			print(metodo)
			f=self.f
			N=self.N
			store=self.store
			epsilon=self.epsilon
			f_value = f(x)
			n = 0
			if store: info = [(x, f_value)]
			while abs(f_value) > epsilon and n <= N:
				dfdx_1 = (Derivative(f,"extrapolada"))
				t=dfdx_1(x)
				x = x - f_value/t
				n += 1
				f_value = f(x)
				if (abs(t) < 1E-140):
					raise ValueError("Newton: f'(%g)=%g" % (x,t))
				if store: info.append((x, f_value))
			if store:
				return x, info
			else:
				return x, f_value,n	
		if(metodo=="Bisectriz"):
			print(metodo)
			f=self.f
			epsilon=self.epsilon
			N=self.N
			t=bisect(f,-10,10,N,epsilon)
			return t
		if(metodo=="Interpolacion"):
			print(metodo)
			f=self.f
			epsilon=self.epsilon
			N=self.N
			r=secante(f,-9,-1,N,epsilon)
			return r
		if(metodo=="newton-sp"):
			print(metodo)
			f=self.f
			epsilon=self.epsilon
			N=self.N
			r = optimize.newton(f,x,fprime=Derivative(_g,"extrapolada"))
			return r
		if(metodo=="solve-sp"):
			print(metodo)
			f=self.f
			epsilon=self.epsilon
			N=self.N
			r = optimize.fsolve(f,x)
			return r
		if(metodo=="brentq-sp"):
			print(metodo)
			f=self.f
			epsilon=self.epsilon
			N=self.N
			r=optimize.brentq(f, -N, N, args=(), xtol=epsilon, rtol=epsilon, maxiter=N, full_output=False, disp=True)
			return r
				
   
class Derivative(object):
	def __init__(self, f,metodo,dx=1E-5):
		self.f=f
		self.dx=float(dx)
		self.metodo=metodo
	def __call__(self, x):
		metodo=self.metodo
		if(metodo=="adelante"):
			print(metodo)
			f=self.f
			dx=self.dx
			return (f(x+dx) - f(x))/dx
		if(metodo=="central"):
			print(metodo)
			f=self.f
			dx=self.dx
			return ((f(x+(dx/2))-f(x-(dx/2)))/dx)
		if(metodo=="extrapolada"):
			f=self.f
			dx=self.dx
			f1=(f(x+(dx/2))-f(x-(dx/2)))/dx
			f2=(f(x+dx/4)-f(x-dx/4))/(dx/2)
			return (((4*f2-f1)/3))
		if(metodo=="segundad"):
			print(metodo)
			f=self.f
			dx=self.dx
			return((f(x+dx)+f(x-dx)-2*f(x))/(dx*dx))

		
def _g(x):
    return x**(3)+5	




def samesign(a, b):
	return a * b > 0

def bisect(func, a, b,N,epsilon):
	assert not samesign(func(a), func(b))
	for i in range(N):
		if(func(x)>=epsilon):
			s = (b + a) / 2.0
			if samesign(func(a), func(s)):
				a = s
			else:
				b = s
	return s,func(s),i
def secante(f,a,b,N,epsilon):
	for i in range(N):
		if(f(a)==f(b)):
			return (t,f(t),R)
		t= b-(f(b))*((b - a))/((f(b))-f(a))
		a=b
		b=t
		if(f(t)<epsilon):
			R=i
			i=N
	return (t,f(t),R)
					
if __name__=="__main__":
	print("DERIVADAS")
	print("---------------------------------------")
	print("extrapolada")
	df1 = Derivative(np.sin,"extrapolada")
	df2 = Derivative(np.sin,"central")
	df3 = Derivative(np.sin,"adelante")
	df4 = Derivative(np.sin,"segundad")
	zf1	= Zeros(_g,"Newton",1E-7,1000,False)
	zf2	= Zeros(_g,"Bisectriz",1E-7,10,False)
	zf3 = Zeros(_g,"Interpolacion",1E-7,1000,False)
	zf4 = Zeros(_g,"newton-sp",1E-7,1000,False)
	zf5 = Zeros(_g,"solve-sp",1E-7,1000,False)
	zf6 = Zeros(_g,"brentq-sp",1E-7,1000,False)
	
	y=np.pi
	x=5
	print(df1(y))
	print(df2(y))
	print(df3(y))
	print(df4(y))
	print("CEROS")
	print("---------------------------------------")
	print(zf1(x))
	print(zf2(x))
	print(zf3(x))
	print(zf4(x))
	print(zf5(x))
	print(zf6(x))
	
	
	
