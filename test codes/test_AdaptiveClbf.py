import osqp
import scipy.sparse as sparse
import scipy.linalg as sl
import numpy as np
import copy
from progress.bar import Bar
import time
import GPy

class Lyapunov():
	def __init__(self,dim,epsilon=1.0):
		self.dim=dim
		self.epsilon = epsilon

	def V(self,x):
		# overload
		return None

	def dV(self,x):
		# overload
		return None

	def signed_angle_dist(self,a,b):
		d = a - b
		if d > np.pi:
			d = d - np.pi*2
		elif d < -np.pi:
			d = d + np.pi*2

		return d

class Barrier():
	def __init__(self,dim,gamma):
		self.dim=dim
		self.gamma = gamma

	def h(self,x):
		# overload
		return None

	def dh(self,x):
		# overload
		return None

	def B(self,x):
		# TODO:  parameterize this, and craft something better.
		hx = self.h(x)
		# return np.exp(-hx+1)
		if hx == 0:
			hx = 1e-3
		return 1.0/hx

	def dB(self,x):
		hx = self.h(x)
		if hx == 0:
			hx = 1e-3
		return -1.0/(hx**2)*self.dh(x)
		# return -np.exp(-hx+1)*self.dh(x)

	def d2B(self,x):
		hx = self.h(x)
		if hx == 0:
			hx = 1e-3
		# return -1.0/(hx**3)*(self.d2h(x) -2*np.matmul(self.dh(x),self.dh(x).T))
		# can ignore d2h because taking trace(G*sig*sig^T*G^T d2B).
		dh = self.dh(x)[2:].T
		return 1.0/(hx**3)*(2*np.outer(dh,dh))

	def get_B_derivatives(self,x):
		hx = self.h(x)
		if hx == 0:
			hx = 1e-3

		dh = self.dh(x)
		return hx, -1.0/(hx*hx)*dh.T, 2.0/(hx*hx*hx)*(np.outer(dh[2:],dh[2:]))


class BarrierAckermannVelocity(Barrier):
	def __init__(self,dim=4,gamma=1.0,bound_from_above=True, v_lim = 1.0):
		self.bound_from_above = bound_from_above
		self.v_lim = v_lim

		Barrier.__init__(self,dim,gamma)

	def h(self,x):
		sgn = 1.0
		if self.bound_from_above:
			sgn = -1.0
		return sgn  * (x[3,:] - self.v_lim)

	def dh(self,x):
		sgn = 1.0
		if self.bound_from_above:
			sgn = -1.0
		return np.array([[0],[0],[0],[sgn]])

	def d2h(self,x):
		return np.array([[0],[0],[0],[0]])

class BarrierAckermannPosition(Barrier):
	def __init__(self,dim=4,gamma=1.0, bound_from_above = True, bound_x = True, pos_lim = 1.0, gamma_p = 1.0):
		self.bound_from_above = bound_from_above
		self.bound_x = bound_x
		self.pos_lim = pos_lim
		self.gamma_p = gamma_p

		Barrier.__init__(self,dim,gamma)

	def h(self,x):
		sgn = 1.0
		if self.bound_from_above:
			sgn = -1.0

		if self.bound_x:
			return sgn * (self.gamma_p * (x[0,:] - self.pos_lim) + x[3,:] * np.cos(x[2,:]))
		else:
			return sgn * (self.gamma_p * (x[1,:] - self.pos_lim) + x[3,:] * np.sin(x[2,:]))
		 
	def dh(self,x):
		sgn = 1.0
		if self.bound_from_above:
			sgn = -1.0

		if self.bound_x:
			return sgn * np.stack((np.array([self.gamma_p]),np.array([0]),-x[3,:] * np.sin(x[2,:]),np.cos(x[2,:])),axis=0)
		else:
			return sgn * np.stack((np.array([0]),np.array([self.gamma_p]),x[3,:] * np.cos(x[2,:]),np.sin(x[2,:])),axis=0)

class BarrierAckermannPositionZ(Barrier):
	def __init__(self,dim=4,gamma=1.0, bound_from_above = True, bound_x = True, pos_lim = 1.0, gamma_p = 1.0):
		self.bound_from_above = bound_from_above
		self.bound_x = bound_x
		self.pos_lim = pos_lim
		self.gamma_p = gamma_p

		Barrier.__init__(self,dim,gamma)

	def h(self,z):
		sgn = 1.0
		if self.bound_from_above:
			sgn = -1.0

		if self.bound_x:
			return sgn * (self.gamma_p * (z[0,:] - self.pos_lim) + z[2,:])
		else:
			return sgn * (self.gamma_p * (z[1,:] - self.pos_lim) + z[3,:])
		 
	def dh(self,z):
		sgn = 1.0
		if self.bound_from_above:
			sgn = -1.0

		if self.bound_x:
			return np.array([[sgn*self.gamma_p],[0],[sgn],[0]])
		else:
			return np.array([[0],[sgn*self.gamma_p],[0],[sgn]])

class BarrierAckermannVelocityZ(Barrier):
	def __init__(self,dim=4,gamma=1.0,bound_from_above=True, v_lim = 1.0):
		self.bound_from_above = bound_from_above
		self.v_lim = v_lim

		Barrier.__init__(self,dim,gamma)

	def h(self,z):
		sgn = 1.0
		if self.bound_from_above:
			sgn = -1.0
		return sgn  * (np.sqrt(z[2]**2 + z[3]**2) * z[4] - self.v_lim)

	def dh(self,z):
		sgn = 1.0
		if self.bound_from_above:
			sgn = -1.0
		z1 = z[0,:]
		z2 = z[1,:]
		z3 = z[2,:]
		z4 = z[3,:]
		z5 = z[4,:]
		v = (np.sqrt(z3**2 + z4**2) + 1.0e-6)
		return np.stack((np.array([0]), np.array([0]), (z3*z5)/v, (z4*z5)/v)) * sgn

	def d2h(self,z):
		sgn = 1.0
		if self.bound_from_above:
			sgn = -1.0
		z1 = z[0,:]
		z2 = z[1,:]
		z3 = z[2,:]
		z4 = z[3,:]
		z5 = z[4,:]
		v = (np.sqrt(z3**2 + z4**2) + 1.0e-6)
		H = np.block([[(z4**2*z5)/(v**3), -(z3*z4*z5)/(v**3)],
			[-(z3*z4*z5)/(v**3),   (z3**2*z5)/(v**3)]])
		return np.block([[np.zeros((2,2)),np.zeros((2,2))],[np.zeros((2,2)),H]], dtype=np.float32) * sgn

class BarrierAckermannPointZ(Barrier):
	def __init__(self,dim=4,gamma=1.0, x=0.0, y=0.0, radius = 1.0, gamma_p = 1.0):
		self.x = x
		self.y = y
		self.radius = radius
		self.gamma_p = gamma_p

		Barrier.__init__(self,dim,gamma)

	def h(self,z):
		sgn = 1.0
		if self.radius < 0.0:
			sgn = -1.0

		d = np.sqrt((z[0] - self.x)*(z[0] - self.x) + (z[1] - self.y)*(z[1] - self.y)) + 1.0e-6
		return sgn * (self.gamma_p * (d - self.radius) + (z[0] - self.x) / d * z[2] + (z[1] - self.y) / d * z[3])
		 
	def dh(self,z):
		sgn = 1.0
		if self.radius < 0.0:
			sgn = -1.0

		d = np.sqrt((z[0] - self.x)*(z[0] - self.x) + (z[1] - self.y)*(z[1] - self.y)) + 1.0e-6
		d_2 = d*d
		d_3 = d*d*d
		y_pos_m_z2 = (self.y - z[1])
		x_pos_m_z1 = (self.x - z[0])
		return sgn * np.array((
			(z[2]*(d_2 - x_pos_m_z1*x_pos_m_z1) - self.gamma_p * d_2 *x_pos_m_z1 - z[3]*x_pos_m_z1*y_pos_m_z2)/d_3,
			(z[3]*(d_2 - y_pos_m_z2*y_pos_m_z2) - self.gamma_p * d_2 *y_pos_m_z2 - z[2]*x_pos_m_z1*y_pos_m_z2)/d_3,
			-x_pos_m_z1/d,
			-y_pos_m_z2/d), dtype=np.float32)

	def d2h(self,z):
		sgn = 1.0
		if self.radius < 0.0:
			sgn = -1.0

		d = np.sqrt((z[0,:] - self.x)**2 + (z[1,:] - self.y)**2) + 1.0e-6
		z1 = z[0,:]
		z2 = z[1,:]
		z3 = z[2,:]
		z4 = z[3,:]
		z5 = z[4,:]
		gamma_p = self.gamma_p
		x_pos = self.x
		y_pos = self.y
		y_pos_m_z2 = (y_pos - z2)
		x_pos_m_z1 = (x_pos - z1)
		d2 = d**2
		d3 = d**3
		d5 = d**5
		z1_2 = z1**2
		z2_2 = z2**2
		x_pos_2 = x_pos**2
		y_pos_2 = y_pos**2
		a11 = (y_pos_m_z2*(gamma_p *(x_pos_2*y_pos - x_pos_2*z2 + 3*y_pos*z2_2 - 2*x_pos*y_pos*z1 + 2*x_pos*z1*z2 + y_pos**3 - 3*y_pos_2*z2 + y_pos*z1_2 - z1_2*z2 - z2**3) - 2*z4*x_pos_2 + 3*z3*x_pos*y_pos + 4*z4*x_pos*z1 - 3*z3*x_pos*z2 + z4*y_pos_2 - 3*z3*y_pos*z1 - 2*z4*y_pos*z2 - 2*z4*z1_2 + 3*z3*z1*z2 + z4*z2_2))/(d5)
		a12 = x_pos_m_z1 * y_pos_m_z2 * (z4 + z3 - gamma_p - (3*z3*x_pos_m_z1)/(d2) - (3*z4*y_pos_m_z2)/(d2))/(d3)
		a13 = y_pos_m_z2**2/(d3)
		a14 = -(x_pos_m_z1*y_pos_m_z2)/(d3)
		a22 = (x_pos_m_z1*(gamma_p*(x_pos**3 - 3*x_pos_2*z1 + x_pos*y_pos_2 - 2*x_pos*y_pos*z2 + 3*x_pos*z1_2 + x_pos*z2_2 - y_pos_2*z1  + 2*y_pos*z1*z2 - z1**3 - z1*z2_2)+ z3*x_pos_2 + 3*z4*x_pos*y_pos - 2*z3*x_pos*z1 - 3*z4*x_pos*z2 - 2*z3*y_pos_2 - 3*z4*y_pos*z1 + 4*z3*y_pos*z2 + z3*z1_2 + 3*z4*z1*z2 - 2*z3*z2_2))/(d5)
		a23 = -(y_pos_m_z2*x_pos_m_z1)/(d3)
		a24 = x_pos_m_z1**2/(d3)

		return sgn * np.block([
			[ a11, a12, a13, a14],
			[ a12, a22, a23, a24],
			[ a13, a23, 0, 0],
			[ a14, a24, 0, 0]], dtype=np.float32)

class LyapunovAckermann(Lyapunov):
	def __init__(self,dim=4,w1=1.0,w2=1.0,w3=1.0,epsilon=1e-6):
		self.w1=w1
		self.w2=w2
		self.w3=w3
		self.epsilon=epsilon
		Lyapunov.__init__(self,dim)

	def convert_to_polar(self,x,x_d):
		if x[3,:] < 0:
			e = -np.sqrt((x_d[0,:]-x[0,:])**2 + (x_d[1,:]-x[1,:])**2)
			phi_t = np.arctan2(x[1,:]-x_d[1,:],x[0,:]-x_d[0,:])
		else:
			e = np.sqrt((x_d[0,:]-x[0,:])**2 + (x_d[1,:]-x[1,:])**2)
			phi_t = np.arctan2(x_d[1,:]-x[1,:],x_d[0,:]-x[0,:])
		alpha = self.signed_angle_dist(phi_t,x[2,:])
		theta = self.signed_angle_dist(phi_t,x_d[2,:])

		# if alpha and theta are opposite signs, we have two confliciting headings.
		# Flip the sign of one to turn in one direction.
		# print('before:',alpha,theta)
		if np.sign(alpha) != np.sign(theta):
			if np.abs(alpha) >= np.abs(theta):
				alpha = alpha - np.sign(alpha) * 2 * np.pi
			else:
				theta = theta - np.sign(theta) * 2 * np.pi
		# print('after:',alpha,theta)
		# print(e,x[2,:]/(np.pi)*180,phi_t/(np.pi)*180,x_d[2,:]/(np.pi)*180,alpha/(np.pi)*180,theta/(np.pi)*180)

		return e,phi_t,alpha,theta

	def V(self,x,x_d):
		e,phi_t,alpha,theta = self.convert_to_polar(x,x_d)
		return 0.5 * (self.w1*alpha**2 + self.w2*theta**2 + self.w3*(x[3,:]-x_d[3,:])**2)

	def dV(self,x,x_d):
		e,phi_t,alpha,theta = self.convert_to_polar(x,x_d)
		c_tmp = (self.w1*alpha + self.w2*theta) / (e + self.epsilon)
		return np.stack((c_tmp*np.sin(phi_t),-c_tmp*np.cos(phi_t),-self.w1*alpha,self.w3*(x[3,:]-x_d[3,:])),axis=0)

class LyapunovAckermannZ(Lyapunov):
	def __init__(self,dim=4,w1=1.0, w2=1.0, w3=1.0, epsilon = 1e-6):
		self.w1=w1
		self.w2=w2
		self.w3=w3
		self.epsilon = epsilon
		Lyapunov.__init__(self,dim)


	def convert_to_polar(self,z,z_d):
		# if z[3,:] < 0:
		# 	e = -np.sqrt((z_d[0,:]-z[0,:])**2 + (z_d[1,:]-z[1,:])**2)
		# 	phi_t = np.arctan2(z[1,:]-z_d[1,:],z[0,:]-z_d[0,:])
		# else:
		e = np.sqrt((z_d[0,:]-z[0,:])**2 + (z_d[1,:]-z[1,:])**2)
		phi_t = np.arctan2(z_d[1,:]-z[1,:],z_d[0,:]-z[0,:])

		phi = np.arctan2(z[3,:],z[2,:])
		phi_d = np.arctan2(z_d[3,:],z_d[2,:])

		alpha = self.signed_angle_dist(phi_t,phi)
		theta = self.signed_angle_dist(phi_t,phi_d)
		# print('angles: ',phi,phi_t,phi_d,alpha,theta)
		# if alpha and theta are opposite signs, we have two confliciting headings.
		# Flip the sign of one to turn in one direction.
		# print('before:',alpha,theta)
		# if np.sign(alpha) != np.sign(theta):
		# 	if np.abs(alpha) >= np.abs(theta):
		# 		alpha = alpha + np.sign(alpha) * 2 * np.pi
		# 	else:
		# 		theta = theta + np.sign(theta) * 2 * np.pi
		# print('after:',alpha,theta)
		# print(e,z[2,:]/(np.pi)*180,phi_t/(np.pi)*180,z_d[2,:]/(np.pi)*180,alpha/(np.pi)*180,theta/(np.pi)*180)

		return e,phi_t,alpha,theta

	def V(self,z,z_d):
		e,phi_t,alpha,theta = self.convert_to_polar(z,z_d)
		v = np.sqrt(z[2,:]**2 + z[3,:]**2)
		v_d = np.sqrt(z_d[2,:]**2 + z_d[3,:]**2)
		return 0.5 * (self.w1*alpha**2 + self.w2*theta**2 + self.w3*(v - v_d)**2)

	def dV(self,z,z_d):
		e,phi_t,alpha,theta = self.convert_to_polar(z,z_d)
		c_tmp = (self.w1*alpha + self.w2*theta) / (e + self.epsilon)
		v = np.sqrt(z[2,:]**2 + z[3,:]**2)
		v_d = np.sqrt(z_d[2,:]**2 + z_d[3,:]**2)

		return np.stack((c_tmp*np.sin(phi_t),
						 -c_tmp*np.cos(phi_t),
						 self.w1*alpha*z[3,:]/(v**2) + self.w3*(v-v_d)*z[2,:]/v,
						 -self.w1*alpha*z[2,:]/(v**2) + self.w3*(v-v_d)*z[3,:]/v),axis=0)

class Dynamics():
	def __init__(self,xdim,udim):
		self.xdim=xdim
		self.udim=udim

	def f(self,x):
		# overload
		return None

	def g(self,x):
		# overload
		return None

	def pseudoinv(self,x):
		return np.matmul(np.linalg.inv(np.matmul(x.T,x)),x.T)

	def convert_z_to_x(self,z):
		v=np.sqrt(z[2,:]**2 + z[3,:]**2) * z[4,:]
		return np.stack((z[0,:],z[1,:],np.arctan2(z[4,:]*z[3,:],z[4,:]*z[2,:]),v),axis=0)

	def convert_x_to_z(self,x):
		v_sign = 1.0
		if x[3,:] < 0:
			v_sign = -1
		return np.stack((x[0,:],x[1,:],x[3,:]*np.cos(x[2,:]),x[3,:]*np.sin(x[2,:]),np.array([v_sign])),axis=0)

class DynamicsAckermann(Dynamics):
	def __init__(self,xdim=4,udim=2):
		Dynamics.__init__(self,xdim=xdim,udim=udim)

	def f(self,x):
		return np.stack((x[3,:] * np.cos(x[2,:]), x[3,:] * np.sin(x[2,:]), np.array([0]), np.array([0])))

	def g(self,x):
		return np.stack((np.array([0,0]),np.array([0,0]),np.append(x[3,:],0),np.array([0,1])))

class DynamicsAckermannModified(Dynamics):
	def __init__(self,xdim=4,udim=2):
		Dynamics.__init__(self,xdim=xdim,udim=udim)

	def f(self,x):
		v_modified = x[3,:]**2 / (2*np.tanh(x[3,:]))
		return np.stack((v_modified * np.cos(x[2,:]), v_modified * np.sin(x[2,:]), x[3,:], -0.5*x[3,:]))

	def g(self,x):
		return np.stack((np.array([0,0]),np.array([0,0]),np.append(x[3,:],0),np.array([0,1])))*1.2


class DynamicsAckermannZ(Dynamics):
	def __init__(self,xdim=4,udim=2,epsilon=1e-6):
		self.epsilon = epsilon
		Dynamics.__init__(self,xdim=xdim,udim=udim)

	def f(self,z):
		return np.stack((np.array([0]),np.array([0]))) * z[4,:]

	def g(self,z):
		v=(np.sqrt(z[2,:]**2 + z[3,:]**2) + self.epsilon) * z[4,:]
		return np.stack((np.concatenate((-z[3,:]*v,z[2,:]/v)),np.concatenate((z[2,:]*v,z[3,:]/v))))


class DynamicsAckermannZModified(Dynamics):
	def __init__(self,xdim=4,udim=2,epsilon=1e-6,disturbance_scale_pos = 0.5,disturbance_scale_vel = 2.0,control_input_scale = 2.0):
		self.epsilon = epsilon
		Dynamics.__init__(self,xdim=xdim,udim=udim)

		self.disturbance_scale_pos = disturbance_scale_pos
		self.disturbance_scale_vel = disturbance_scale_vel
		self.control_input_scale = control_input_scale
		
	def f(self,z):
		v=np.sqrt(z[2,:]**2 + z[3,:]**2) * z[4,:]
		theta = np.arctan2(z[3]*z[4],z[2]*z[4])
		v_disturbance_body = [np.tanh(v**2)*self.disturbance_scale_vel, (0.1+v)*self.disturbance_scale_vel]
		v_disturbance_world = [v_disturbance_body[0] * np.cos(theta) - v_disturbance_body[1] * np.sin(theta),
							   v_disturbance_body[0] * np.sin(theta) + v_disturbance_body[1] * np.cos(theta)]
		return np.array([v_disturbance_world[0], v_disturbance_world[1]])
		
	def g(self,z):
		v=(np.sqrt(z[2,:]**2 + z[3,:]**2)  + self.epsilon) * z[4,:]
		return self.control_input_scale * np.stack((np.concatenate((-z[3,:]*v,z[2,:]/v)),np.concatenate((z[2,:]*v,z[3,:]/v))))

	def step(self,z,u,dt):
		v=np.sqrt(z[2,:]**2 + z[3,:]**2)
		znext = np.zeros((self.xdim+1,1))
		znext[0:2,:] = z[0:2,:] + dt*(z[2:-1,:] + np.array([np.sin(v**3)*self.disturbance_scale_pos, np.cos(-v)*self.disturbance_scale_pos]))
		znext[2:-1,:] = z[2:-1,:] + dt*(self.f(z) + np.matmul(self.g(z),u))
		znext[-1] = z[-1]
		return znext


class QPSolve():
    def __init__(self,dyn,clf,cbf_list,u_lim,u_cost=0.0,u_prev_cost=1.0,p1_cost=1.0e10,p2_cost=1.0e10,verbose=True):
        self.xdim = dyn.xdim
        self.udim = dyn.udim
        self.dyn = dyn
        self.clf = clf
        self.cbf_list = cbf_list
        self.u_lim = u_lim
        self.p1_cost = p1_cost
        self.p2_cost = p2_cost
        self.verbose = verbose
        self.u_cost = u_cost
        self.u_prev_cost = u_prev_cost
        self.K = 0.0
        self.ksig = 1.0
        self.max_var = 1.0
        self.mu_qp_prev = np.zeros((int(self.xdim/2),1), dtype=np.float32)
        self.P = np.eye(self.xdim, dtype=np.float32)
        self.A = np.zeros((self.xdim,self.xdim), dtype=np.float32)
        self.A0 = np.block([[np.zeros((int(self.xdim/2),int(self.xdim/2))),np.eye(int(self.xdim/2))],[np.zeros((int(self.xdim/2),int(self.xdim/2))),np.zeros((int(self.xdim/2),int(self.xdim/2)))]]).astype(np.float32)
        self.G = np.block([[np.zeros((int(self.xdim/2),int(self.xdim/2)))],[np.eye(int(self.xdim/2))]]).astype(np.float32)
        self.res = None
        self.max_error = 1.0

    def update_ricatti(self,A):
        self.A = A
        Q = np.eye(self.xdim, dtype=np.float32)
        self.P = sl.solve_continuous_are(self.A,np.zeros((self.xdim,self.xdim), dtype=np.float32),Q,np.eye(self.xdim, dtype=np.float32))

    def solve(self,x,x_d,mu_d,sigDelta):
        sigDelta = sigDelta * self.ksig
        sigDelta = np.clip(sigDelta,0.0,self.max_var)
        # sigDelta = np.ones((self.xdim/2,1)) * self.max_var # for testing

        # build Q and p matrices to specify minimization expression
        Q = np.diag(np.append(np.append(np.ones(int(self.xdim/2), dtype=np.float32)*(self.u_cost + self.u_prev_cost),self.p1_cost),self.p2_cost))
        self.Q = sparse.csc_matrix(Q)
        self.p = 2*np.append(np.append(-self.mu_qp_prev*self.u_prev_cost,0),0)

        #error dynamics for clf
        e = x[:-1,:]-x_d[:-1,:]
        e = np.clip(e,-self.max_error,self.max_error)
        eTPG = np.matmul(e.T,np.matmul(self.P,self.G))
        G_dyn = np.expand_dims(np.append(np.append(eTPG,1),0),0) #append 1 for clf < d
        Gsig = np.matmul(self.G,sigDelta)
        GssG = np.matmul(Gsig,Gsig.T)
        self.trGssGP = np.trace(np.matmul(GssG,self.P))
        h_dyn = -1 * ( -0.5*np.matmul(e.T,np.matmul(Q,e))
                    + 0.5*np.matmul(e.T,np.matmul(self.P,e)) / self.clf.epsilon
                    + 0.5*self.trGssGP)

        # build constraints for barriers
        N_cbf = len(self.cbf_list)
        G_cbf = np.zeros((N_cbf,int(self.xdim/2)+2), dtype=np.float32)
        h_cbf = np.zeros((N_cbf,1), dtype=np.float32)
        A0x_Gmud = np.matmul(self.A0,x[:-1,:]) + np.matmul(self.G,mu_d)
        GssG_22 = GssG[2:,2:]
        for i, cbf in enumerate(self.cbf_list):
            h_x, dB, d2B = cbf.get_B_derivatives(x)
            G_cbf[i,:] = np.append(np.append(np.einsum('ij,jk',dB,self.G),0),1)
            trGssGd2B = np.einsum('ii',np.einsum('ij,jk',GssG_22,d2B))
            h_cbf[i,:] = -1 * (np.einsum('ij,jk',dB,A0x_Gmud)
                                - cbf.gamma * h_x
                                + 0.5*trGssGd2B)

        # build constraints for control limits
        ginv = np.linalg.inv(self.dyn.g(x))
        l_ctrl = np.matmul(ginv, mu_d - self.dyn.f(x))
        A_ctrl = ginv

        G_ctrl = np.zeros((self.udim*2,int(self.xdim/2)+2), dtype=np.float32)
        h_ctrl = np.zeros((self.udim*2,1), dtype=np.float32)
        for i in range(self.udim):
            G_ctrl[i*2,:int(self.xdim/2)] = - A_ctrl[i,:]
            h_ctrl[i*2] = - self.u_lim[i,0] + l_ctrl[i]
            G_ctrl[i*2+1,:int(self.xdim/2)] = A_ctrl[i,:]
            h_ctrl[i*2+1] = self.u_lim[i,1] - l_ctrl[i]

        # stack into one matrix and vector
        G = np.concatenate((G_dyn,G_cbf,G_ctrl),axis=0)
        h = np.concatenate((h_dyn,h_cbf,h_ctrl),axis=0)

        self.G_csc = sparse.csc_matrix(G)
        self.h = h

        # dummy lower bound
        l = np.ones(h.shape, dtype=np.float32)*np.inf * -1

        #Solve QP
        self.prob = osqp.OSQP()
        exception_called = False
        mu_bar = np.zeros((self.xdim+1), dtype=np.float32)
        # try:
        self.prob.setup(P=self.Q, q=self.p, A=self.G_csc, l=l, u=self.h, verbose=self.verbose)
        self.res = self.prob.solve()
        # except:
            # exception_called = True
        # else:
        mu_bar = self.res.x
        # if exception_called or u_bar[0] is None or np.isnan(u_bar).any():
        if mu_bar[0] is None or np.isnan(mu_bar).any():
            mu_bar = np.zeros((self.xdim+1))
            self.res = None
            print("QP failed!")
        

        self.mu_qp_prev = np.expand_dims(mu_bar[0:int(self.xdim/2)],axis=0).T

        self.V =self.clf.V(x,x_d)

        if self.verbose:
        # if True:
            print('z_ref: ', x.T)
            print('z_des: ', x_d.T)
            print('u_lim', self.u_lim)
            print('V: ', self.V)
            print('Q:', Q)
            print('p:', np.array(self.p))
            print('G_dyn:', G_dyn)
            print('h_dyn:', h_dyn)
            print('trGssGP',self.trGssGP)
            if h_cbf.shape[0] < 10:
                print('G_cbf:', G_cbf)
                print('h_cbf:', h_cbf)
            print('G_ctrl:', G_ctrl)
            print('h_ctrl:', h_ctrl)
            print('result:', mu_bar)
            print('G_all:',self.G_csc)
            print('h_all:',self.h)

        return self.mu_qp_prev
	
import osqp
import scipy.sparse as sparse
import scipy.linalg as sl
import numpy as np

class Lyapunov():
	def __init__(self,dim,epsilon=1.0):
		self.dim=dim
		self.epsilon = epsilon

	def V(self,x):
		# overload
		return None

	def dV(self,x):
		# overload
		return None

	def signed_angle_dist(self,a,b):
		d = a - b
		if d > np.pi:
			d = d - np.pi*2
		elif d < -np.pi:
			d = d + np.pi*2

		return d

class Barrier():
	def __init__(self,dim,gamma):
		self.dim=dim
		self.gamma = gamma

	def h(self,x):
		# overload
		return None

	def dh(self,x):
		# overload
		return None

	def B(self,x):
		# TODO:  parameterize this, and craft something better.
		hx = self.h(x)
		# return np.exp(-hx+1)
		if hx == 0:
			hx = 1e-3
		return 1.0/hx

	def dB(self,x):
		hx = self.h(x)
		if hx == 0:
			hx = 1e-3
		return -1.0/(hx**2)*self.dh(x)
		# return -np.exp(-hx+1)*self.dh(x)

	def d2B(self,x):
		hx = self.h(x)
		if hx == 0:
			hx = 1e-3
		# return -1.0/(hx**3)*(self.d2h(x) -2*np.matmul(self.dh(x),self.dh(x).T))
		# can ignore d2h because taking trace(G*sig*sig^T*G^T d2B).
		dh = self.dh(x)[2:].T
		return 1.0/(hx**3)*(2*np.outer(dh,dh))

	def get_B_derivatives(self,x):
		hx = self.h(x)
		if hx == 0:
			hx = 1e-3

		dh = self.dh(x)
		return hx, -1.0/(hx*hx)*dh.T, 2.0/(hx*hx*hx)*(np.outer(dh[2:],dh[2:]))


class BarrierAckermannVelocity(Barrier):
	def __init__(self,dim=4,gamma=1.0,bound_from_above=True, v_lim = 1.0):
		self.bound_from_above = bound_from_above
		self.v_lim = v_lim

		Barrier.__init__(self,dim,gamma)

	def h(self,x):
		sgn = 1.0
		if self.bound_from_above:
			sgn = -1.0
		return sgn  * (x[3,:] - self.v_lim)

	def dh(self,x):
		sgn = 1.0
		if self.bound_from_above:
			sgn = -1.0
		return np.array([[0],[0],[0],[sgn]])

	def d2h(self,x):
		return np.array([[0],[0],[0],[0]])

class BarrierAckermannPosition(Barrier):
	def __init__(self,dim=4,gamma=1.0, bound_from_above = True, bound_x = True, pos_lim = 1.0, gamma_p = 1.0):
		self.bound_from_above = bound_from_above
		self.bound_x = bound_x
		self.pos_lim = pos_lim
		self.gamma_p = gamma_p

		Barrier.__init__(self,dim,gamma)

	def h(self,x):
		sgn = 1.0
		if self.bound_from_above:
			sgn = -1.0

		if self.bound_x:
			return sgn * (self.gamma_p * (x[0,:] - self.pos_lim) + x[3,:] * np.cos(x[2,:]))
		else:
			return sgn * (self.gamma_p * (x[1,:] - self.pos_lim) + x[3,:] * np.sin(x[2,:]))
		 
	def dh(self,x):
		sgn = 1.0
		if self.bound_from_above:
			sgn = -1.0

		if self.bound_x:
			return sgn * np.stack((np.array([self.gamma_p]),np.array([0]),-x[3,:] * np.sin(x[2,:]),np.cos(x[2,:])),axis=0)
		else:
			return sgn * np.stack((np.array([0]),np.array([self.gamma_p]),x[3,:] * np.cos(x[2,:]),np.sin(x[2,:])),axis=0)

class BarrierAckermannPositionZ(Barrier):
	def __init__(self,dim=4,gamma=1.0, bound_from_above = True, bound_x = True, pos_lim = 1.0, gamma_p = 1.0):
		self.bound_from_above = bound_from_above
		self.bound_x = bound_x
		self.pos_lim = pos_lim
		self.gamma_p = gamma_p

		Barrier.__init__(self,dim,gamma)

	def h(self,z):
		sgn = 1.0
		if self.bound_from_above:
			sgn = -1.0

		if self.bound_x:
			return sgn * (self.gamma_p * (z[0,:] - self.pos_lim) + z[2,:])
		else:
			return sgn * (self.gamma_p * (z[1,:] - self.pos_lim) + z[3,:])
		 
	def dh(self,z):
		sgn = 1.0
		if self.bound_from_above:
			sgn = -1.0

		if self.bound_x:
			return np.array([[sgn*self.gamma_p],[0],[sgn],[0]])
		else:
			return np.array([[0],[sgn*self.gamma_p],[0],[sgn]])

class BarrierAckermannVelocityZ(Barrier):
	def __init__(self,dim=4,gamma=1.0,bound_from_above=True, v_lim = 1.0):
		self.bound_from_above = bound_from_above
		self.v_lim = v_lim

		Barrier.__init__(self,dim,gamma)

	def h(self,z):
		sgn = 1.0
		if self.bound_from_above:
			sgn = -1.0
		return sgn  * (np.sqrt(z[2]**2 + z[3]**2) * z[4] - self.v_lim)

	def dh(self,z):
		sgn = 1.0
		if self.bound_from_above:
			sgn = -1.0
		z1 = z[0,:]
		z2 = z[1,:]
		z3 = z[2,:]
		z4 = z[3,:]
		z5 = z[4,:]
		v = (np.sqrt(z3**2 + z4**2) + 1.0e-6)
		return np.stack((np.array([0]), np.array([0]), (z3*z5)/v, (z4*z5)/v)) * sgn

	def d2h(self,z):
		sgn = 1.0
		if self.bound_from_above:
			sgn = -1.0
		z1 = z[0,:]
		z2 = z[1,:]
		z3 = z[2,:]
		z4 = z[3,:]
		z5 = z[4,:]
		v = (np.sqrt(z3**2 + z4**2) + 1.0e-6)
		H = np.block([[(z4**2*z5)/(v**3), -(z3*z4*z5)/(v**3)],
			[-(z3*z4*z5)/(v**3),   (z3**2*z5)/(v**3)]])
		return np.block([[np.zeros((2,2)),np.zeros((2,2))],[np.zeros((2,2)),H]], dtype=np.float32) * sgn

class BarrierAckermannPointZ(Barrier):
	def __init__(self,dim=4,gamma=1.0, x=0.0, y=0.0, radius = 1.0, gamma_p = 1.0):
		self.x = x
		self.y = y
		self.radius = radius
		self.gamma_p = gamma_p

		Barrier.__init__(self,dim,gamma)

	def h(self,z):
		sgn = 1.0
		if self.radius < 0.0:
			sgn = -1.0

		d = np.sqrt((z[0] - self.x)*(z[0] - self.x) + (z[1] - self.y)*(z[1] - self.y)) + 1.0e-6
		return sgn * (self.gamma_p * (d - self.radius) + (z[0] - self.x) / d * z[2] + (z[1] - self.y) / d * z[3])
		 
	def dh(self,z):
		sgn = 1.0
		if self.radius < 0.0:
			sgn = -1.0

		d = np.sqrt((z[0] - self.x)*(z[0] - self.x) + (z[1] - self.y)*(z[1] - self.y)) + 1.0e-6
		d_2 = d*d
		d_3 = d*d*d
		y_pos_m_z2 = (self.y - z[1])
		x_pos_m_z1 = (self.x - z[0])
		return sgn * np.array((
			(z[2]*(d_2 - x_pos_m_z1*x_pos_m_z1) - self.gamma_p * d_2 *x_pos_m_z1 - z[3]*x_pos_m_z1*y_pos_m_z2)/d_3,
			(z[3]*(d_2 - y_pos_m_z2*y_pos_m_z2) - self.gamma_p * d_2 *y_pos_m_z2 - z[2]*x_pos_m_z1*y_pos_m_z2)/d_3,
			-x_pos_m_z1/d,
			-y_pos_m_z2/d), dtype=np.float32)

	def d2h(self,z):
		sgn = 1.0
		if self.radius < 0.0:
			sgn = -1.0

		d = np.sqrt((z[0,:] - self.x)**2 + (z[1,:] - self.y)**2) + 1.0e-6
		z1 = z[0,:]
		z2 = z[1,:]
		z3 = z[2,:]
		z4 = z[3,:]
		z5 = z[4,:]
		gamma_p = self.gamma_p
		x_pos = self.x
		y_pos = self.y
		y_pos_m_z2 = (y_pos - z2)
		x_pos_m_z1 = (x_pos - z1)
		d2 = d**2
		d3 = d**3
		d5 = d**5
		z1_2 = z1**2
		z2_2 = z2**2
		x_pos_2 = x_pos**2
		y_pos_2 = y_pos**2
		a11 = (y_pos_m_z2*(gamma_p *(x_pos_2*y_pos - x_pos_2*z2 + 3*y_pos*z2_2 - 2*x_pos*y_pos*z1 + 2*x_pos*z1*z2 + y_pos**3 - 3*y_pos_2*z2 + y_pos*z1_2 - z1_2*z2 - z2**3) - 2*z4*x_pos_2 + 3*z3*x_pos*y_pos + 4*z4*x_pos*z1 - 3*z3*x_pos*z2 + z4*y_pos_2 - 3*z3*y_pos*z1 - 2*z4*y_pos*z2 - 2*z4*z1_2 + 3*z3*z1*z2 + z4*z2_2))/(d5)
		a12 = x_pos_m_z1 * y_pos_m_z2 * (z4 + z3 - gamma_p - (3*z3*x_pos_m_z1)/(d2) - (3*z4*y_pos_m_z2)/(d2))/(d3)
		a13 = y_pos_m_z2**2/(d3)
		a14 = -(x_pos_m_z1*y_pos_m_z2)/(d3)
		a22 = (x_pos_m_z1*(gamma_p*(x_pos**3 - 3*x_pos_2*z1 + x_pos*y_pos_2 - 2*x_pos*y_pos*z2 + 3*x_pos*z1_2 + x_pos*z2_2 - y_pos_2*z1  + 2*y_pos*z1*z2 - z1**3 - z1*z2_2)+ z3*x_pos_2 + 3*z4*x_pos*y_pos - 2*z3*x_pos*z1 - 3*z4*x_pos*z2 - 2*z3*y_pos_2 - 3*z4*y_pos*z1 + 4*z3*y_pos*z2 + z3*z1_2 + 3*z4*z1*z2 - 2*z3*z2_2))/(d5)
		a23 = -(y_pos_m_z2*x_pos_m_z1)/(d3)
		a24 = x_pos_m_z1**2/(d3)

		return sgn * np.block([
			[ a11, a12, a13, a14],
			[ a12, a22, a23, a24],
			[ a13, a23, 0, 0],
			[ a14, a24, 0, 0]], dtype=np.float32)

class LyapunovAckermann(Lyapunov):
	def __init__(self,dim=4,w1=1.0,w2=1.0,w3=1.0,epsilon=1e-6):
		self.w1=w1
		self.w2=w2
		self.w3=w3
		self.epsilon=epsilon
		Lyapunov.__init__(self,dim)

	def convert_to_polar(self,x,x_d):
		if x[3,:] < 0:
			e = -np.sqrt((x_d[0,:]-x[0,:])**2 + (x_d[1,:]-x[1,:])**2)
			phi_t = np.arctan2(x[1,:]-x_d[1,:],x[0,:]-x_d[0,:])
		else:
			e = np.sqrt((x_d[0,:]-x[0,:])**2 + (x_d[1,:]-x[1,:])**2)
			phi_t = np.arctan2(x_d[1,:]-x[1,:],x_d[0,:]-x[0,:])
		alpha = self.signed_angle_dist(phi_t,x[2,:])
		theta = self.signed_angle_dist(phi_t,x_d[2,:])

		# if alpha and theta are opposite signs, we have two confliciting headings.
		# Flip the sign of one to turn in one direction.
		# print('before:',alpha,theta)
		if np.sign(alpha) != np.sign(theta):
			if np.abs(alpha) >= np.abs(theta):
				alpha = alpha - np.sign(alpha) * 2 * np.pi
			else:
				theta = theta - np.sign(theta) * 2 * np.pi
		# print('after:',alpha,theta)
		# print(e,x[2,:]/(np.pi)*180,phi_t/(np.pi)*180,x_d[2,:]/(np.pi)*180,alpha/(np.pi)*180,theta/(np.pi)*180)

		return e,phi_t,alpha,theta

	def V(self,x,x_d):
		e,phi_t,alpha,theta = self.convert_to_polar(x,x_d)
		return 0.5 * (self.w1*alpha**2 + self.w2*theta**2 + self.w3*(x[3,:]-x_d[3,:])**2)

	def dV(self,x,x_d):
		e,phi_t,alpha,theta = self.convert_to_polar(x,x_d)
		c_tmp = (self.w1*alpha + self.w2*theta) / (e + self.epsilon)
		return np.stack((c_tmp*np.sin(phi_t),-c_tmp*np.cos(phi_t),-self.w1*alpha,self.w3*(x[3,:]-x_d[3,:])),axis=0)

class LyapunovAckermannZ(Lyapunov):
	def __init__(self,dim=4,w1=1.0, w2=1.0, w3=1.0, epsilon = 1e-6):
		self.w1=w1
		self.w2=w2
		self.w3=w3
		self.epsilon = epsilon
		Lyapunov.__init__(self,dim)


	def convert_to_polar(self,z,z_d):
		# if z[3,:] < 0:
		# 	e = -np.sqrt((z_d[0,:]-z[0,:])**2 + (z_d[1,:]-z[1,:])**2)
		# 	phi_t = np.arctan2(z[1,:]-z_d[1,:],z[0,:]-z_d[0,:])
		# else:
		e = np.sqrt((z_d[0,:]-z[0,:])**2 + (z_d[1,:]-z[1,:])**2)
		phi_t = np.arctan2(z_d[1,:]-z[1,:],z_d[0,:]-z[0,:])

		phi = np.arctan2(z[3,:],z[2,:])
		phi_d = np.arctan2(z_d[3,:],z_d[2,:])

		alpha = self.signed_angle_dist(phi_t,phi)
		theta = self.signed_angle_dist(phi_t,phi_d)
		# print('angles: ',phi,phi_t,phi_d,alpha,theta)
		# if alpha and theta are opposite signs, we have two confliciting headings.
		# Flip the sign of one to turn in one direction.
		# print('before:',alpha,theta)
		# if np.sign(alpha) != np.sign(theta):
		# 	if np.abs(alpha) >= np.abs(theta):
		# 		alpha = alpha + np.sign(alpha) * 2 * np.pi
		# 	else:
		# 		theta = theta + np.sign(theta) * 2 * np.pi
		# print('after:',alpha,theta)
		# print(e,z[2,:]/(np.pi)*180,phi_t/(np.pi)*180,z_d[2,:]/(np.pi)*180,alpha/(np.pi)*180,theta/(np.pi)*180)

		return e,phi_t,alpha,theta

	def V(self,z,z_d):
		e,phi_t,alpha,theta = self.convert_to_polar(z,z_d)
		v = np.sqrt(z[2,:]**2 + z[3,:]**2)
		v_d = np.sqrt(z_d[2,:]**2 + z_d[3,:]**2)
		return 0.5 * (self.w1*alpha**2 + self.w2*theta**2 + self.w3*(v - v_d)**2)

	def dV(self,z,z_d):
		e,phi_t,alpha,theta = self.convert_to_polar(z,z_d)
		c_tmp = (self.w1*alpha + self.w2*theta) / (e + self.epsilon)
		v = np.sqrt(z[2,:]**2 + z[3,:]**2)
		v_d = np.sqrt(z_d[2,:]**2 + z_d[3,:]**2)

		return np.stack((c_tmp*np.sin(phi_t),
						 -c_tmp*np.cos(phi_t),
						 self.w1*alpha*z[3,:]/(v**2) + self.w3*(v-v_d)*z[2,:]/v,
						 -self.w1*alpha*z[2,:]/(v**2) + self.w3*(v-v_d)*z[3,:]/v),axis=0)

class Dynamics():
	def __init__(self,xdim,udim):
		self.xdim=xdim
		self.udim=udim

	def f(self,x):
		# overload
		return None

	def g(self,x):
		# overload
		return None

	def pseudoinv(self,x):
		return np.matmul(np.linalg.inv(np.matmul(x.T,x)),x.T)

	def convert_z_to_x(self,z):
		v=np.sqrt(z[2,:]**2 + z[3,:]**2) * z[4,:]
		return np.stack((z[0,:],z[1,:],np.arctan2(z[4,:]*z[3,:],z[4,:]*z[2,:]),v),axis=0)

	def convert_x_to_z(self,x):
		v_sign = 1.0
		if x[3,:] < 0:
			v_sign = -1
		return np.stack((x[0,:],x[1,:],x[3,:]*np.cos(x[2,:]),x[3,:]*np.sin(x[2,:]),np.array([v_sign])),axis=0)

class DynamicsAckermann(Dynamics):
	def __init__(self,xdim=4,udim=2):
		Dynamics.__init__(self,xdim=xdim,udim=udim)

	def f(self,x):
		return np.stack((x[3,:] * np.cos(x[2,:]), x[3,:] * np.sin(x[2,:]), np.array([0]), np.array([0])))

	def g(self,x):
		return np.stack((np.array([0,0]),np.array([0,0]),np.append(x[3,:],0),np.array([0,1])))

class DynamicsAckermannModified(Dynamics):
	def __init__(self,xdim=4,udim=2):
		Dynamics.__init__(self,xdim=xdim,udim=udim)

	def f(self,x):
		v_modified = x[3,:]**2 / (2*np.tanh(x[3,:]))
		return np.stack((v_modified * np.cos(x[2,:]), v_modified * np.sin(x[2,:]), x[3,:], -0.5*x[3,:]))

	def g(self,x):
		return np.stack((np.array([0,0]),np.array([0,0]),np.append(x[3,:],0),np.array([0,1])))*1.2


class DynamicsAckermannZ(Dynamics):
	def __init__(self,xdim=4,udim=2,epsilon=1e-6):
		self.epsilon = epsilon
		Dynamics.__init__(self,xdim=xdim,udim=udim)

	def f(self,z):
		return np.stack((np.array([0]),np.array([0]))) * z[4,:]

	def g(self,z):
		v=(np.sqrt(z[2,:]**2 + z[3,:]**2) + self.epsilon) * z[4,:]
		return np.stack((np.concatenate((-z[3,:]*v,z[2,:]/v)),np.concatenate((z[2,:]*v,z[3,:]/v))))


class DynamicsAckermannZModified(Dynamics):
	def __init__(self,xdim=4,udim=2,epsilon=1e-6,disturbance_scale_pos = 0.5,disturbance_scale_vel = 2.0,control_input_scale = 2.0):
		self.epsilon = epsilon
		Dynamics.__init__(self,xdim=xdim,udim=udim)

		self.disturbance_scale_pos = disturbance_scale_pos
		self.disturbance_scale_vel = disturbance_scale_vel
		self.control_input_scale = control_input_scale
		
	def f(self,z):
		v=np.sqrt(z[2,:]**2 + z[3,:]**2) * z[4,:]
		theta = np.arctan2(z[3]*z[4],z[2]*z[4])
		v_disturbance_body = [np.tanh(v**2)*self.disturbance_scale_vel, (0.1+v)*self.disturbance_scale_vel]
		v_disturbance_world = [v_disturbance_body[0] * np.cos(theta) - v_disturbance_body[1] * np.sin(theta),
							   v_disturbance_body[0] * np.sin(theta) + v_disturbance_body[1] * np.cos(theta)]
		return np.array([v_disturbance_world[0], v_disturbance_world[1]])
		
	def g(self,z):
		v=(np.sqrt(z[2,:]**2 + z[3,:]**2)  + self.epsilon) * z[4,:]
		return self.control_input_scale * np.stack((np.concatenate((-z[3,:]*v,z[2,:]/v)),np.concatenate((z[2,:]*v,z[3,:]/v))))

	def step(self,z,u,dt):
		v=np.sqrt(z[2,:]**2 + z[3,:]**2)
		znext = np.zeros((self.xdim+1,1))
		znext[0:2,:] = z[0:2,:] + dt*(z[2:-1,:] + np.array([np.sin(v**3)*self.disturbance_scale_pos, np.cos(-v)*self.disturbance_scale_pos]))
		znext[2:-1,:] = z[2:-1,:] + dt*(self.f(z) + np.matmul(self.g(z),u))
		znext[-1] = z[-1]
		return znext


class QPSolve():
    def __init__(self,dyn,clf,cbf_list,u_lim,u_cost=0.0,u_prev_cost=1.0,p1_cost=1.0e10,p2_cost=1.0e10,verbose=True):
        self.xdim = dyn.xdim
        self.udim = dyn.udim
        self.dyn = dyn
        self.clf = clf
        self.cbf_list = cbf_list
        self.u_lim = u_lim
        self.p1_cost = p1_cost
        self.p2_cost = p2_cost
        self.verbose = verbose
        self.u_cost = u_cost
        self.u_prev_cost = u_prev_cost
        self.K = 0.0
        self.ksig = 1.0
        self.max_var = 1.0
        self.mu_qp_prev = np.zeros((int(self.xdim/2),1), dtype=np.float32)
        self.P = np.eye(self.xdim, dtype=np.float32)
        self.A = np.zeros((self.xdim,self.xdim), dtype=np.float32)
        self.A0 = np.block([[np.zeros((int(self.xdim/2),int(self.xdim/2))),np.eye(int(self.xdim/2))],[np.zeros((int(self.xdim/2),int(self.xdim/2))),np.zeros((int(self.xdim/2),int(self.xdim/2)))]]).astype(np.float32)
        self.G = np.block([[np.zeros((int(self.xdim/2),int(self.xdim/2)))],[np.eye(int(self.xdim/2))]]).astype(np.float32)
        self.res = None
        self.max_error = 1.0

    def update_ricatti(self,A):
        self.A = A
        Q = np.eye(self.xdim, dtype=np.float32)
        self.P = sl.solve_continuous_are(self.A,np.zeros((self.xdim,self.xdim), dtype=np.float32),Q,np.eye(self.xdim, dtype=np.float32))

    def solve(self,x,x_d,mu_d,sigDelta):
        sigDelta = sigDelta * self.ksig
        sigDelta = np.clip(sigDelta,0.0,self.max_var)
        # sigDelta = np.ones((self.xdim/2,1)) * self.max_var # for testing

        # build Q and p matrices to specify minimization expression
        Q = np.diag(np.append(np.append(np.ones(int(self.xdim/2), dtype=np.float32)*(self.u_cost + self.u_prev_cost),self.p1_cost),self.p2_cost))
        self.Q = sparse.csc_matrix(Q)
        self.p = 2*np.append(np.append(-self.mu_qp_prev*self.u_prev_cost,0),0)

        #error dynamics for clf
        e = x[:-1,:]-x_d[:-1,:]
        e = np.clip(e,-self.max_error,self.max_error)
        eTPG = np.matmul(e.T,np.matmul(self.P,self.G))
        G_dyn = np.expand_dims(np.append(np.append(eTPG,1),0),0) #append 1 for clf < d
        Gsig = np.matmul(self.G,sigDelta)
        GssG = np.matmul(Gsig,Gsig.T)
        self.trGssGP = np.trace(np.matmul(GssG,self.P))
        h_dyn = -1 * ( -0.5*np.matmul(e.T,np.matmul(Q,e))
                    + 0.5*np.matmul(e.T,np.matmul(self.P,e)) / self.clf.epsilon
                    + 0.5*self.trGssGP)

        # build constraints for barriers
        N_cbf = len(self.cbf_list)
        G_cbf = np.zeros((N_cbf,int(self.xdim/2)+2), dtype=np.float32)
        h_cbf = np.zeros((N_cbf,1), dtype=np.float32)
        A0x_Gmud = np.matmul(self.A0,x[:-1,:]) + np.matmul(self.G,mu_d)
        GssG_22 = GssG[2:,2:]
        for i, cbf in enumerate(self.cbf_list):
            h_x, dB, d2B = cbf.get_B_derivatives(x)
            G_cbf[i,:] = np.append(np.append(np.einsum('ij,jk',dB,self.G),0),1)
            trGssGd2B = np.einsum('ii',np.einsum('ij,jk',GssG_22,d2B))
            h_cbf[i,:] = -1 * (np.einsum('ij,jk',dB,A0x_Gmud)
                                - cbf.gamma * h_x
                                + 0.5*trGssGd2B)

        # build constraints for control limits
        ginv = np.linalg.inv(self.dyn.g(x))
        l_ctrl = np.matmul(ginv, mu_d - self.dyn.f(x))
        A_ctrl = ginv

        G_ctrl = np.zeros((self.udim*2,int(self.xdim/2)+2), dtype=np.float32)
        h_ctrl = np.zeros((self.udim*2,1), dtype=np.float32)
        for i in range(self.udim):
            G_ctrl[i*2,:int(self.xdim/2)] = - A_ctrl[i,:]
            h_ctrl[i*2] = - self.u_lim[i,0] + l_ctrl[i]
            G_ctrl[i*2+1,:int(self.xdim/2)] = A_ctrl[i,:]
            h_ctrl[i*2+1] = self.u_lim[i,1] - l_ctrl[i]

        # stack into one matrix and vector
        G = np.concatenate((G_dyn,G_cbf,G_ctrl),axis=0)
        h = np.concatenate((h_dyn,h_cbf,h_ctrl),axis=0)

        self.G_csc = sparse.csc_matrix(G)
        self.h = h

        # dummy lower bound
        l = np.ones(h.shape, dtype=np.float32)*np.inf * -1

        #Solve QP
        self.prob = osqp.OSQP()
        exception_called = False
        mu_bar = np.zeros((self.xdim+1), dtype=np.float32)
        # try:
        self.prob.setup(P=self.Q, q=self.p, A=self.G_csc, l=l, u=self.h, verbose=self.verbose)
        self.res = self.prob.solve()
        # except:
            # exception_called = True
        # else:
        mu_bar = self.res.x
        # if exception_called or u_bar[0] is None or np.isnan(u_bar).any():
        if mu_bar[0] is None or np.isnan(mu_bar).any():
            mu_bar = np.zeros((self.xdim+1))
            self.res = None
            print("QP failed!")
        

        self.mu_qp_prev = np.expand_dims(mu_bar[0:int(self.xdim/2)],axis=0).T

        self.V =self.clf.V(x,x_d)

        if self.verbose:
        # if True:
            print('z_ref: ', x.T)
            print('z_des: ', x_d.T)
            print('u_lim', self.u_lim)
            print('V: ', self.V)
            print('Q:', Q)
            print('p:', np.array(self.p))
            print('G_dyn:', G_dyn)
            print('h_dyn:', h_dyn)
            print('trGssGP',self.trGssGP)
            if h_cbf.shape[0] < 10:
                print('G_cbf:', G_cbf)
                print('h_cbf:', h_cbf)
            print('G_ctrl:', G_ctrl)
            print('h_ctrl:', h_ctrl)
            print('result:', mu_bar)
            print('G_all:',self.G_csc)
            print('h_all:',self.h)

        return self.mu_qp_prev

#---------------------------------------------------------------------
class ScaledGP:
	def __init__(self,xdim=1,ydim=1):
		self.xdim=xdim
		self.ydim=ydim
		self.ystd = np.ones(ydim)
		self.ymean = np.zeros(ydim)
		self.xstd = np.ones(xdim)
		self.xmean = np.zeros(xdim)
		self.m = GPy.models.GPRegression(np.zeros((1,xdim)),np.zeros((1,ydim)))

	def optimize(self,x,y,update_scaling=True, num_inducing=50):
		assert(x.shape[1] == self.xdim and y.shape[1] == self.ydim)
		assert(x.shape[0] > 0 and y.shape[0] > 0)
		
		xmean = self.xmean
		xstd = self.xstd
		ymean = self.ymean
		ystd = self.ystd
		if update_scaling:
			xmean,xstd = self.update_xscale(x)
			ymean,ystd = self.update_yscale(y)

		x = self.scalex(x,xmean,xstd)
		y = self.scaley(y,ymean,ystd)
		updated_model = GPy.models.GPRegression(x,y)
		# self.m = GPy.models.SparseGPRegression(x,y,num_inducing=num_inducing)
		updated_model.optimize('bfgs')
		self.m = updated_model

		self.xmean = xmean
		self.xstd = xstd
		self.ymean = ymean
		self.ystd = ystd

	def predict(self,x):
		x = self.scalex(x,self.xmean,self.xstd)
		mean,var = self.m.predict_noiseless(x)
		mean = self.unscaley(mean,self.ymean,self.ystd)
		var = var * self.ystd
		if mean.size == 1:
			mean = mean[0,0]
			var = var[0,0]
		return mean,var

	def update_xscale(self,x):
		xmean = np.mean(x,axis=0)
		xstd = np.std(x,axis=0)

		return xmean,xstd

	def update_yscale(self,y):
		ymean = np.mean(y,axis=0)
		ystd = np.std(y,axis=0)

		return ymean, ystd

	def scalex(self,x,xmean,xstd):
		if (xstd == 0).any():
			return (x-xmean)
		else:
			return (x - xmean) / xstd

	def scaley(self,y,ymean,ystd):
		if (ystd == 0).any():
			return (y-ymean)
		else:
			return (y - ymean) / ystd
		
	def unscalex(self,x,xmean,xstd):
		if (xstd == 0).any():
			return x + xmean
		else:
			return x * xstd + xmean

	def unscaley(self,y,ymean,ystd):
		if (ystd == 0).any():
			return y + ymean
		else:
			return y * ystd + ymean

class ModelService(object):
	_train_result = controller_adaptiveclbf.msg.TrainModelResult()

	def __init__(self,xdim,odim,use_obs, use_service = True):

		self.xdim=xdim
		self.odim=odim
		self.use_obs = use_obs
		self.verbose = True

	def predict(self,req):
		# overload
		return None

	def train(self,goal):
		# overload
		return None

	def add_data(self,req):
		# overload
		return None

class ModelGPService(ModelService):
	def __init__(self,xdim,odim,use_obs=False,use_service=True):
		ModelService.__init__(self,xdim,odim,use_obs,use_service)
		# note:  use use_obs and observations with caution.  model may overfit to this input.
		model_xdim=self.xdim/2
		if self.use_obs:
			 model_xdim += self.odim
		model_ydim=self.xdim/2

		self.m = ScaledGP(xdim=model_xdim,ydim=model_ydim)
		self.y = np.zeros((0,model_ydim))
		self.Z = np.zeros((0,model_xdim))
		self.N_data = 400

	def rotate(self,x,theta):
		x_body = np.zeros((2,1))
		x_body[0] = x[0] * np.cos(theta) + x[1] * np.sin(theta)
		x_body[1] = -x[0] * np.sin(theta) + x[1] * np.cos(theta)
		return x_body

	def make_input(self,x,obs):
		# format input vector
		theta = obs[0]
		x_body = self.rotate(x[2:-1,:],theta)
		if self.use_obs:
			Z = np.concatenate((x_body,obs[1:,:])).T
		else:
			Z = np.concatenate((x_body)).T

		return Z

	def predict(self,x,obs):
		# format the input and use the model to make a prediction.
		Z = self.make_input(x,obs)
		y, var = self.m.predict(Z)
		# theta = np.arctan2(x[3]*x[4],x[2]*x[4])
		theta=obs[0]
		y_out = self.rotate(y.T,-theta)

		return y_out,var

	def train(self):
		self.m.optimize(self.Z,self.y)
				
	def add_data(self,x_next,x,mu_model,obs,dt):
		# add a sample to the history of data
		x_dot = (x_next[2:-1,:]-x[2:-1,:])/dt
		ynew = x_dot - mu_model
		Znew = self.make_input(x,obs)
		# theta = np.arctan2(x[3]*x[4],x[2]*x[4])
		theta=obs[0]
		ynew_rotated = self.rotate(ynew,theta)
		self.y = np.concatenate((self.y,ynew_rotated.T))
		self.Z = np.concatenate((self.Z,Znew))

		# throw away old samples if too many samples collected.
		if self.y.shape[0] > self.N_data:
			self.y = self.y[-self.N_data:,:]
			self.Z = self.Z[-self.N_data:,:]
#---------------------------------------------------------------------
class AdaptiveClbf(object):
	def __init__(self,odim=2, use_service = True, ):
		self.xdim = 4
		self.udim = 2
		self.odim = odim
		self.use_service = use_service

		self.u_lim = np.array([[-2.0,2.0],[-1.0,1.0]])
		self.K=np.block([[np.zeros((2,2)), np.zeros((2,2))],[np.eye(2), np.eye(2)]])

		self.dyn = DynamicsAckermannZ()

		self.model_trained = False
		
		if self.use_service:
			## train model action
			self.train_model_action_client = actionlib.SimpleActionClient('train_model_service', controller_adaptiveclbf.msg.TrainModelAction)
			self.train_model_action_client.wait_for_server()
			self.train_model_goal = controller_adaptiveclbf.msg.TrainModelGoal()

			## add data srv
			rospy.wait_for_service('add_data_2_model')
			self.model_add_data_srv = rospy.ServiceProxy('add_data_2_model', AddData2Model)

			## predict srv
			rospy.wait_for_service('predict_model')
			self.model_predict_srv = rospy.ServiceProxy('predict_model', PredictModel)
		else:
			# setup non-service model object
			# from model_service import ModelVanillaService
			# self.model = ModelVanillaService(self.xdim,self.odim,use_obs=True,use_service=False)

			# from model_service import ModelGPService
			self.model = ModelGPService(self.xdim,self.odim,use_obs=True,use_service=False)

			# from model_service import ModelALPaCAService
			# self.model = ModelALPaCAService(self.xdim,self.odim,use_obs=True,use_service=False)

		self.clf = LyapunovAckermannZ(w1=10.0,w2=1.0,w3=1.0,epsilon=1.0)
		self.qpsolve = QPSolve(dyn=self.dyn,cbf_list=[],clf=self.clf,u_lim=self.u_lim,u_cost=0.0,u_prev_cost=1.0,p1_cost=1.0e8,p2_cost=1.0e8,verbose=False)

		x_init = np.zeros((self.xdim,1))
		x_init[3] = 0.01
		self.z_ref = self.dyn.convert_x_to_z(x_init)
		self.z = copy.copy(self.z_ref)
		self.z_ref_dot = np.zeros((self.xdim,1))
		self.z_dot = np.zeros((int(self.xdim/2),1))
		self.z_prev = copy.copy(self.z)
		self.y_out = np.zeros((int(self.xdim/2),1))
		self.mu_prev = np.zeros((self.xdim,1))
		self.u_prev = np.zeros((self.udim,1))
		self.u_prev_prev = np.zeros((self.udim,1))
		self.obs_prev = np.zeros((self.odim,1))
		self.dt = 0.1
		self.max_error = 1.0

		self.barrier_locations={}
		self.barrier_locations["x"]=np.array([])
		self.barrier_locations["y"]=np.array([])
		self.barrier_radius = 1.0

		self.measurement_noise = 1.0
		self.true_dyn = None
		self.true_predict_error = 0.0
		self.predict_error = 0
		self.predict_var = np.zeros((int(self.xdim/2),1))

		self.debug={}

	def wrap_angle(self,a):
		return (a + np.pi) % (2 * np.pi) - np.pi

	def saturate(self,u,ulim):
		return np.array([np.clip(u[0],ulim[0][0],ulim[0][1]), np.clip(u[1],ulim[1][0],ulim[1][1])])

	def update_params(self,params):
		self.params = params
		print("updated params!")
		self.vehicle_length = self.params["vehicle_length"]
		self.steering_limit = self.params["steering_limit"]
		self.max_accel = self.params["max_accel"]
		self.min_accel = self.params["min_accel"]
		self.u_lim = np.array([[np.tan(-self.steering_limit)/self.vehicle_length,np.tan(self.steering_limit)/self.vehicle_length],
								[self.min_accel,self.max_accel]],dtype=np.float32)
		self.qpsolve.u_lim = self.u_lim

		self.k1 = self.params["kp_z"]
		self.k2 = self.params["kd_z"]
		self.A=np.block([[np.zeros((2,2),dtype=np.float32), np.eye(2,dtype=np.float32)],[-self.k1*np.eye(2,dtype=np.float32), -self.k2*np.eye(2,dtype=np.float32)]])
		self.qpsolve.update_ricatti(self.A)
		self.K=np.block([[self.k1*np.eye(2,dtype=np.float32), self.k2*np.eye(2,dtype=np.float32)]])
		self.max_error = self.params["max_error"]

		self.clf.epsilon = self.params["clf_epsilon"]
		self.measurement_noise = self.params["measurement_noise"]

		self.qpsolve.u_cost = self.params["qp_u_cost"]
		self.qpsolve.u_prev_cost = self.params["qp_u_prev_cost"]
		self.qpsolve.p1_cost = self.params["qp_p1_cost"]
		self.qpsolve.p2_cost = self.params["qp_p2_cost"]
		self.qpsolve.verbose = self.params["qp_verbose"]
		self.qpsolve.ksig = self.params["qp_ksig"]
		self.qpsolve.max_var = self.params["qp_max_var"]
		self.verbose = self.params["verbose"]

		self.dt = self.params["dt"]

	def update_barrier_locations(self,x,y,radius):
		self.barrier_locations["x"] = x
		self.barrier_locations["y"] = y
		self.barrier_radius = radius

	def update_barriers(self):
		cbf_list = []
		bar_loc_x = copy.copy(self.barrier_locations["x"])
		bar_loc_y = copy.copy(self.barrier_locations["y"])
		bar_rad = self.barrier_radius

		if self.params["use_barrier_vel"]:
			cbf_list = cbf_list + \
							[BarrierAckermannVelocityZ(bound_from_above=True, v_lim = self.params["max_velocity"], gamma=self.params["barrier_vel_gamma"]),
							BarrierAckermannVelocityZ(bound_from_above=False, v_lim = self.params["min_velocity"], gamma=self.params["barrier_vel_gamma"])]

		if self.params["use_barrier_pointcloud"]:
			cbf_list = cbf_list + \
							[BarrierAckermannPointZ(x=bar_loc_x[i],y=bar_loc_y[i], radius=bar_rad, gamma_p=self.params["barrier_pc_gamma_p"], gamma=self.params["barrier_pc_gamma"]) for i in range(bar_loc_x.size)]
		self.qpsolve.cbf_list = cbf_list

	def get_control(self,z,z_ref,z_ref_dot,dt,obs=None,train=False,use_model=False,add_data=True,check_model=True,use_qp=True):
		assert z.shape[0] == self.xdim + 1
		assert z_ref_dot.shape[0] == self.xdim + 1
		assert z_ref.shape[0] == self.xdim + 1

		self.update_barriers()

		self.z = copy.copy(z.astype(np.float32))
		self.z_ref = copy.copy(z_ref.astype(np.float32))
		self.obs = np.array(obs,dtype=np.float32)

		mu_ad = np.zeros((int(self.xdim/2),1),dtype=np.float32)
		mDelta = np.zeros((int(self.xdim/2),1),dtype=np.float32)
		sigDelta = np.zeros((int(self.xdim/2),1),dtype=np.float32)
		rho = np.zeros((int(self.xdim/2),1),dtype=np.float32)
		trueDelta = np.zeros((int(self.xdim/2),1),dtype=np.float32)

		e = self.z_ref[:-1,:]-self.z[:-1,:]
		mu_pd = np.matmul(self.K,e)
		mu_pd = np.clip(mu_pd,-self.max_error,self.max_error)

		# self.z_ref_dot = (self.z_ref_next[:-1,:] - self.z_ref[:-1,:]) / dt # use / self.dt for test_adaptive_clbf
		self.z_ref_dot = copy.copy(z_ref_dot)
		mu_rm = self.z_ref_dot[2:-1]

		mu_model = np.matmul(self.dyn.g(self.z_prev),self.u_prev) + self.dyn.f(self.z_prev)
		mu_model = np.clip(mu_model,-self.max_error,self.max_error)

		if add_data:
			if self.use_service:
				try:
					self.model_add_data_srv(self.z.flatten(),self.z_prev.flatten(),mu_model.flatten(),self.obs_prev.flatten(),dt)
				except:
					print("add data service unavailable")
			else:
				req = AddData2Model()
				req.x_next = self.z.flatten()
				req.x = self.z_prev.flatten()
				req.mu_model = mu_model.flatten()
				req.obs = self.obs_prev.flatten()
				req.dt = dt
				self.model.add_data(req)

		self.z_dot = (self.z[2:-1,:]-self.z_prev[2:-1,:])/dt - mu_model

		# if check_model and self.model.model_trained:
		if check_model and self.model_trained:
			# check how the model is doing.  compare the model's prediction with the actual sampled data.
			predict_service_success = False
			result = None
			# if self.use_service:
				# try:
					# result = self.model_predict_srv(self.z_prev.flatten(),self.obs_prev.flatten())
					# if result.result:
						# predict_service_success = True
				# except:
					# print("predict service unavailable")
			# else:
			self.y_out,var = self.model.predict(z_prev,obs_prev)
			predict_service_success = True

			if predict_service_success:

				if self.verbose:
					print("predicted y_out: ", self.y_out)
					print("predicted ynew: ", self.z_dot)
					print("predicted var: ", var)

				self.predict_error = np.linalg.norm(self.y_out - self.z_dot)
				self.predict_var = var

		if use_model and self.model_trained:
			predict_service_success = False
			result = None
			if self.use_service:
				try:
					result = self.model_predict_srv(self.z.flatten(),self.obs.flatten())
					if result.result:
						predict_service_success = True
				except:
					print("predict service unavailable")
			else:
				mDelta,sigDelta = self.model.predict(z_prev,obs_prev)
				predict_service_success = True

			if predict_service_success:

				# log error if true system model is available
				if self.true_dyn is not None:
					trueDelta = self.true_dyn.f(self.z) - self.dyn.f(self.z)
					self.true_predict_error = np.linalg.norm(trueDelta - mDelta)

				# rho = self.measurement_noise / (self.measurement_noise + (sigDelta - 1.0) + 1e-6)
				rho[0] = self.measurement_noise / (self.measurement_noise + np.linalg.norm(sigDelta) + 1e-6)
				#rho = np.linalg.norm(sigDelta) * self.measurement_noise
				mu_ad = mDelta * rho[0]
				sigDelta = sigDelta * (1.0 - rho[0]) 
		else:
			sigDelta = np.ones((int(self.xdim/2),1))

		mu_d = mu_rm + mu_pd - mu_ad

		self.mu_qp = np.zeros((int(self.xdim/2),1))
		if use_qp:
			self.mu_qp = self.qpsolve.solve(self.z,self.z_ref,mu_d,sigDelta)

		self.mu_new = mu_d + self.mu_qp
		self.u_new = np.matmul(np.linalg.inv(self.dyn.g(self.z)), (self.mu_new-self.dyn.f(self.z)))
		# print(np.linalg.inv(self.dyn.g(self.z)))
		
		u_new_unsaturated = copy.copy(self.u_new)
		#saturate in case constraints are not met
		self.u_new = self.saturate(self.u_new,self.u_lim)

		if self.verbose:
			print('z: ', self.z.T)
			print('z_ref: ', self.z_ref.T)
			print('mu_rm', mu_rm)
			print('mu_pd', mu_pd)
			print('mu_ad', mu_ad)
			print('mu_d', mu_d)
			print('mu_model', mu_model)
			print('rho', rho)
			print('mu_qp', self.mu_qp)
			print('mu',self.mu_new)
			print('u_new', self.u_new)
			print('u_unsat', u_new_unsaturated)
			print('trueDelta',trueDelta)
			print('true predict error', self.true_predict_error)
			print('mDelta', mDelta)
			print('sigDelta', sigDelta)

		self.debug["z"] = self.z.flatten().tolist()
		self.debug["z_ref"] = self.z_ref.flatten().tolist()
		self.debug["z_dot"] = self.z_dot.flatten().tolist()
		self.debug["y_out"] = self.y_out.flatten().tolist()
		self.debug["mu_rm"] = self.z_ref_dot.flatten().tolist()
		self.debug["mu_pd"] = mu_pd.flatten().tolist()
		self.debug["mu_ad"] = mu_ad.flatten().tolist()
		self.debug["mu_model"] = mu_model.flatten().tolist()
		self.debug["rho"] = rho.flatten().tolist()
		self.debug["mu_qp"] = self.mu_qp.flatten().tolist()
		self.debug["mu"] = self.mu_new.flatten().tolist()
		self.debug["u_new"] = self.u_new.flatten().tolist()
		self.debug["u_unsat"] = u_new_unsaturated.flatten().tolist()
		self.debug["trueDelta"] = trueDelta.flatten().tolist()
		self.debug["true_predict_error"] = self.true_predict_error
		self.debug["mDelta"] = mDelta.flatten().tolist()
		self.debug["sigDelta"] = sigDelta.flatten().tolist()

		self.mu_prev = copy.copy(self.mu_new)
		self.u_prev_prev = copy.copy(self.u_prev)
		self.u_prev = copy.copy(self.u_new)
		self.obs_prev = copy.copy(self.obs)
		self.z_prev = copy.copy(self.z)

		self.controls = np.zeros(self.udim)
		self.controls[0] = np.arctan(self.u_new[0] * self.vehicle_length)
		self.controls[1] = self.u_new[1]
		return self.controls

odim = 2
np.random.seed(0)
adaptive_clbf = AdaptiveClbf(odim=odim, use_service = False)
adaptive_clbf_qp = AdaptiveClbf(odim=odim, use_service = False)
adaptive_clbf_ad = AdaptiveClbf(odim=odim, use_service = False)
adaptive_clbf_pd = AdaptiveClbf(odim=odim, use_service = False)

params={}
params["vehicle_length"] = 0.25
params["steering_limit"] = 0.75
params["max_accel"] = 1.0
params["min_accel"] = -1.0
params["kp_z"] = 1.0
params["kd_z"] = 1.0
params["clf_epsilon"] = 100.0


params["qp_u_cost"] = 100.0
params["qp_u_prev_cost"] = 1.0
params["qp_p1_cost"] = 1.0
params["qp_p2_cost"] = 1.0e12
params["qp_max_var"] = 1.5
params["qp_verbose"] = False
params["max_velocity"] = 2.0
params["min_velocity"] = 0.5
params["barrier_vel_gamma"] = 10.0
params["use_barrier_vel"] = True
params["use_barrier_pointcloud"] = True
params["barrier_radius"] = 1.0
params["barrier_radius_velocity_scale"] = 0.0
params["barrier_pc_gamma_p"] = 5.0
params["barrier_pc_gamma"] = 1.0
params["verbose"] = False
params["dt"] = 0.1
params["max_error"] = 10.0

# alpaca params
# params["qp_ksig"] = 2.0e3
# params["measurement_noise"] = 1.0e-3

# gp params
# params["qp_ksig"] = 1.0e5
# params["measurement_noise"] = 1.0

# vanilla nn params
params["qp_ksig"] = 1.0e2
params["measurement_noise"] = 1.0

params["N_data"] = 600
params["learning_verbose"] = False
params["N_updates"] = 50
params["meta_batch_size"] = 50
params["data_horizon"] = 20
params["test_horizon"] = 30
params["learning_rate"] = 0.001
params["min_datapoints"] = 50
params["save_data_interval"] = 10000

true_dyn = DynamicsAckermannZModified(disturbance_scale_pos = 0.0, disturbance_scale_vel = -1.0, control_input_scale = 1.0)

adaptive_clbf.update_params(params)
adaptive_clbf_qp.update_params(params)
adaptive_clbf_ad.update_params(params)
adaptive_clbf_pd.update_params(params)
adaptive_clbf.true_dyn = true_dyn
adaptive_clbf_ad.true_dyn = true_dyn

barrier_x = np.array([5,15,25,35,45,55])
barrier_y = np.array([0,-0.5,0.5,-0.5,0.5,0])
# barrier_x = np.linspace(0,1,100)
# barrier_y = np.linspace(0,1,100)
# barrier_x = np.array([])
# barrier_y = np.array([])
adaptive_clbf.update_barrier_locations(barrier_x,barrier_y,params["barrier_radius"])
adaptive_clbf_qp.update_barrier_locations(barrier_x,barrier_y,params["barrier_radius"])
adaptive_clbf_ad.update_barrier_locations(barrier_x,barrier_y,params["barrier_radius"])
adaptive_clbf_pd.update_barrier_locations(barrier_x,barrier_y,params["barrier_radius"])

x0=np.array([[0.0],[0.0],[0.0],[0.0001]])
z0 = true_dyn.convert_x_to_z(x0)

T = 1
dt = 0.1
N = int(round(T/dt))
t = np.linspace(0,T-2*dt,N-1)
xdim=4
udim=2

train_interval = 40
start_training = 100

width = 1.0
speed = 1.0
freq = 1.0/10
x_d = np.stack((t * speed, width * np.sin(2 * np.pi * t * freq),np.zeros(N-1), np.zeros(N-1)))
x_d[2,:-1] = np.arctan2(np.diff(x_d[1,:]),np.diff(x_d[0,:]))
x_d[3,:-1] = np.sqrt(np.diff(x_d[0,:])**2 + np.diff(x_d[1,:])**2)/dt
x_d[2,-1]=x_d[2,-2]
x_d[3,-1]=x_d[3,-2]

u = np.zeros((udim,N))
x = np.zeros((xdim,N-1))
x[:,0:1] = x0

u_qp = np.zeros((udim,N))
x_qp = np.zeros((xdim,N-1))
x_qp[:,0:1] = x0

u_pd = np.zeros((udim,N))
x_pd = np.zeros((xdim,N-1))
x_pd[:,0:1] = x0

u_ad = np.zeros((udim,N))
x_ad = np.zeros((xdim,N-1))
x_ad[:,0:1] = x0

z_d = np.zeros((xdim+1,N-1))
z = np.zeros((xdim+1,N-1))
z[:,0:1] = z0
z_qp = np.zeros((xdim+1,N-1))
z_qp[:,0:1] = z0
z_pd = np.zeros((xdim+1,N-1))
z_pd[:,0:1] = z0
z_ad = np.zeros((xdim+1,N-1))
z_ad[:,0:1] = z0

z_d_dot = np.zeros((xdim+1,1))

prediction_error_ad = np.zeros(N)
prediction_error_true_ad = np.zeros(N)
prediction_error = np.zeros(N)
prediction_error_true = np.zeros(N)
prediction_var = np.zeros((int(xdim/2),N))
prediction_var_ad = np.zeros((int(xdim/2),N))
trGssGP = np.zeros(N)

i=0
z_d[:,i+1:i+2] = true_dyn.convert_x_to_z(x_d[:,i+1:i+2])

bar = Bar(max=N-1)
for i in range(N-2):
	bar.next()
	start = time.time()

	if i < N-3:
		z_d[:,i+2:i+3] = true_dyn.convert_x_to_z(x_d[:,i+2:i+3])
		z_d_dot = (z_d[:,i+2:i+3] - z_d[:,i+1:i+2])/dt

	if i == 0:
		add_data = False
	else:
		add_data = True

	# u[:,i+1] = adaptive_clbf.get_control(z[:,i:i+1],z_d[:,i+1:i+2],z_d_dot,dt=dt,obs=np.concatenate([x_ad[2,i:i+1],u_ad[:,i]]),use_model=True,add_data=add_data,use_qp=True)
	# if (i - start_training -1 ) % train_interval == 0 and i > start_training:
			# adaptive_clbf.model.train()
			# adaptive_clbf.model_trained = True
			
	u_ad[:,i+1] = adaptive_clbf_ad.get_control(z_ad[:,i:i+1],z_d[:,i+1:i+2],z_d_dot,dt=dt,obs=np.concatenate([x_ad[2,i:i+1],u_ad[:,i]]),use_model=True,add_data=add_data,use_qp=False)
	if (i - start_training - 1) % train_interval == 0 and i > start_training:
		adaptive_clbf_ad.model.train()
		adaptive_clbf_ad.model_trained = True
		
	# print(z_qp[:,i:i+1],z_d[:,i+1:i+2],z_d_dot)
	u_qp[:,i+1] = adaptive_clbf_qp.get_control(z_qp[:,i:i+1],z_d[:,i+1:i+2],z_d_dot,dt=dt,obs=[],use_model=False,add_data=False,use_qp=True)
	print(u_qp[:,i+1])
	# u_pd[:,i+1] = adaptive_clbf_pd.get_control(z_pd[:,i:i+1],z_d[:,i+1:i+2],z_d_dot,dt=dt,obs=[],use_model=False,add_data=False,use_qp=False)

	# dt = np.random.uniform(0.05,0.15)
	# c = copy.copy(u[:,i+1:i+2])
	c_ad = copy.copy(u_ad[:,i+1:i+2])
	c_qp = copy.copy(u_qp[:,i+1:i+2])
	# c_pd = copy.copy(u_pd[:,i+1:i+2])

	# c[0] = np.tan(c[0])/params["vehicle_length"]
	c_ad[0] = np.tan(c_ad[0])/params["vehicle_length"]
	c_qp[0] = np.tan(c_qp[0])/params["vehicle_length"]
	# c_pd[0] = np.tan(c_pd[0])/params["vehicle_length"]

	# z[:,i+1:i+2] = true_dyn.step(z[:,i:i+1],c,dt)
	z_ad[:,i+1:i+2] = true_dyn.step(z_ad[:,i:i+1],c_ad,dt)
	z_qp[:,i+1:i+2] = true_dyn.step(z_qp[:,i:i+1],c_qp,dt)
	# z_pd[:,i+1:i+2] = true_dyn.step(z_pd[:,i:i+1],c_pd,dt)

	# x[:,i+1:i+2] = true_dyn.convert_z_to_x(z[:,i+1:i+2])
	x_ad[:,i+1:i+2] = true_dyn.convert_z_to_x(z_ad[:,i+1:i+2])
	x_qp[:,i+1:i+2] = true_dyn.convert_z_to_x(z_qp[:,i+1:i+2])
	# x_pd[:,i+1:i+2] = true_dyn.convert_z_to_x(z_pd[:,i+1:i+2])

	print('Iteration ', i, ', Time elapsed (ms): ', (time.time() - start)*1000)


