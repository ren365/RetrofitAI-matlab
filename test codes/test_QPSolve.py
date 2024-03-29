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
	
dyn = DynamicsAckermannZ()
clf = LyapunovAckermannZ(w1=10.0,w2=1.0,w3=1.0,epsilon=1.0)
u_lim = np.array([[-2.0,2.0],[-1.0,1.0]])
# u_lim = np.array([[np.tan(-steering_limit)/vehicle_length,np.tan(steering_limit)/vehicle_length],
								# [min_accel,max_accel]],dtype=np.float32)
# cbf list
cbf_list = []
cbf_list = cbf_list + \
			[BarrierAckermannVelocityZ(bound_from_above=True, v_lim = 2.0, gamma=10.0),
			BarrierAckermannVelocityZ(bound_from_above=False, v_lim = 0.5, gamma=10.0)]
cbf_list = cbf_list + \
			[BarrierAckermannPointZ(x=1,y=0, radius=3.0, gamma_p=5.0, gamma=1.0)]
# qpsolve
qpsolve = QPSolve(dyn=dyn,cbf_list=cbf_list,clf=clf,u_lim=u_lim,u_cost=0.0,u_prev_cost=1.0,p1_cost=1.0e8,p2_cost=1.0e8,verbose=False)
qpsolve.u_cost = 100.0
qpsolve.u_prev_cost = 1.0
qpsolve.p1_cost = 1.0
qpsolve.p2_cost = 1.0e12
qpsolve.verbose = False
qpsolve.ksig = 1.0e2
qpsolve.max_var = 1.5
qpsolve.u_lim = u_lim
qpsolve.cbf_list = cbf_list

A=np.block([[np.zeros((2,2),dtype=np.float32), np.eye(2,dtype=np.float32)],[-np.eye(2,dtype=np.float32), -np.eye(2,dtype=np.float32)]])
qpsolve.update_ricatti(A)

sigDelta = np.array([[0.35],[0.6]])
x = np.zeros((5,4))
x[0,0]=0.1
x[1,0]=0.2
x[2,0]=0.3
x[3,0]=0.4
x[4,0]=0.5
# x[:,0:1]
x[0,1]=0.4
x[1,1]=0.3
x[2,1]=0.2
x[3,1]=0.1
x[4,1]=0.5
z_ref_dot = np.zeros((int(4/2),1))
mu_qp = qpsolve.solve(x[:,0:1],x[:,1:2],z_ref_dot,sigDelta)
print(mu_qp)

