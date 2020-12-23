import numpy as np

class DynamicsAckermannZModified():
	def __init__(self,xdim=4,udim=2,epsilon=1e-6,disturbance_scale_pos = 0.5,disturbance_scale_vel = 2.0,control_input_scale = 2.0):
		self.epsilon = epsilon

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
		
	def convert_z_to_x(self,z):
		v=np.sqrt(z[2,:]**2 + z[3,:]**2) * z[4,:]
		return np.stack((z[0,:],z[1,:],np.arctan2(z[4,:]*z[3,:],z[4,:]*z[2,:]),v),axis=0)

	def convert_x_to_z(self,x):
		v_sign = 1.0
		if x[3,:] < 0:
			v_sign = -1
		return np.stack((x[0,:],x[1,:],x[3,:]*np.cos(x[2,:]),x[3,:]*np.sin(x[2,:]),np.array([v_sign])),axis=0)

		
true_dyn = DynamicsAckermannZModified(disturbance_scale_pos = 0.0, disturbance_scale_vel = -1.0, control_input_scale = 1.0)

x0=np.array(((0.0),(0.0),(0.0),(0.0001)))
print(x0)
z0 = true_dyn.convert_x_to_z(x0)
print(z0)

T = 40
dt = 0.1
N = int(round(T/dt))
t = np.linspace(0,T-2*dt,N-1)

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

z_d = np.zeros((xdim+1,N-1))
z = np.zeros((xdim+1,N-1))
z[:,0:1] = z0

z_d_dot = np.zeros((xdim+1,1))

i=0
z_d[:,i+1:i+2] = true_dyn.convert_x_to_z(x_d[:,i+1:i+2])









