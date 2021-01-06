import numpy as np
import GPy

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




