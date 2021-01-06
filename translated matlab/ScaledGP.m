classdef  ScaledGP
	% untest
	properties
		
	end
	methods (Access = public)
		function obj = ScaledGP(xdim,ydim)
			obj.xdim = xdim;
			obj.ydim = ydim;
			obj.ystd = ones(1,ydim);
			obj.ymean = zeros(1,ydim);
			obj.xstd = ones(1,xdim);
			obj.xmean = zeros(1,xdim);
			obj.m = fitrgp(zeros(1,xdim),zeros(1,ydim));
		end
		
		function obj = optimize(obj,x,y)

		end
		
		def optimize(obj,x,y,update_scaling=True, num_inducing=50):
		
		xmean = obj.xmean
		xstd = obj.xstd
		ymean = obj.ymean
		ystd = obj.ystd
		if update_scaling:
			xmean,xstd = obj.update_xscale(x)
			ymean,ystd = obj.update_yscale(y)

		x = obj.scalex(x,xmean,xstd)
		y = obj.scaley(y,ymean,ystd)
		updated_model = GPy.models.GPRegression(x,y)
		# obj.m = GPy.models.SparseGPRegression(x,y,num_inducing=num_inducing)
		updated_model.optimize('bfgs')
		obj.m = updated_model

		obj.xmean = xmean
		obj.xstd = xstd
		obj.ymean = ymean
		obj.ystd = ystd
				
		function result = V(obj,z,z_d)
		
		end

	end
end


	

	def predict(obj,x):
		x = obj.scalex(x,obj.xmean,obj.xstd)
		mean,var = obj.m.predict_noiseless(x)
		mean = obj.unscaley(mean,obj.ymean,obj.ystd)
		var = var * obj.ystd
		if mean.size == 1:
			mean = mean[0,0]
			var = var[0,0]
		return mean,var