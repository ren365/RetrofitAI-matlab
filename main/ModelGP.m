classdef  ModelGP
	% untest
	properties
		xdim
		udim
		odim
		use_obs
		m
		y
		z
		N_data

		model_trained
	end
	methods (Access = public)
		function obj = ModelGP(xdim,udim,odim,use_obs,N_data)
			% # note:  use use_obs and observations with caution.  model may overfit to this input.
			obj.xdim=xdim;%4
			obj.udim=udim;
			obj.odim=odim;%2
			obj.use_obs = use_obs;
			model_xdim=obj.xdim/2;
			if obj.use_obs
				 model_xdim = model_xdim + obj.odim;
			end
			model_ydim=obj.xdim/2;
			obj.m = ScaledGP(model_xdim,model_ydim);
			obj.y = [];
			obj.z = [];
			obj.N_data = N_data;

			obj.model_trained = false;
		
		end
		
		function result = rotate(obj,x,theta)
			result = zeros(2,1);
			result(1) = x(1) * cos(theta) + x(2) * sin(theta);
			result(2) = -x(1) * sin(theta) + x(2) * cos(theta);
		end

		function result = make_input(obj,x,obs)
			% format input vector
			theta = obs(1);
			x_body = obj.rotate(x(3:4,:),theta);
			if obj.use_obs
				z = [x_body',obs(2:end)]';
			else
				z = [x_body];
			end
			result = z;
		end
		
		
		function obj = add_data(obj,req)
			% get values
			x_next = req.x_next;  % req.x_next' ?
			x = req.x;
			mu_model = req.mu_model;
			obs = req.obs;
			dt = req.dt;
			% add a sample to the history of data
			x_dot = (x_next(3:4,:)-x(3:4,:))/dt;
			ynew = x_dot - mu_model;
			Znew = obj.make_input(x,obs);
			
			% theta = arctan2(x[3]*x[4],x[2]*x[4])
			theta=obs(1);
			ynew_rotated = obj.rotate(ynew,theta);
			obj.y = [obj.y;ynew_rotated'];
			obj.z = [obj.z;Znew'];
			% throw away old samples if too many samples collected.
			if length(obj.y) > obj.N_data
				obj.y = obj.y((end-obj.N_data):end,:);
				obj.z = obj.z((end-obj.N_data):end,:);
			end
			
		end
		
		
		function obj = train(obj)
			if length(obj.z) > 0
				obj.m = obj.m.optimize(obj.z,obj.y);
				obj.model_trained = true;
			end
		end
		
		function result = predict(obj,req)
			% input
			x = req.x;
			obs = req.obs;
			
			% # format the input and use the model to make a prediction.
			z = obj.make_input(x,obs);
			[y,var] = obj.m.predict(z);
			% # theta = arctan2(x[3]*x[4],x[2]*x[4])
			theta=obs(1);
			y_out = obj.rotate(y',-theta);
			
			% results
			result = {};
			result.y_out = y_out;
			result.var = var';
			result.result = true;

		end

	end
end

