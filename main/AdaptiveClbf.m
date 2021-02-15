classdef AdaptiveClbf
	properties
		xdim
		udim
		odim
		u_lim
		K
		dyn
		model_trained

		model
		
		clf
		qpsolve
		
		z_ref
		z
		z_ref_dot
		z_dot
		z_prev
		y_out
		mu_prev
		u_prev
		u_prev_prev
		obs_prev
		dt
		max_error
		
		barrier_locations
		barrier_radius
		
		measurement_noise
		true_dyn
		true_predict_error
		predict_error
		predict_var
		params
		vehicle_length
		steering_limit
		max_accel
		min_accel
		k1
		k2
		clf_epsilon
		
		mu_qp
		mu_new
		u_new
		controls
		
		A
		obs
	end
	methods (Access = public)
	
		function obj = AdaptiveClbf(odim)
			obj.xdim = 4;
			obj.udim = 2;
			obj.odim = odim;
			obj.u_lim = [-2.0,2.0;-1.0,1.0];
			obj.K = [zeros(2), zeros(2);eye(2), eye(2)];
			obj.dyn = DynamicsAckermannZ();
			obj.model_trained = false;
			%%%% USE ROSE - UNTRASLATE
				% from model_service import ModelVanillaService
				obj.model = ModelGP(obj.xdim,obj.udim,obj.odim,true);
			%%%% USE ROSE - UNTRASLATE - END
			obj.clf = LyapunovAckermannZ(10.0,1.0,1.0,1.0);
			obj.qpsolve = QPSolve(obj.dyn,[],obj.clf,obj.u_lim,0.0,1.0,1.0e8,1.0e8,false);
			
			x_init = zeros(obj.xdim,1);
			x_init(4) = 0.01;
			obj.z_ref = obj.dyn.convert_x_to_z(x_init);
			obj.z = obj.z_ref;
			obj.z_ref_dot = zeros(obj.xdim,1);
			obj.z_dot = zeros(obj.xdim/2,1);
			obj.z_prev = obj.z;
			obj.y_out = zeros(obj.xdim/2,1);
			obj.mu_prev = zeros(obj.xdim,1);
			obj.u_prev = zeros(obj.udim,1);
			obj.u_prev_prev = zeros(obj.udim,1);
			obj.obs_prev = zeros(obj.odim+1,1);
			obj.dt = 0.1;
			obj.max_error = 1.0;
			
			obj.barrier_locations = containers.Map;
			obj.barrier_locations('x') = [];
			obj.barrier_locations('y') = [];
			obj.barrier_radius = 1.0;
			
			obj.measurement_noise = 1.0;
			obj.true_dyn = NaN;
			obj.true_predict_error = 0.0;
			obj.predict_error = 0;
			obj.predict_var = zeros(obj.xdim/2,1);
			obj.params = NaN;
			obj.vehicle_length = 0.25;
			obj.steering_limit = 0.75;
			obj.max_accel = 1;
			obj.min_accel = -1;
			obj.k1 = 1;
			obj.k2 = 1;
			obj.clf_epsilon = 100;
			
			obj.A = NaN;
			obj.obs = zeros(obj.odim+1,1);
			
			obj.mu_qp = NaN;
			obj.mu_new = NaN;
			obj.u_new = NaN;
			obj.controls = zeros(obj.udim,1);
		end
		
		% untested
		function obj = update_params(obj,params)
			obj.params = params;
			disp("updated params!")
			obj.vehicle_length = params.vehicle_length;
			obj.steering_limit = params.steering_limit;
			obj.max_accel = params.max_accel;
			obj.min_accel = params.min_accel;
			obj.u_lim = [tan(-obj.steering_limit)/obj.vehicle_length,tan(obj.steering_limit)/obj.vehicle_length;...
						obj.min_accel,obj.max_accel];
			obj.qpsolve.u_lim = obj.u_lim;

			obj.k1 = params.kp_z;
			obj.k2 = params.kd_z;
			obj.A  = [zeros(2,'double'), eye(2,'double');-obj.k1*eye(2,'double'), -obj.k2*eye(2,'double')];
			obj.qpsolve.update_ricatti(obj.A);
			obj.K  = [obj.k1*eye(2,'double'), obj.k2*eye(2,'double')];
			obj.max_error = params.max_error;

			obj.clf.epsilon = params.clf_epsilon;
			obj.measurement_noise = params.measurement_noise;

			obj.qpsolve.u_cost = params.qp_u_cost;
			obj.qpsolve.u_prev_cost = params.qp_u_prev_cost;
			obj.qpsolve.p1_cost = params.qp_p1_cost;
			obj.qpsolve.p2_cost = params.qp_p2_cost;
			obj.qpsolve.verbose = params.qp_verbose;
			obj.qpsolve.ksig = params.qp_ksig;
			obj.qpsolve.max_var = params.qp_max_var;

			obj.dt = params.dt;

			%%%% USE ROSE - UNTRASLATE
				% obj.model.N_data = params.N_data;
				% obj.model.verbose = params.learning_verbose;
				% obj.model.N_updates = params.N_updates;
				% obj.model.config.meta_batch_size = params.meta_batch_size;
				% obj.model.config.data_horizon = params.data_horizon;
				% obj.model.config.test_horizon = params.test_horizon;
				% obj.model.config.learning_rate = params.learning_rate;
				% obj.model.config.min_datapoints = params.min_datapoints;
				% obj.model.config.save_data_interval = params.save_data_interval;
			%%%% USE ROSE - UNTRASLATE - END
		end
		
		function obj = update_barrier_locations(obj,x,y,radius)
			obj.barrier_locations('x') = x;
			obj.barrier_locations('y') = y;
			obj.barrier_radius = radius;
		end
		
		function result = saturate(obj,u,ulim)
			u0=u(1);
			u0(u0<ulim(1,1))=ulim(1,1);
			u0(u0>ulim(1,2))=ulim(1,2);
			
			u1=u(2);
			u1(u1<ulim(2,1))=ulim(2,1);
			u1(u1>ulim(2,2))=ulim(2,2);
		
			result = [u0,u1];
		end


		
		% untested
		function obj = update_barriers(obj)
			if obj.params.use_barrier_vel
				cbf_list = {BarrierAckermannVelocityZ(true, obj.params.max_velocity, obj.params.barrier_vel_gamma), ...
							BarrierAckermannVelocityZ(false, obj.params.min_velocity, obj.params.barrier_vel_gamma)};
			else
				cbf_list = {};
			end
			

			if obj.params.use_barrier_pointcloud
				bar_loc_x = obj.barrier_locations('x');
				bar_loc_y = obj.barrier_locations('y');
				bar_rad = obj.barrier_radius;
				for i = 1:length(bar_loc_x)
					cbf_list = {cbf_list{:},BarrierAckermannPointZ(bar_loc_x(i),bar_loc_y(i), bar_rad, obj.params.barrier_pc_gamma_p, obj.params.barrier_pc_gamma)};
				end
			end
			obj.qpsolve.cbf_list = cbf_list;
			
		end
		
		% untested
		function obj = get_control(obj,z,z_ref,z_ref_dot,dt,obs,use_model,add_data,use_qp)
			% z,z_ref,z_ref_dot,dt,obs=None,train=False,use_model=False,add_data=True,check_model=True,use_qp=True
			train = false;
			check_model = true;
			obj = obj.update_barriers();

			obj.z = z;
			obj.z_ref = z_ref;
			obj.obs = obs;

			mu_ad = zeros(obj.xdim/2,1);
			mDelta = zeros(obj.xdim/2,1);
			sigDelta = zeros(obj.xdim/2,1);
			rho = zeros(obj.xdim/2,1);
			trueDelta = zeros(obj.xdim/2,1);

			e = obj.z_ref(1:end-1,:)-obj.z(1:end-1,:);
			mu_pd = obj.K*e;
			mu_pd(mu_pd>obj.max_error)=obj.max_error;
			mu_pd(mu_pd<-obj.max_error)=-obj.max_error;

			obj.z_ref_dot = z_ref_dot;
			mu_rm = obj.z_ref_dot(3:end-1);
			
			mu_model = obj.dyn.g(obj.z_prev) * obj.u_prev + obj.dyn.f(obj.z_prev);
			mu_model(mu_model>obj.max_error)=obj.max_error;
			mu_model(mu_model<-obj.max_error)=-obj.max_error;
			
			%%%% USE ROSE - UNTRASLATE
				if add_data
					req = {};
					req.x_next = obj.z;
					req.x = obj.z_prev;
					req.mu_model = mu_model;
					req.obs = obj.obs_prev;
					req.dt = dt;
					obj.model = obj.model.add_data(req);
				end
			%%%% USE ROSE - UNTRASLATE - END
			
			obj.z_dot = (obj.z(3:end-1,:)-obj.z_prev(3:end-1,:))/dt - mu_model;
			
			% # if check_model and obj.model.model_trained:
			if check_model && obj.model_trained
				% # check how the model is doing.  compare the model's prediction with the actual sampled data.
				predict_service_success = false;
				result = NaN;
				%%%% USE ROSE - UNTRASLATE
					req = {};
					req.x = obj.z_prev;
					req.obs = obj.obs_prev;
					result = obj.model.predict(req);
					predict_service_success = true;
				%%%% USE ROSE - UNTRASLATE - END

				if predict_service_success
					obj.y_out = [result.y_out]';
					obj.predict_var = [result.var]';
					obj.predict_error = norm((obj.y_out - obj.z_dot),'fro');
				end
			end
			
			if use_model && obj.model_trained
				predict_service_success = false;
				result = NaN;
				%%%% USE ROSE - UNTRASLATE
					req = {};
					req.x = obj.z;
					req.obs = obj.obs;
					result = obj.model.predict(req);
					predict_service_success = true;
				%%%% USE ROSE - UNTRASLATE - END


				if predict_service_success
					mDelta = [result.y_out]';
					sigDelta = [result.var]';
					% log error if true system model is available
					trueDelta = obj.true_dyn.f(obj.z) - obj.dyn.f(obj.z);
					obj.true_predict_error = norm((trueDelta - mDelta),'fro');

					rho(1) = obj.measurement_noise / (obj.measurement_noise + norm(sigDelta,'fro') + 1e-6);
					mu_ad = (mDelta * rho(1))';
					sigDelta = (sigDelta * (1.0 - rho(1)))';
				end
			else
				sigDelta = ones(obj.xdim/2,1);
			end

			mu_d = mu_rm + mu_pd - mu_ad;

			obj.mu_qp = zeros(obj.xdim/2,1);
			if use_qp
				obj.qpsolve = obj.qpsolve.solve(obj.z,obj.z_ref,mu_d,sigDelta);
                obj.mu_qp = obj.qpsolve.result;
			end
		
			obj.mu_new = mu_d + obj.mu_qp;
			
			obj.u_new = inv(obj.dyn.g(obj.z)) * (obj.mu_new-obj.dyn.f(obj.z));

			u_new_unsaturated = obj.u_new;
			% #saturate in case constraints are not met
			obj.u_new = obj.saturate(obj.u_new,obj.u_lim);

			obj.mu_prev = obj.mu_new;
			obj.u_prev_prev = obj.u_prev;
			obj.u_prev = obj.u_new';
			obj.obs_prev = obj.obs;
			obj.z_prev = obj.z;

			obj.controls = zeros(obj.udim,1);
			obj.controls(1) = atan(obj.u_new(1) * obj.vehicle_length);
			obj.controls(2) = obj.u_new(2);
			
			% return controls
		end	
	end
end
	
			
			
		
		