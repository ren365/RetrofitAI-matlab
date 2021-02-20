function result = System_settings()
	% para settings
	% feel free to modify the system
	params={};
	params.vehicle_length = 0.25;
	params.steering_limit = 0.75;
	params.max_accel = 1.0;
	params.min_accel = -1.0;
	params.kp_z = 1.0;
	params.kd_z = 1.0;
	params.clf_epsilon = 100.0;


	params.qp_u_cost = 100.0;
	params.qp_u_prev_cost = 1.0;
	params.qp_p1_cost = 1.0;
	params.qp_p2_cost = 1.0e12;
	params.qp_max_var = 1.5;
	params.qp_verbose = false;
	params.max_velocity = 2.0;
	params.min_velocity = 0.5;
	params.barrier_vel_gamma = 10.0;
	params.use_barrier_vel = true;
	params.use_barrier_pointcloud = true;
	params.barrier_radius = 1.0;
	params.barrier_radius_velocity_scale = 0.0;
	params.barrier_pc_gamma_p = 5.0;
	params.barrier_pc_gamma = 1.0;
	params.verbose = false;
	params.dt = 0.1;
	params.max_error = 10.0;

	params.qp_ksig = 1.0e2;
	params.measurement_noise = 1.0;
	
	params.train_interval = 40;
	params.start_training = 100;

end