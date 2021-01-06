clear

odim = 2;
adaptive_clbf = AdaptiveClbf(odim);
adaptive_clbf_qp = AdaptiveClbf(odim);
adaptive_clbf_ad = AdaptiveClbf(odim);
adaptive_clbf_pd = AdaptiveClbf(odim);

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

% # vanilla nn params
params.qp_ksig = 1.0e2;
params.measurement_noise = 1.0;

% params.N_data = 600;
% params.learning_verbose = false;
% params.N_updates = 50;
% params.meta_batch_size = 50;
% params.data_horizon = 20;
% params.test_horizon = 30;
% params.learning_rate = 0.001;
% params.min_datapoints = 50;
% params.save_data_interval = 10000;

true_dyn = DynamicsAckermannZModified(0.0, -1.0, 1.0);

adaptive_clbf = adaptive_clbf.update_params(params);
adaptive_clbf_qp = adaptive_clbf_qp.update_params(params);
adaptive_clbf_ad = adaptive_clbf_ad.update_params(params);
adaptive_clbf_pd = adaptive_clbf_pd.update_params(params);
adaptive_clbf.true_dyn = true_dyn;
adaptive_clbf_ad.true_dyn = true_dyn;

barrier_x = [5,15,25,35,45,55];
barrier_y = [0,-0.5,0.5,-0.5,0.5,0];
adaptive_clbf = adaptive_clbf.update_barrier_locations(barrier_x,barrier_y,params.barrier_radius);
adaptive_clbf_qp = adaptive_clbf_qp.update_barrier_locations(barrier_x,barrier_y,params.barrier_radius);
adaptive_clbf_ad = adaptive_clbf_ad.update_barrier_locations(barrier_x,barrier_y,params.barrier_radius);
adaptive_clbf_pd = adaptive_clbf_pd.update_barrier_locations(barrier_x,barrier_y,params.barrier_radius);

x0=[[0.0];[0.0];[0.0];[0.0001]];
z0 = true_dyn.convert_x_to_z(x0);

T = 40;
dt = 0.1;
N = floor(round(T/dt));
t = linspace(0,T-2*dt,N-1);
xdim=4;
udim=2;

train_interval = 40;
start_training = 100;

width = 1.0;
speed = 1.0;
freq = 1.0/10;
x_d = [t * speed; width * sin(2 * pi * t * freq);zeros(1,N-1); zeros(1,N-1)];
x_d(3,1:end-1) = atan2(diff(x_d(2,:)),diff(x_d(1,:)));
x_d(4,1:end-1) = sqrt(diff(x_d(1,:)).^2 + diff(x_d(2,:)).^2)/dt;
x_d(3,end)=x_d(3,end-1);
x_d(4,end)=x_d(4,end-1);

u = zeros(udim,N);
x = zeros(xdim,N-1);
x(:,1) = x0;

u_qp = zeros(udim,N);
x_qp = zeros(xdim,N-1);
x_qp(:,1) = x0;

u_pd = zeros(udim,N);
x_pd = zeros(xdim,N-1);
x_pd(:,1) = x0;

u_ad = zeros(udim,N);
x_ad = zeros(xdim,N-1);
x_ad(:,1) = x0;

z_d = zeros(xdim+1,N-1);
z = zeros(xdim+1,N-1);
z(:,1) = z0;
z_qp = zeros(xdim+1,N-1);
z_qp(:,1) = z0;
z_pd = zeros(xdim+1,N-1);
z_pd(:,1) = z0;
z_ad = zeros(xdim+1,N-1);
z_ad(:,1) = z0;

z_d_dot = zeros(xdim+1,1);

prediction_error_ad = zeros(1,N);
prediction_error_true_ad = zeros(1,N);
prediction_error = zeros(1,N);
prediction_error_true = zeros(1,N);
prediction_var = zeros(xdim/2,N);
prediction_var_ad = zeros(xdim/2,N);
trGssGP = zeros(1,N);

i=0;
z_d(:,i+2) = true_dyn.convert_x_to_z(x_d(:,i+2));

% waitBar = waitbar(0,'Simulating ...');;
tic;
for i = 1:N-2
	% waitbar(i/(N-2),waitBar);

	if i < N-3
		z_d(:,i+2) = true_dyn.convert_x_to_z(x_d(:,i+2));
		z_d_dot = (z_d(:,i+2) - z_d(:,i+1))/dt;
	end

	if i == 0
		add_data = false;
	else
		add_data = true;
    end
	adaptive_clbf = adaptive_clbf.get_control(z(:,i),z_d(:,i+1),z_d_dot,dt,[],true,add_data,true);
	u(:,i+1) = adaptive_clbf.controls;
	% if (i - start_training -1 ) 
		% adaptive_clbf.model.train();
		% adaptive_clbf.model_trained = true;
	% prediction_error(i) = adaptive_clbf.predict_error;
	% prediction_error_true(i) = adaptive_clbf.true_predict_error;
	% prediction_var(:,i:i+1) = np.clip(adaptive_clbf.predict_var,0,params.qp_max_var);
	% trGssGP(i) = adaptive_clbf.qpsolve.trGssGP;

	adaptive_clbf_ad = adaptive_clbf_ad.get_control(z_ad(:,i),z_d(:,i+1),z_d_dot,dt,[],true,add_data,false);
	u_ad(:,i+1) = adaptive_clbf_ad.controls;
	% if (i - start_training - 1) % train_interval == 0 and i > start_training:;
		% adaptive_clbf_ad.model.train();
		% adaptive_clbf_ad.model_trained = true;
	% prediction_error_ad(i) = adaptive_clbf_ad.predict_error;
	% prediction_error_true_ad(i) = adaptive_clbf_ad.true_predict_error;
	% prediction_var_ad(:,i:i+1) = np.clip(adaptive_clbf_ad.predict_var,0,params.qp_max_var);

	adaptive_clbf_qp = adaptive_clbf_qp.get_control(z_qp(:,i),z_d(:,i+1),z_d_dot,dt,[],true,add_data,true);
	u_qp(:,i+1) = adaptive_clbf_qp.controls;
	
	adaptive_clbf_pd = adaptive_clbf_pd.get_control(z_pd(:,i),z_d(:,i+1),z_d_dot,dt,[],false,false,false);
	u_pd(:,i+1) = adaptive_clbf_pd.controls;

	c = u(:,i+1);
	c_ad = u_ad(:,i+1);
	c_qp = u_qp(:,i+1);
	c_pd = u_pd(:,i+1);

	c(1) = tan(c(1))/params.vehicle_length;
	c_ad(1) = tan(c_ad(1))/params.vehicle_length;
	c_qp(1) = tan(c_qp(1))/params.vehicle_length;
	c_pd(1) = tan(c_pd(1))/params.vehicle_length;

	z(:,i+1) = true_dyn.step(z(:,i),c,dt);
	z_ad(:,i+1) = true_dyn.step(z_ad(:,i),c_ad,dt);
	z_qp(:,i+1) = true_dyn.step(z_qp(:,i),c_qp,dt);
	z_pd(:,i+1) = true_dyn.step(z_pd(:,i),c_pd,dt);

	x(:,i+1) = true_dyn.convert_z_to_x(z(:,i+1));
	x_ad(:,i+1) = true_dyn.convert_z_to_x(z_ad(:,i+1));
	x_qp(:,i+1) = true_dyn.convert_z_to_x(z_qp(:,i+1));
	x_pd(:,i+1) = true_dyn.convert_z_to_x(z_pd(:,i+1));

	toc % calculate the time;

end

% fig = plt.figure();
% plt.rcParams.update({'font.size': 12});
% plt.plot(x_d(0,:),x_d(1,:),'k-',label='ref');
% plt.plot(x_ad(0,:),x_ad(1,:),'m--',alpha=0.9,label='ad');
% plt.plot(x_qp(0,:),x_qp(1,:),'b-',alpha=0.9,label='qp');
% plt.plot(x_pd(0,:),x_pd(1,:),'y:',alpha=0.9,label='pd');
% plt.plot(x(0,:),x(1,:),'g-',alpha=0.9,label='balsa',linewidth=3.0);
% plt.legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower center", ncol=5);
% ax = fig.gca();
% for i in range(barrier_x.size):;
	% circle = plt.Circle((barrier_x(i),barrier_y(i)), params.barrier_radius, color='r');
	% ax.add_artist(circle);
% plt.xlabel('X Position');
% plt.ylabel('Y Position');