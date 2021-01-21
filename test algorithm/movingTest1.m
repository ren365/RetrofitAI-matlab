% prototype of initial barrier dynamics / 
S = dbstack();
filename = S(1).file(1:end-2);
% x,y
barriers = [0 ,   3, ... % upper vehicle 
			0 ,-  3, ... % lower vehicle 
			10, -1, ... % pedestrian 1
			20,  2, ... % pedestrian 2
			30, -3];	 % pedestrian 3
barriers_radius = [1,1,1,1,1];
% barriers_radius = [1.2,1.2,1.2,1.2,1.2];
barriers_numbers = length(barriers_radius);
% dx, dy
dynamics = [0.9,0, ...
			1.1,0, ...
			0  ,0.1, ...
			0  ,-0.1, ...
			0  ,0.1];
% barriers location recording
barriers_locations = barriers;

% original codes
odim = 2;
adaptive_clbf = AdaptiveClbf(odim);

params={};
params.vehicle_length = 0.25;
params.steering_limit = 0.75;% 0.75
params.max_accel = 1.0; % 1.0
params.min_accel = -1.0;% -1.0
params.kp_z = 1.0;
params.kd_z = 1.0;
params.clf_epsilon = 100.0;


params.qp_u_cost = 100.0;
params.qp_u_prev_cost = 1.0;
params.qp_p1_cost = 1.0;
params.qp_p2_cost = 1.0e12;
params.qp_max_var = 1.5;
params.qp_verbose = false;
params.max_velocity = 2.0; % used for CBF, BarrierAckermannVelocityZ
params.min_velocity = 0.5; % used for CBF, BarrierAckermannVelocityZ
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

true_dyn = DynamicsAckermannZModified(0.0, -1.0, 1.0);

adaptive_clbf = adaptive_clbf.update_params(params);
adaptive_clbf.true_dyn = true_dyn;

x0=[[0.0];[0.0];[0.0];[0.0001]];
z0 = true_dyn.convert_x_to_z(x0);

T = 40;
dt = 0.1;
N = floor(round(T/dt));
t = linspace(0,T-2*dt,N-1);
xdim=4;
udim=2;

train_interval = 40;
start_training = 30;

width = 1.0;
speed = 1.0;
freq = 1.0/10;
x_d = [t * speed; zeros(1,N-1);zeros(1,N-1); zeros(1,N-1)];
x_d(3,1:end-1) = atan2(diff(x_d(2,:)),diff(x_d(1,:)));
x_d(4,1:end-1) = sqrt(diff(x_d(1,:)).^2 + diff(x_d(2,:)).^2)/dt;
x_d(3,end)=x_d(3,end-1);
x_d(4,end)=x_d(4,end-1);

u = zeros(udim,N);
x = zeros(xdim,N-1);
x(:,1) = x0;

z_d = zeros(xdim+1,N-1);
z = zeros(xdim+1,N-1);
z(:,1) = z0;
z_d_dot = zeros(xdim+1,1);


prediction_error = zeros(1,N);
prediction_error_true = zeros(1,N);
prediction_var = zeros(xdim/2,N);
trGssGP = zeros(1,N);

i=0;
z_d(:,i+2) = true_dyn.convert_x_to_z(x_d(:,i+2));
add_data = false;

waitBar = waitbar(0,'Simulating ...');

distance_to_barriers = [];
qpStatus = [];

tic;
for i = 1:N-2
	
	waitbar(i/(N-2),waitBar);

	% update the barriers
	barriers_locations_next = barrier_dynamics_linear(barriers_locations(end,:),dynamics,dt);
	barriers_locations =[barriers_locations;barriers_locations_next];
	% current time step the barrier location
	adaptive_clbf = adaptive_clbf.update_barrier_locations(barriers_locations(i,1:2:end-1), ...
		barriers_locations(i,2:2:end),params.barrier_radius);

	% end udpate barriers

	if i < N-2
		z_d(:,i+2) = true_dyn.convert_x_to_z(x_d(:,i+2));
		z_d_dot = (z_d(:,i+2) - z_d(:,i+1))/dt;
	end

	if i > 1
		add_data = true;
    end
	
	adaptive_clbf = adaptive_clbf.get_control(z(:,i),z_d(:,i+1),z_d_dot,dt,[x(3,i),u(:,i)'],true,add_data,true);
	u(:,i+1) = adaptive_clbf.controls;
	qpStatus = [qpStatus,adaptive_clbf.qpsolve.res.info.status_val];
	if mod((i-1 - start_training - 1) , train_interval) == 0 && i > start_training
		adaptive_clbf.model = adaptive_clbf.model.train();
		adaptive_clbf.model_trained = true;
	end

	c = u(:,i+1);
	c(1) = tan(c(1))/params.vehicle_length;
	z(:,i+1) = true_dyn.step(z(:,i),c,dt);
	x(:,i+1) = true_dyn.convert_z_to_x(z(:,i+1));
	
	% next time step the barrier location
	current_distance = nearest_distance(barriers_locations(i+1,:),x(:,i), ...
									barriers_radius,barriers_numbers); % barriers_radius-0.2 = actual circle
	distance_to_barriers = [distance_to_barriers,current_distance];

	toc % calculate the time;

end

% distance to barrier
subplot(2,1,1)
plot(t(1:end-1),distance_to_barriers,'g-')
hold on
plot([t(1),t(end)],[0,0],'r--')
ylabel("Dist to barrier (m)")
subplot(2,1,2)
plot(t(1:end-1),qpStatus,'k-')
ylabel("QP status")
xlabel("Time(s)")
saveas(gcf,[filename,'-distance.png'])
close all

figure
subplot(2,1,1)
plot(t,u(1,1:end-1),'b-')
hold on
plot([t(1),t(end)],[params.steering_limit,params.steering_limit],'r--')
plot([t(1),t(end)],[-params.steering_limit,-params.steering_limit],'r--')
axis([t(1),t(end),-params.steering_limit-0.5,params.steering_limit+0.5]);
ylabel('Steering Angle (rad)')
subplot(2,1,2)
plot(t,u(2,1:end-1),'b-')
hold on
plot([t(1),t(end)],[params.min_accel,params.min_accel],'r--')
plot([t(1),t(end)],[params.max_accel,params.max_accel],'r--')
axis([t(1),t(end),params.min_accel-0.5,params.max_accel+0.5]);
ylabel('Throttle (m/s^2)')
xlabel('Time (s)')
saveas(gcf,[filename,'-control.png'])
close all

figure
if true
	% figure('visible','off');
	waitBar = waitbar(0,'generating GIF ...');
	
	xlabel('X Position');
	ylabel('Y Position');
	box on;
	for i = 1:N-1
		waitbar(i/(N-1),waitBar);
		cla
		axis([0,40,-5,5]);
		hold on
		% circles(barrier_x,barrier_y,radius,'facecolor','cyan','edgecolor','cyan');
		plot(x_d(1,:),x_d(2,:),'k--');
		plot(x(1,i),x(2,i),'go');
		% barriers
		barrier_x=[];barrier_y=[];radius=[];
		for j=1:barriers_numbers
			xtmp = barriers_locations(i,j*2-1);
			ytmp = barriers_locations(i,j*2);
			barrier_x = [barrier_x;xtmp];
			barrier_y = [barrier_y;ytmp];
			radius	  = [radius   ;barriers_radius(j)];
		end

		circles(barrier_x,barrier_y,radius,'facecolor','red','edgecolor','red');
		yline(-1.5)
		yline(-4.5)
		yline(1.5)
		yline(4.5)
		legend({'reference path','paper''s method','cars'},'location','northeastoutside');

		frame=getframe(gcf);
		imind=frame2im(frame);
		[imind,cm] = rgb2ind(imind,256);
		if i==1
			 imwrite(imind,cm,[filename,'-move.gif'],'gif', 'Loopcount',inf,'DelayTime',0);
		else
			 imwrite(imind,cm,[filename,'-move.gif'],'gif','WriteMode','append','DelayTime',0);
		end
	end
	close all
end