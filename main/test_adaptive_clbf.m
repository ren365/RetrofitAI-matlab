clear

% Barriers
barrier_x = [10,10,10,10]; % feel free to change the barrier location
barrier_y = [3,1,-1,-3];   % if needed

% Time 
T = 30;
dt = 0.1;
N = floor(round(T/dt));
t = linspace(0,T-2*dt,N-1);

% settings updates
odim = 2;
xdim = 4;
udim = 2;
params = System_settings();

adaptive_clbf = AdaptiveClbf(odim);
true_dyn = DynamicsAckermannZModified(0.0, -1.0, 1.0);
adaptive_clbf = adaptive_clbf.update_params(params);
adaptive_clbf.true_dyn = true_dyn;
adaptive_clbf = adaptive_clbf.update_barrier_locations(barrier_x,barrier_y,params.barrier_radius);


useFSM = false; % change this line to true to test your FSM function
if useFSM
	% FSM
	start_position = [0,0];
	end_position = [10,0];
	map = Nan;% create your own map based on Barriers
	x_d = FSM(start_position,end_position,barrier_x,barrier_y);
	% end FSM
else
	% original settings for reference
	width = 1.0;
	speed = 1.0;
	freq = 1.0/10;

	x_d = [t * speed; width * sin(2 * pi * t * freq);zeros(1,N-1); zeros(1,N-1)];
	x_d(3,1:end-1) = atan2(diff(x_d(2,:)),diff(x_d(1,:)));
	x_d(4,1:end-1) = sqrt(diff(x_d(1,:)).^2 + diff(x_d(2,:)).^2)/dt;
	x_d(3,end)=x_d(3,end-1);
	x_d(4,end)=x_d(4,end-1);
end

% location & Path & other paras
x0=[[0.0];[0.0];[0.0];[0.0001]];
z0 = true_dyn.convert_x_to_z(x0);
u = zeros(udim,N);
x = zeros(xdim,N-1);
x(:,1) = x0;

z_d = zeros(xdim+1,N-1);
z = zeros(xdim+1,N-1);
z(:,1) = z0;

z_d_dot = zeros(xdim+1,1);


% main to run
waitBar = waitbar(0,'Simulating ...');

z_d(:,2) = true_dyn.convert_x_to_z(x_d(:,2)); % init
tic;

for i = 1:N-2
	
	waitbar(i/(N-2),waitBar);

	% corner case
	if i < N-2
		z_d(:,i+2) = true_dyn.convert_x_to_z(x_d(:,i+2));
		z_d_dot = (z_d(:,i+2) - z_d(:,i+1))/dt;
	end

	if i == 1
		add_data = false;
	else
		add_data = true;
    end
	
	% get control
	adaptive_clbf = adaptive_clbf.get_control(z(:,i),z_d(:,i+1),z_d_dot,dt,[x(3,i),u(:,i)'],true,add_data,true);
	u(:,i+1) = adaptive_clbf.controls;
	
	% train the adaptation model
	if mod((i-1 - params.start_training - 1) , params.train_interval) == 0 && i > params.start_training
		adaptive_clbf.model = adaptive_clbf.model.train();
		adaptive_clbf.model_trained = true;
	end

	% went over a step
	c = u(:,i+1);
	c(1) = tan(c(1))/params.vehicle_length;
	z(:,i+1) = true_dyn.step(z(:,i),c,dt);
	x(:,i+1) = true_dyn.convert_z_to_x(z(:,i+1));

	toc % calculate the time;

end

% figure

plot(x_d(1,:),x_d(2,:),'k-','LineWidth',2);
hold on
plot(x(1,:),x(2,:),'g-','LineWidth',3);
radius = ones(length(barrier_x),1)*params.barrier_radius;
circles(barrier_x',barrier_y',radius,'facecolor','red','edgecolor','red');
legend({'reference','only adaptive','circles'},...
		'location','northeastoutside')
xlabel('X Position');
ylabel('Y Position');
