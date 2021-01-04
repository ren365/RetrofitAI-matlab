clear
dyn = DynamicsAckermannZ();
clf = LyapunovAckermannZ(10.0,1.0,1.0,1.0);
u_lim = [-2.0,2.0;-1.0,1.0];

% cbf list
cbf_list = 	{BarrierAckermannVelocityZ(true, 2.0, 10.0), ...
			BarrierAckermannVelocityZ(false, 0.5, 10.0), ...
			BarrierAckermannPointZ(1,0, 3.0, 5.0, 1.0)};
% cbf_list = 	[BarrierAckermannPointZ(1,0, 3.0, 5.0, 1.0)];
% qpsolve
qpsolve = QPSolve(dyn,cbf_list,clf,u_lim,0.0,1.0,1.0e8,1.0e8,false);
qpsolve.u_cost = 100.0;
qpsolve.u_prev_cost = 1.0;
qpsolve.p1_cost = 1.0;
qpsolve.p2_cost = 1.0e12;
qpsolve.verbose = false;
qpsolve.ksig = 1.0e2;
qpsolve.max_var = 1.5;

A=[zeros(2,2), eye(2);-eye(2), -eye(2)];
qpsolve = qpsolve.update_ricatti(A);

sigDelta = [0.35;0.6];
x = zeros(5,4);
x(1,1)=0.1;
x(2,1)=0.2;
x(3,1)=0.3;
x(4,1)=0.4;
x(5,1)=0.5;
% # x[:,0:1];
x(1,2)=0.4;
x(2,2)=0.3;
x(3,2)=0.2;
x(4,2)=0.1;
x(5,2)=0.5;

z_ref_dot = zeros(2,1);
qpsolve = qpsolve.solve(x(:,1),x(:,2),z_ref_dot,sigDelta);
qpsolve.mu_qp_prev