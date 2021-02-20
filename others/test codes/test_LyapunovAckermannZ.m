% test dynamics matlab
true_dyn = DynamicsAckermannZModified(0.0, -1.0, 1.0);
clf = LyapunovAckermannZ(10.0,1.0,1.0,1.0);
x0=[[0.0];[0.0];[0.0];[0.0001]];
% z0 = true_dyn.convert_x_to_z(x0);
% x0 = true_dyn.convert_z_to_x(z0);

N=5;
xdim=4;

x0=[[0.0];[0.0];[0.0];[0.0001]];
z0 = true_dyn.convert_x_to_z(x0);
z = zeros(xdim+1,N-1);
z(:,1) = z0;
x = zeros(xdim,N-1);
x(:,1) = x0;
z_ref = ones(xdim+1,N-1);
z_ref(:,1) = z0;

clf.V(z(:,3),z_ref(:,3))