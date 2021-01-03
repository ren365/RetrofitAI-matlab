% test dynamics matlab

cbf=BarrierAckermannVelocityZ(true, 2.0, 10.0);

[a,b,c] = cbf.get_B_derivatives([0.1;0.2;0.3;0.4;0.5])


cbf = BarrierAckermannPointZ(3,5, 1.0, 5.0, 1.0)
% x=3,y=5, radius=1.0, gamma_p=5.0, gamma=1.0
[a,b,c] = cbf.get_B_derivatives([0.1;0.2;0.3;0.4;0.5])