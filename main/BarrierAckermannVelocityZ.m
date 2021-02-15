classdef BarrierAckermannVelocityZ
	properties
		dim
		gamma
		bound_from_above
		v_lim
	end
	methods
		function obj = BarrierAckermannVelocityZ(bound_from_above,v_lim,gamma)
			obj.dim=4;
			obj.gamma=gamma;
			obj.bound_from_above=bound_from_above;
			obj.v_lim=v_lim;
		end
		
		function hx = h(obj,z)
			sgn = 1.0;
			if obj.bound_from_above
				sgn = -1.0;
			end
			hx = sgn * (sqrt(z(3)^2 + z(4)^2) * z(5) - obj.v_lim);
		end
		
		function dh = dh(obj,z)
			sgn = 1.0;
			if obj.bound_from_above
				sgn = -1.0;
			end
			z1 = z(1,:);
			z2 = z(2,:);
			z3 = z(3,:);
			z4 = z(4,:);
			z5 = z(5,:);
			v = sqrt(z3^2 + z4^2) + 1.0e-6;
			dh = [0; 0; (z3*z5)/v; (z4*z5)/v] * sgn;
			
		end
		
		function [hx, dB, d2B]=get_B_derivatives(obj,x)
			hx = obj.h(x);
			if hx == 0
				hx = 1e-3;
			end
			dh = obj.dh(x);
			dB = -1.0/(hx*hx)*dh';
			d2B= 2.0/(hx*hx*hx)*(dh(3:length(dh))*dh(3:length(dh))');		
		end

	end
end