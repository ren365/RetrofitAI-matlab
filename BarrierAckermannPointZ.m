classdef BarrierAckermannPointZ
	properties
		dim
		gamma
		x
		y
		radius
		gamma_p
	end
	methods
		function obj = BarrierAckermannPointZ(x,y,radius,gamma_p,gamma)
			obj.dim=4;
			obj.gamma=gamma;
			obj.x=x;
			obj.y=y;
			obj.radius=radius;
			obj.gamma_p=gamma_p;
		end
		
		function hx = h(obj,z)
			sgn = 1.0;
			if obj.radius < 0.0
				sgn = -1.0;
			end
			d = sqrt((z(1) - obj.x)*(z(1) - obj.x) + (z(2) - obj.y)*(z(2) - obj.y)) + 1.0e-6;
			hx = sgn * (obj.gamma_p * (d - obj.radius) + (z(1) - obj.x) / d * z(3) + (z(2) - obj.y) / d * z(4));
		end
		
		function dh = dh(obj,z)
			sgn = 1.0;
			if obj.radius < 0.0
				sgn = -1.0;
			end
			d = sqrt((z(1) - obj.x)*(z(1) - obj.x) + (z(2) - obj.y)*(z(2) - obj.y)) + 1.0e-6;
			d_2 = d*d;
			d_3 = d*d*d;
			y_pos_m_z2 = (obj.y - z(2));
			x_pos_m_z1 = (obj.x - z(1));
			dh = sgn * double([
				(z(3)*(d_2 - x_pos_m_z1*x_pos_m_z1) - obj.gamma_p * d_2 *x_pos_m_z1 - z(4)*x_pos_m_z1*y_pos_m_z2)/d_3,
				(z(4)*(d_2 - y_pos_m_z2*y_pos_m_z2) - obj.gamma_p * d_2 *y_pos_m_z2 - z(3)*x_pos_m_z1*y_pos_m_z2)/d_3,
				-x_pos_m_z1/d,
				-y_pos_m_z2/d]);
			
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