classdef  LyapunovAckermannZ
	% untest
	properties
		w1
		w2
		w3
        epsilon
		dim
	end
	methods (Access = public)
		function obj = LyapunovAckermannZ(w1,w2,w3,epsilon)
			obj.dim= 4;
			obj.w1 = w1;
			obj.w2 = w2;
			obj.w3 = w3;
			obj.epsilon = epsilon;
		end
		
		function result = signed_angle_dist(obj,a,b)
			d = a - b;
			if d > pi
				d = d - pi*2;
            elseif d < -pi
				d = d + pi*2;
			end
			result = d;
		end
		
		function [e,phi_t,alpha,theta] = convert_to_polar(obj,x,x_d)

			e = sqrt((x_d(1,:)-x(1,:))^2 + (x_d(2,:)-x(2,:))^2);
			phi_t = atan2(x_d(2,:)-x(2,:),x_d(1,:)-x(1,:));

			phi = atan2(x(4,:),x(3,:));
			phi_d = atan2(x_d(4,:),x_d(3,:));
			
			alpha = obj.signed_angle_dist(phi_t,phi);
			theta = obj.signed_angle_dist(phi_t,phi_d);

		end
		
		
		function result = V(obj,z,z_d)
			[e,phi_t,alpha,theta] = obj.convert_to_polar(z,z_d);
			v = sqrt(z(3,:)^2 + z(4,:)^2);
			v_d = sqrt(z_d(3,:)^2 + z_d(4,:)^2);
			result = 0.5 * (obj.w1*alpha^2 + obj.w2*theta^2 + obj.w3*(v - v_d)^2);
		end

	end
end
