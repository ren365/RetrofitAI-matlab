classdef  DynamicsAckermannZ
	properties
		xdim
        udim
		epsilon
	end
	methods
		% test successful
		function obj = DynamicsAckermannZ()
			obj.xdim=4;obj.udim=2;
			obj.epsilon = 1e-6;
		end
		
		%
		function result = f(obj,z) 
            result = [0; 0]*z(5,:);
		end
		
		% 
		function result = g(obj,z)
			v=(sqrt(z(3,:)^2 + z(4,:)^2)  + obj.epsilon) * z(5,:);
			result = [-z(4,:)*v,z(3,:)/v;z(3,:)*v,z(4,:)/v];
		end
		
		% test successful
		function result = convert_z_to_x(obj,z) 
			v = sqrt(z(3,:)^2 + z(4,:)^2) * z(5,:);
			result = [z(1,:);z(2,:);atan2(z(5,:)*z(4,:),z(5,:)*z(3,:));v];
		end
		
		% test successful
		function result = convert_x_to_z(obj,x)
			v_sign = 1.0;
			if x(4,:) < 0
				v_sign = -1;
			end
			result = [x(1,:);x(2,:);x(4,:)*cos(x(3,:));x(4,:)*sin(x(3,:));v_sign];
		end	
	end
end