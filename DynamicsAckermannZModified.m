classdef  DynamicsAckermannZModified
	properties
		xdim
        udim
		epsilon
		disturbance_scale_pos
		disturbance_scale_vel
		control_input_scale
	end
	methods
		% test successful
		function obj = DynamicsAckermannZModified(disturbance_scale_pos,disturbance_scale_vel,control_input_scale)
			obj.xdim=4;obj.udim=2;
			obj.epsilon = 1e-6;
			obj.disturbance_scale_pos = disturbance_scale_pos;
			obj.disturbance_scale_vel = disturbance_scale_vel;
			obj.control_input_scale = control_input_scale;
		end
		
		% test successful
		function result = f(obj,z) 
			v = sqrt(z(2,:)^2 + z(3,:)^2) * z(4,:);
			theta = atan2(z(4)*z(5),z(3)*z(5));
			v_disturbance_body = [tanh(v^2)*obj.disturbance_scale_vel, (0.1+v)*obj.disturbance_scale_vel];
			v_disturbance_world = [v_disturbance_body(1) * cos(theta) - v_disturbance_body(2) * sin(theta), ...
								   v_disturbance_body(1) * sin(theta) + v_disturbance_body(2) * cos(theta)];
            result = [v_disturbance_world(1); v_disturbance_world(2)];
		end
		
		% test successful
		function result = g(obj,z)
			v=(sqrt(z(3,:)^2 + z(4,:)^2)  + obj.epsilon) * z(5,:);
			result = obj.control_input_scale * [-z(4,:)*v,z(3,:)/v;z(3,:)*v,z(4,:)/v];
		end
		
		% test successful
		function znext = step(obj,z,u,dt) 
			v=sqrt(z(3,:)^2 + z(4,:)^2);
			znext = zeros(obj.xdim+1,1);
			znext(1:2,:) = z(1:2,:) + dt*(z(3:4,:) + [sin(v^3)*obj.disturbance_scale_pos; cos(-v)*obj.disturbance_scale_pos]);
			znext(3:4,:) = z(3:4,:) + dt*(obj.f(z) + obj.g(z) * u);
			znext(length(znext)) = z(length(z));
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