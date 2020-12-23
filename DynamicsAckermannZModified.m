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
		function obj = DynamicsAckermannZModified(disturbance_scale_pos,disturbance_scale_vel,control_input_scale)
			obj.xdim=4;obj.udim=2;
			obj.epsilon = 1e-6;
			obj.disturbance_scale_pos = disturbance_scale_pos;
			obj.disturbance_scale_vel = disturbance_scale_vel;
			obj.control_input_scale = control_input_scale;
		end
		
		function result = f(obj,z) % what is start of z ???
			v = sqrt(z(2,:)^2 + z(3,:)^2) * z(4,:);% sqrt(z(3,:)^2 + z(4,:)^2) * z(5,:);
			theta = arctan2(z(3)*z(4),z(2)*z(4));% arctan2(z(4)*z(5),z(3)*z(5))
			v_disturbance_body = [tanh(v^2)*obj.disturbance_scale_vel, (0.1+v)*obj.disturbance_scale_vel];
			v_disturbance_world = [v_disturbance_body(1) * cos(theta) - v_disturbance_body(2) * sin(theta), ...
								   v_disturbance_body(1) * sin(theta) + v_disturbance_body(2) * cos(theta)];
            result = [v_disturbance_world(0), v_disturbance_world(1)];
		end
			
		function result = g(obj,z) % what is start of z ???
			v=(sqrt(z(3,:)^2 + z(4,:)^2)  + obj.epsilon) * z(5,:);
			result = obj.control_input_scale * [[-z(3,:)*v,z(2,:)/v],[z(3,:)*v,z(4,:)/v]];
		end
		
		function znext = step(obj,z,u,dt) % what is start of z ???
			v=sqrt(z(3,:)^2 + z(4,:)^2);
			znext = zeros(obj.xdim+1,1);
			znext(0:2,:) = z(1:3,:) + dt*(z(3:5,:) + [sin(v^3)*obj.disturbance_scale_pos, cos(-v)*obj.disturbance_scale_pos]);
			znext(2:-1,:) = z(2:-1,:) + dt*(obj.f(z) + obj.g(z) * u);
			znext(length(znext)) = z(length(z));
		end
		
		function result = convert_z_to_x(obj,z) % what is start of z ???
			v = sqrt(z(3,:)^2 + z(4,:)^2) * z(5,:);
			result = np.stack((z(1,:),z(2,:),arctan2(z(5,:)*z(4,:),z(5,:)*z(3,:)),v),axis=0)
		end
			
		def convert_x_to_z(self,x):
			v_sign = 1.0
			if x[3,:] < 0:
				v_sign = -1
			return np.stack((x[0,:],x[1,:],x[3,:]*np.cos(x[2,:]),x[3,:]*np.sin(x[2,:]),np.array([v_sign])),axis=0)
			
	end
end