function [barrier_x,barrier_y]=moving_barrier(...
			barrier_x,barrier_y,barrier_vx,barrier_vy,dt)
	barrier_x = barrier_x + barrier_vx*dt;
	barrier_y = barrier_y + barrier_vy*dt;
end