% dynamics
function result = barrier_dynamics_linear(barriers,dynamics,time)
	result = barriers + dynamics * time;
end