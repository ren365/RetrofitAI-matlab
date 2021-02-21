function x_d = FSM(start_position,end_position,barrier_x,barrier_y)
	% What is FSM programming
	% For example, if there is two status 1 and 2, the input can be 1 or 2.
	% If we are in status 1 and now the input is 1, then status change to 2 and output 1.
	% If we are in status 1 and now the input is 2, then status stays  in 1 and output 0.
	% If we are in status 2 and now the input is 1, then status change to 1 and output 1.
	% If we are in status 2 and now the input is 2, then status stays  in 2 and output 0.
	% 
	% We can draw a FSM graph about the above logics.
	%
	% Obviously, the first way to program the FSM is by "IF ELSE" method
	% Alternatively, we can use matrix to "hard" program the same logics
	% For example, we still use the above cases.
	% In the status table "S", we have:
	% input     =1 , =2
	% state=1: [ 2 , 1 ]
	% state=2: [ 1 , 2 ]
	% next_status = S(current_status,input_value)
	%
	% For the output table "O", we have:
	% input     =1 , =2
	% state=1: [ 1 , 0 ]
	% state=2: [ 1 , 0 ]
	% current_output = O(current_status,input_value)
	
	% what you need to do
	
	% step 1: generate proper map by the location of barriers and start/end point
	% not too large or small!
	todo = Nan;
	% step 2: build your own FSM to run logicly, the only things is that 
	% avoid barriers and reach the end!
	todo = Nan;
	% step 3: now we have the FSM to run, however, there is no difference among 
	% non-barrier neighbors. In such case, we need to find out a way to go to the next
	% right direction towards the end
	todo = Nan;
	% step 4: After simulating the above FSM, we have a reference path now.
	% add direction "theta" and velocity to fill up the x_d!
	todo = Nan;
end