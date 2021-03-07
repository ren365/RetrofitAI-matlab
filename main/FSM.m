function x_d = FSM(start_position,end_position,barrier_x,barrier_y,dt)
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
    offset = min(barrier_y);
    barrier_y_all_positive = barrier_y;
    if offset < 0
        offset = offset * (-1);
        for i = 1:length(barrier_y)
            barrier_y_all_positive(i) = barrier_y + offset;
        end
    else
        offset = 0;
    end
            
    row_count = max(range(barrier_y)+1, start_position(1), end_position(1));
    col_count = max(range(barrier_x)+1, start_position(1), end_position(1));
    matrix = ones(row_count, col_count)
    for i = 1:length(barrier_x)
        matrix(barrier_y_all_positive(i), barrier_x(i)) = 0;
    end
	% step 2: build your own FSM to run logicly, the only things is that 
	% avoid barriers and reach the end!
    path = bfs_search(matrix, start_position, end_position);
    
	% step 3: now we have the FSM to run, however, there is no difference among 
	% non-barrier neighbors. In such case, we need to find out a way to go to the next
	% right direction towards the end
    [minlength, minindex] = min(cellfun('size', path, 2));
    shortest_path = path{minindex}; % cell array
    % change position back to the +/- range
    path_len = size(shortest_path, 2);
    x_pos = ones(1,path_len);
    y_pos = ones(1,path_len);
    for i=1:path_len
        x_pos(1,i) = shortest_path{i}(1);
        y_pos(1,i) = shortest_path{i}(2) - offset;
    end
    
	% step 4: After simulating the above FSM, we have a reference path now.
	% add direction "theta" and velocity to fill up the x_d!
	todo = Nan;
end