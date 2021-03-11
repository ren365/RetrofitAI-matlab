function x_d = FSM(start_position,end_position,barrier_x,barrier_y,barrier_radius,dt, N)
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
    row_count = max(range(barrier_y)+1, start_position(1), end_position(1));
    col_count = max(range(barrier_x)+1, start_position(1), end_position(1));
    matrix = ones(row_count, col_count)
    
    matrix = set_barrier(matrix, barrier_x, barrier_y, barrier_radius)

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
    
    % original settings for reference
	width = 1.0;
	speed = 1.0;
	freq = 1.0/10;
    
    [pathX pathY] = find_referenceXY(x_pos, y_pos, dt, radius,x,N);
    
	x_d = [pathX; pathY;zeros(1,N-1); zeros(1,N-1)];
	x_d(3,1:end-1) = atan2(diff(x_d(2,:)),diff(x_d(1,:)));
	x_d(4,1:end-1) = sqrt(diff(x_d(1,:)).^2 + diff(x_d(2,:)).^2)/dt;
	x_d(3,end)=x_d(3,end-1);
	x_d(4,end)=x_d(4,end-1);

end

% credit: function to draw circle in the matrix
% https://matlab.fandom.com/wiki/FAQ#How_do_I_create_a_circle.3F
function graph = set_barrier(matrix, barrier_x, barrier_y, barrier_radius)
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
    col_size = size(matrix,2);
    row_size = size(matrix,1);
    radius = barrier_radius;
    [columnsInImage rowsInImage] = meshgrid(1:col_size, 1:row_size);

    graph = zeros(sizeof(matrix));
    for i = 1:length(barrier_x)
        centerX = barrier_x(i);
        centerY = barrier_y_all_positive(i);
        graph = graph | ((rowsInImage - centerY).^2 ...
            + (columnsInImage - centerX).^2 <= radius.^2);
    end
    % negate the graph, since in the search algorithm, 1 can go, 0 cannot
    % however, before the following step, anywhere belongs to barrier
    % is marked as 1.
    graph = ~graph;
end

% credit of how to draw the circle path in matrix
% https://www.mathworks.com/matlabcentral/answers/568719-how-to-plot-a-quarter-of-a-circle
function [pathX, pathY] find_referenceXY(x_pos, y_pos, dt, radius,x, N)
    speed = 1;
    pathX = ones(1,N) * x_pos(end);
    pathY = ones(1,N) * y_pos(end);
    
    point_count = 0;
    for i=1:size(x_pos,2)-1
        % point is include the right/next step, but not the beginning
        if (y_pos(i+1) - y_pos(i) == 0)
            left = min(x_pos(i+1),x_pos(i));
            right = max(x_pos(i+1),x_pos(i));
            dist = right - left;
            point_num = floor(dist/(v*dt));
            
            for p =point_count+1:point_count+point_num
                temp = p-point_count+1;
                pathX(p) = (left+temp*dt);
                pathY(p) = y_pos(i);
            end
            point_count = point_count + point_num;
        elseif (x_pos(i+1) - x_pos(i) == 0)
            up = max(y_pos(i+1), y_pos(i));
            down = min(y_pos(i+1), y_pos(i));
            dist = down - up;
            point_num = floor(dist/(v*dt));
            
            for p=point_count+1:point_count+point_num
                temp = p-point_count+1;
                pathX(p) = x_pos(i);
                pathY(p) = (down+temp*dt);
            end
            point_count = point_count + point_num;
        elseif (abs(x_pos(i+1)-x_pos(i))==radius & abs(y_pos(i+1)-y_pos(i))==radius)
            dist = (pi*radius)/2;
            point_num = floor(dist/(v*dt));
            iteration_interval = 90/(point_num-1);
            if (x_pos(i+1)>x_pos(i)&y_pos(i+1)>y_pos(i))
                % Center
                C = [x_pos(i+1),y_pos(i)];
                % angle 
                th = -90:iteration_interval:0;
            elseif (x_pos(i+1)>x_pos(i)&y_pos(i+1)<y_pos(i))
                % Center
                C = [x_pos(i),y_pos(i+1)];
                % angle 
                th = 0:iteration_interval:90;
            elseif (x_pos(i+1)<x_pos(i)&y_pos(i+1)>y_pos(i))
                % Center
                C = [x_pos(i+1),y_pos(i)];
                % angle 
                th = 90:-iteration_interval:0;
            elseif (x_pos(i+1)<x_pos(i)&y_pos(i+1)<y_pos(i))
                % Center
                C = [x_pos(i),y_pos(i+1)];
                % angle 
                th = 0:-iteration_interval:-90;
            end
            % points of cricle 
            xc = C(1)+R*sind(th) ; 
            yc = C(2)+R*cosd(th) ;
            for p=1:size(xc,2)
                index = point_count+p;
                pathX(index) = xc(p);
                pathY(index) = yc(p);
            end
            point_count = point_count + size(xc,2);
            
            if size(xc,2) < size(point_num)
                diff = size(point_num) - size(xc,2);
                for d=1:diff
                    pathX(point_count+d) = x_pos(i+1);
                    pathY(point_count+d) = y_pos(i+1);
                end
                point_count = point_count + diff;
            end
        end

    end

   
    