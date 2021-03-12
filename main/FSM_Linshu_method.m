start_position = [0,0];
end_position = [10,0];

barrier_x = [3,3,3,5]; % feel free to change the barrier location
barrier_y = [5,-5,0,3];   % if needed

% Time 
T = 30;
dt = 0.1;
N = floor(round(T/dt));
t = linspace(0,T-2*dt,N-1);

x_d = FSM(start_position,end_position,barrier_x,barrier_y,1,dt,N);


function x_d = FSM(start_pos,end_pos,barrier_x,barrier_y,barrier_radius,dt,N)
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
    start_position = start_pos + 1;
    end_position = end_pos + 1;
    
    row_temp = max(range(barrier_y)+1, start_position(2));
    row_count = max(row_temp, end_position(2)); % position index start at 0
    col_temp = max(range(barrier_x)+1, start_position(1));
    col_count = max(col_temp, end_position(1)); % position index start at 0
    matrix = ones(row_count, col_count);
    [matrix, y_axis] = set_barrier(matrix, barrier_x, barrier_y, barrier_radius);
    
    % convert original matrix to fit the property of graph
    % barrier is Inf, target is -Inf, all the rest are 0
    formatted_matrix = ~matrix;
    formatted_matrix = formatted_matrix * Inf;
    formatted_matrix(isnan(formatted_matrix)) = 0;
    row_position = y_axis - end_pos(2);
    
    
	% step 2: build your own FSM to run logicly, the only things is that 
	% avoid barriers and reach the end!
    % input start and target is [col, row], the position input and output 
    % of bfs search is [row, col]
    
    s = [start_position(2), start_position(1)];
    e = [y_axis - end_pos(2), end_position(1)];
    
    tic
    path = bfs_search(matrix, s, e);
    X = ["for matrix with size ", size(matrix), "excute time for the bfs search is: "]
    disp(X)
    toc
    
	% step 3: now we have the FSM to run, however, there is no difference among 
	% non-barrier neighbors. In such case, we need to find out a way to go to the next
	% right direction towards the end
    [minlength, minindex] = min(cellfun('size', path, 2));
    shortest_path = path{minindex}; % cell array
    disp("Shortest Path is: [row, col] format");
    celldisp(shortest_path);
    % change position back to the +/- range
    path_len = size(shortest_path, 2);
    x_pos = ones(1,path_len);
    y_pos = ones(1,path_len);
    for i=1:path_len
        x_pos(1,i) = shortest_path{i}(2);
        y_pos(1,i) = y_axis - shortest_path{i}(1);
        formatted_matrix(shortest_path{i}(1),shortest_path{i}(2)) = 1;
    end
    formatted_matrix(row_position,end_position(1)) = log(0);
    
    %plot(matrix)
    %hold on
    disp(formatted_matrix);
    disp("x_pos");
    disp(x_pos);
    disp("y_pos");
    disp(y_pos)
    plot(x_pos,y_pos,'k-','LineWidth',2);
    circles(barrier_x',barrier_y',barrier_radius,'facecolor','red','edgecolor','red');
legend({'reference', 'barrier'},...
		'location','northeastoutside')
    
	% step 4: After simulating the above FSM, we have a reference path now.
	% add direction "theta" and velocity to fill up the x_d!
    
    [pathX pathY] = find_referenceXY(x_pos, y_pos, barrier_radius, dt, N);

    
	x_d = [pathX; pathY;zeros(1,N-1); zeros(1,N-1)];
	x_d(3,1:end-1) = atan2(diff(x_d(2,:)),diff(x_d(1,:)));
	x_d(4,1:end-1) = sqrt(diff(x_d(1,:)).^2 + diff(x_d(2,:)).^2)/dt;
	x_d(3,end)=x_d(3,end-1);
	x_d(4,end)=x_d(4,end-1);

end

% credit: function to draw circle in the matrix
% https://matlab.fandom.com/wiki/FAQ#How_do_I_create_a_circle.3F
function [graph, y_axis] = set_barrier(matrix, barrier_x, barrier_y, barrier_radius)
    barrier_y_all_positive = barrier_y;
    barrier_x = barrier_x+1;
    
    y_axis = ceil(size(matrix,1)/2);
    for i = 1:length(barrier_y)
        barrier_y_all_positive(i) = y_axis - barrier_y(i);
    end
    col_size = size(matrix,2);
    row_size = size(matrix,1);
    radius = barrier_radius;
    [columnsInImage, rowsInImage] = meshgrid(1:col_size, 1:row_size);

    graph = zeros(size(matrix));
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
function [pathX, pathY] = find_referenceXY(x_pos, y_pos, radius, dt, N)
    % original settings for reference
	width = 1.0;
	speed = 1.0;
	freq = 1.0/10;

    pathX = ones(1,N) * x_pos(end);
    pathY = ones(1,N) * y_pos(end);
    
    point_count = 0;
    for i=1:size(x_pos,2)-1
        % point is include the right/next step, but not the beginning
        if (y_pos(i+1) - y_pos(i) == 0)
            left = min(x_pos(i+1),x_pos(i));
            right = max(x_pos(i+1),x_pos(i));
            dist = right - left;
            point_num = floor(dist/(speed*dt));
            
            for p = (point_count+1):(point_count+point_num)
                disp("p is")
                disp(p)
                
                temp = p-point_count+1;
                pathX(p) = (left+temp*dt);
                pathY(p) = y_pos(i);
            end
            point_count = point_count + point_num;
        elseif (x_pos(i+1) - x_pos(i) == 0)
            up = max(y_pos(i+1), y_pos(i));
            down = min(y_pos(i+1), y_pos(i));
            dist = down - up;
            point_num = floor(dist/(speed*dt));
            
            for p=point_count+1:point_count+point_num
                disp("p is")
                disp(p)
                
                temp = p-point_count+1;
                pathX(p) = x_pos(i);
                pathY(p) = (down+temp*dt);
            end
            point_count = point_count + point_num;
        elseif (abs(x_pos(i+1)-x_pos(i))==radius & abs(y_pos(i+1)-y_pos(i))==radius)
            dist = (pi*radius)/2;
            point_num = floor(dist/(speed*dt));
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
                disp("p is")
                disp(p)
                
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
end

function all_path = bfs_search(graph, start, target)
% graph: N x M Matrix
% start: [row, col]: 1x2 matrix contain the position of start position
% end: [row, col] 1x2 matrix contain the position of target position
disp("graph is")
disp(graph)
disp("start is")
disp(start)
disp("end is")
disp(target)


visited = graph;
all_path = {};
queue = {{start}};
visited(start(1),start(2)) = 0;
while size(queue,2) ~= 0
    % pop the first path from queue
    path = queue{1};
    queue = queue(:,2:end);
    % get the last position in this path
    last_node = path{end};
    % check if last_node is target
    if (isequal(last_node, target))
        all_path{end+1} = path;
        continue;
    end
    adjacent = find_adjacent(visited,last_node(1),last_node(2),true);
    
    for adj=adjacent
        visited(adj{1}(1),adj{1}(2)) = 0;
        new_path = path;
        new_path{end+1} = [adj{1}(1), adj{1}(2)];
        queue{end+1} = new_path;
    end
end
end


function adjacent = find_adjacent(graph,row,col,only_valid)
% assume valid is 1, invalid is 0
% only_valid (bool) indicate whether only will return the valid adjacent position
adjacent = {};
col_max = size(graph,1);
row_max = size(graph,2);
if (row<1 || row>row_max || col<1 || col>col_max)
    return;
end
if (row-1 >= 1)
    if only_valid
        if graph(row-1, col) == 1
            adjacent{end+1} = [row-1,col];
        end
    else
        adjacent{end+1} = [row-1,col];
    end
end

if (row+1 <= row_max)
    if only_valid
        if graph(row+1, col) == 1
            adjacent{end+1} = [row+1,col];
        end
    else
        adjacent{end+1} = [row+1,col];
    end
end

if (col-1 >= 1)
    if only_valid
        if graph(row,col-1) == 1
            adjacent{end+1} = [row,col-1];
        end
    else
        adjacent{end+1} = [row,col-1];
    end
end

if (col+1 <= col_max)
    if only_valid
        if graph(row,col+1) == 1
            adjacent{end+1} = [row,col+1];
        end
    else
        adjacent{end+1} = [row,col+1];
    end
end
end


function surrounding = find_surrounding(graph,row,col, only_valid)
surrounding = {};
col_max = size(graph,1);
row_max = size(graph,2);
if (row<1 || row>row_max || col<1 || col>col_max)
    return;
end

col_range = {col};
if (col-1 >= 1)
    col_range{end+1} = col-1;
end
if (col+1 <= col_max)
    col_range{end+1} = col+1;
end

row_range = {row};
if (row-1 >= 1)
    row_range{end+1} = row-1;
end
if (row+1 <= row_max)
    row_range{end+1} = row+1;
end

for r = row_range
    for c = col_range
        if ((r{:} == row) & (c{:} == col))
                continue
        end
        if only_valid
            if graph(r{:},c{:}) == 1
                surrounding{end+1} = [r,c];
            end
        else
            surrounding{end+1} = [r,c];
        end    
    end
end
end