function x_d = FSM(start_position,end_position,barrier_x,barrier_y)
    map = create_map(barrier_x, barrier_y, start_position,end_position);

    input = sense([1,1],map);


    curr_position = start_position;
    curr_value = get_value_by_position(curr_position,map);



    counter = 1;

    trace_x(counter) = start_position(1);
    trace_y(counter) = start_position(2);

    %{
    for i = 1:11
        input = sense_d(curr_position,map);
        next_position = get_next_state(input, curr_position);
        map(6-next_position(2),1+next_position(1)) = map(6-next_position(2),1+next_position(1)) + 1;
        curr_position = next_position;
        curr_value = get_value_by_position(curr_position,map);
    end
    %}



    tic
    while(curr_value ~= log(0))
        input = sense(curr_position,map);
        next_position = get_next_state(input, curr_position);
        counter = counter+1;
        trace_x(counter) = next_position(1);
        trace_y(counter) = next_position(2);
        map(6-next_position(2),1+next_position(1)) = map(6-next_position(2),1+next_position(1)) + 1;
        curr_position = next_position;
        curr_value = get_value_by_position(curr_position,map);
    end
    toc

    discretization = 5;% every step is divided intro 2 steps
    %T = length(trace_x)-1;
    %t = linspace(1,length(trace_x),discretization*T*2);

    trace_2x(1) = start_position(1);
    trace_2y(1) = start_position(2);


    for i = 1:(length(trace_x)-1)*2
        if mod(i,2) == 1
            trace_2x(i) = trace_x((i+1)/2);
            trace_2y(i) = trace_y((i+1)/2);
        else
            trace_2x(i) = (trace_x(i/2) + trace_x(i/2+1))/2.0;
            trace_2y(i) = (trace_y(i/2) + trace_y(i/2+1))/2.0;
        end
    end

    trace_2x(2*length(trace_x)-1) = trace_x(length(trace_x));
    trace_2y(2*length(trace_y)-1) = trace_y(length(trace_y));


    dt = 1/2/discretization;
    %construct x_d
    xdx(1) = start_position(1);
    xdy(1) = start_position(2);
    %theta = 0;
    %velo 


itr = 1;
while itr <= length(trace_2x) - 2
    curr_vector_x = trace_2x(itr+1) - trace_2x(itr);
    curr_vector_y = trace_2y(itr+1) - trace_2y(itr);
    
    next_vector_x = trace_2x(itr+2) - trace_2x(itr+1);
    next_vector_y = trace_2y(itr+2) - trace_2y(itr+1);
    

    turn = curr_vector_x * next_vector_x + curr_vector_y * next_vector_y;

    if turn ~= 0
        for j = 1:discretization
            if curr_vector_x > 0 
                xdx(end+1) = trace_2x(itr) + dt*j;
            end
            
            if curr_vector_x == 0
                xdx(end+1) = trace_2x(itr);
            end
            
            if curr_vector_y > 0
                xdy(end+1) = trace_2y(itr) + dt*j;
            end
            
            if curr_vector_y < 0
                xdy(end+1) = trace_2y(itr) - dt*j;
            end
            
            if curr_vector_y == 0
                if (itr == 1)
                    xdy(end+1) = start_position(2);
                else 
                    xdy(end+1) = trace_2y(itr);
                end
            end
        end
        itr = itr+1;
    else
        for j = 1:discretization
            if curr_vector_x > 0 && next_vector_y > 0
                xdx(end+1) = trace_2x(itr) + dt*j;
                xdy(end+1) = trace_2y(itr) + 0.5 - sqrt(0.5*0.5 - j*dt * j*dt);
            end
            if curr_vector_x > 0 && next_vector_y < 0
                xdx(end+1) = trace_2x(itr) + dt*j;
                xdy(end+1) = trace_2y(itr) - (0.5 - sqrt(0.5*0.5 - j*dt * j*dt));
                
                
            end
            if curr_vector_y > 0 
                xdx(end+1) = trace_2x(itr) + 0.5 - sqrt(0.5*0.5 - j*dt * j*dt);
                xdy(end+1) = trace_2y(itr) + dt*j;
            end
            if curr_vector_y < 0
                xdx(end+1) = trace_2x(itr) + (0.5 - sqrt(0.5*0.5 - j*dt * j*dt));
                xdy(end+1) = trace_2y(itr) - dt*j;
            end
            
            
        end
        itr = itr+2;
    end
    
    
            
    
end

x_d = [xdx; xdy; zeros(1,length(xdx)); zeros(1,length(xdx))];
x_d(3,1:end-1) = atan2(diff(x_d(2,:)),diff(x_d(1,:)));
x_d(4,1:end-1) = sqrt(diff(x_d(1,:)).^2 + diff(x_d(2,:)).^2)/dt;
x_d(3,end)=x_d(3,end-1);
x_d(4,end)=x_d(4,end-1);


%{
plot(xdx,xdy,'o')
hold on
plot(trace_2x,trace_2y)
ylim([0,5]);
yticks(0:0.2:5);
xlim([0,10]);
xticks(0:0.2:10);
%}


function map = create_map(bx, by, start_pos, end_pos)
    map = zeros(11,11);
    for i = 1:length(bx)
        bx_map = 6-by(i);
        
        by_map = bx(i)+1;
       
        map(bx_map,by_map) = exp(1000);
        %surroundings should also be inf
        %side
        if (bx_map) == 1
            map(bx_map+1,by_map) = exp(1000);
            map(bx_map,by_map-1) = exp(1000);
            map(bx_map,by_map+1) = exp(1000);
        elseif bx_map == length(map(1,:))
            map(bx_map-1,by_map) = exp(1000);
            map(bx_map,by_map-1) = exp(1000);
            map(bx_map,by_map+1) = exp(1000);
        else
            map(bx_map+1,by_map) = exp(1000);
            map(bx_map-1,by_map) = exp(1000);
            map(bx_map,by_map-1) = exp(1000);
            map(bx_map,by_map+1) = exp(1000);
        end
        
        
    end
    map(6-start_pos(2),start_pos(1)+1) = 1;
    map(6-end_pos(2), end_pos(1)+1) = log(0);
end

function input_value = sense(curr_pos,map)
    x = 6 - curr_pos(2);
    y = curr_pos(1) + 1;
    %[up right down left]
    %corner
    if (x == 1 && y == 1)
        input_value = [exp(1000), map(x,y+1),map(x+1,y), exp(1000)];
    elseif (x == 1 && y == length(map(1,:)))
        input_value = [exp(1000),exp(1000),map(x+1,y),map(x,y-1)];
    elseif (x == length(map(1,:)) && y == 1)
        input_value = [map(x-1,y), map(x,y+1), exp(1000), exp(1000)];
    elseif (x == length(map(1,:)) && y == length(map(1,:)))
        input_value = [map(x-1,y), exp(1000), exp(1000), map(x,y-1)];
    %side
    elseif (x == 1 && y ~= 1 && y ~= length(map(1,:)))
        input_value = [exp(1000), map(x,y+1),map(x+1,y), map(x,y-1)];
    elseif (x == length(map(1,:)) && y ~= 1 && y ~= length(map(1,:)))
        input_value = [map(x-1,y), map(x,y+1), exp(1000), map(x,y-1)];
    elseif (x~=1 && x ~= length(map(1,:)) && y == 1)
        input_value = [map(x-1,y), map(x,y+1),map(x+1,y), exp(1000)];
    elseif (x~=1 && x ~= length(map(1,:)) && y == length(map(1,:)))
        input_value = [map(x-1,y),exp(1000),map(x+1,y), map(x,y-1)];
    %default
    else
        input_value = [map(x-1,y), map(x,y+1), map(x+1,y) map(x,y-1)];   
    end
end

function input_value = sense_d(curr_pos,map)
    x = 6 - curr_pos(2);
    y = curr_pos(1) + 1;
    %[up right down left]
    %corner
    if (x == 1 && y == 1)
        input_value = [exp(1000), d(x,y+1,map),d(x+1,y,map), exp(1000)];
    elseif (x == 1 && y == length(map(1,:)))
        input_value = [exp(1000),exp(1000),d(x+1,y,map),d(x,y-1,map)];
    elseif (x == length(map(1,:)) && y == 1)
        input_value = [d(x-1,y,map), d(x,y+1,map), exp(1000), exp(1000)];
    elseif (x == length(map(1,:)) && y == length(map(1,:)))
        input_value = [d(x-1,y,map), exp(1000), exp(1000), d(x,y-1,map)];
    %side
    elseif (x == 1 && y ~= 1 && y ~= length(map(1,:)))
        input_value = [exp(1000), d(x,y+1,map), d(x+1,y,map), d(x,y-1,map)];
    elseif (x == length(map(1,:)) && y ~= 1 && y ~= length(map(1,:)))
        input_value = [d(x-1,y,map), d(x,y+1,map), exp(1000), d(x,y-1,map)];
    elseif (x~=1 && x ~= length(map(1,:)) && y == 1)
        input_value = [d(x-1,y,map), d(x,y+1,map), d(x+1,y,map), exp(1000)];
    elseif (x~=1 && x ~= length(map(1,:)) && y == length(map(1,:)))
        input_value = [d(x-1,y,map), exp(1000), d(x+1,y,map), d(x,y-1,map)];
    %default
    else
        input_value = [d(x-1,y,map), d(x,y+1,map), d(x+1,y,map), d(x,y-1,map)];   
    end
end


function nxt_pos = get_next_state(input_value, curr_pos)
    [v, i] = min(input_value);
    switch i
        case 1
            nxt_pos = [curr_pos(1), curr_pos(2)+1];
        case 2
            nxt_pos = [curr_pos(1)+1, curr_pos(2)];
        case 3
            nxt_pos = [curr_pos(1), curr_pos(2)-1];
        otherwise
            nxt_pos = [curr_pos(1)-1, curr_pos(2)];
    end
end

function value = get_value_by_position(pos,map)
    value = map(6 - pos(2),1+pos(1));
end

function distance = d(x,y,map)
    end_position = [10,0];
    end_x = 6 - end_position(2);
    end_y = end_position(1) + 1;
    

    if map(x,y) == exp(1000)
        distance = exp(1000);
    else
        distance = abs(x - end_x) + abs(y - end_y);
    end
end
end

 