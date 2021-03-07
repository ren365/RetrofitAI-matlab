X = zeros(4,4);
adja = find_adjacent(X, 2, 2, false);
%disp(size(adja));

graph = [1,1,1,1; 1,0,1,0; 1,1,1,1; 1,1,1,1];

tic
all_path = bfs_search(graph, [1,1], [4,4])
disp("excute time for the bfs search is: ")
toc
for path = all_path
    disp("A new path");
    celldisp(path)
end
toc


function all_path = bfs_search(graph, start, target)
% graph: N x M Matrix
% start: [row, col]: 1x2 matrix contain the position of start position
% end: [row, col] 1x2 matrix contain the position of target position
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
    disp("last node is ");
    disp(last_node);
    adjacent = find_adjacent(visited,last_node(1),last_node(2),true);
    
    for adj=adjacent
        visited(adj{1}(1),adj{1}(2)) = 0;
        new_path = path;
        new_path{end+1} = [adj{1}(1), adj{1}(2)];
        queue{end+1} = new_path;
    end
    disp("intermediate path");
    for pa=path
        celldisp(pa);
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