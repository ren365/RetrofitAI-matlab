function result = nearest_distance(barriers,current_position,radius,barriers_numbers)
	dist = [];
	for j=1:barriers_numbers
		xtmp = barriers(j*2-1);
		ytmp = barriers(j*2);
		dis = sqrt((xtmp-current_position(1))^2+(ytmp-current_position(2))^2) - radius(j);
		dist= [dist,dis];
	end
	result = min(dist);
end