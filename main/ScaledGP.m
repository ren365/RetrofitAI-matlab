classdef  ScaledGP
	% untest
	properties
		xdim
		ydim
		ystd
		ymean
		xstd
		xmean
		m
	end
	methods (Access = public)
		function obj = ScaledGP(xdim,ydim)
			obj.xdim = xdim;
			obj.ydim = ydim;
			obj.ystd = ones(1,ydim);
			obj.ymean = zeros(1,ydim);
			obj.xstd = ones(1,xdim);
			obj.xmean = zeros(1,xdim);
			% obj.m = [fitrgp(zeros(1,xdim),zeros(1,ydim));
			obj.m = {fitrgp(zeros(1,xdim),zeros(1,1)),fitrgp(zeros(1,xdim),zeros(1,1))};
		end
		
		function obj = optimize(obj,x,y)

			[xmean,xstd] = obj.update_xscale(x);
			[ymean,ystd] = obj.update_yscale(y);

			x = obj.scalex(x,xmean,xstd);
			y = obj.scaley(y,ymean,ystd);
			
			% obj.m = fitrgp(x,y(:,1),'Optimizer','lbfgs')
			obj.m = {fitrgp(x,y(:,1)),fitrgp(x,y(:,2))};
			
			obj.xmean = xmean;
			obj.xstd = xstd;
			obj.ymean = ymean;
			obj.ystd = ystd;
			
		end
				
		function [mean0,var] = predict(obj,x)
			x = obj.scalex(x',obj.xmean,obj.xstd);
			% [mean,var] = predict(obj.m,x);
			
			[mean1,var1] = predict(obj.m{1},x);
			[mean2,var2] = predict(obj.m{2},x);
			mean0 = [mean1,mean2];
			% simplified version for std calculation, due to no multi-gp function in matlab
			% the origianl value of "var" is between 10^-7 to 10^-11
			var = var1*var2;
			mean0 = obj.unscaley(mean0,obj.ymean,obj.ystd);
			var = var * obj.ystd;
		end

		function [xmean,xstd] = update_xscale(obj,x)
			xmean = mean(x);
			xstd = std(x);
		end

		function [ymean,ystd] = update_yscale(obj,y)
			ymean = mean(y);
			ystd = std(y);
		end
		
		function result = scalex(obj,x,xmean,xstd)
			if any((xstd == 0))
				result = x-xmean;
			else
				result = (x - xmean) ./ xstd;
			end
		end
		
		function result = scaley(obj,y,ymean,ystd)
			if any((ystd == 0))
				result = y-ymean;
			else
				result = (y - ymean) ./ ystd;
			end
		end
		
		function result = unscalex(obj,x,xmean,xstd)
			if any((xstd == 0))
				result = x + xmean;
			else
				result = x .* xstd + xmean;
			end
		end
		
		function result = unscaley(obj,y,ymean,ystd)
			if any((ystd == 0))
				result = y + ymean;
			else
				result = y .* ystd + ymean;
			end
		end

	end
end
