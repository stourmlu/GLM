function [beta_star] = estimate_OLS(data)
	
	g = compute_Xt_Y(data, data.Y)';  % NumX x 1
	
	if isfield(data, 'Xparts')
		NumObs = data.dims.NumObs;
	else
		NumObs = size(data.X,1);
	end
	H = compute_Xt_X_Y(data, ones(NumObs,1)); % NumX x NumX
	
	beta_star = H\g;
end
