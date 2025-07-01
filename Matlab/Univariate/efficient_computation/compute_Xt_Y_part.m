function [res] = compute_Xt_Y_part(Xpart, Yval)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% No description available for function compute_Xt_Y_part.
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Inputs:
	% Xpart:				object:
	%	.X:						dim1 x NumX
	%	.X_FEs:					dim1 x NumX_FEs
	%	.NumX_FE_vals:			1 x NumX_FEs: gives integer
	% Yval:					dim1 x 1
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Outputs:
	% res:					object:
	%	.X_Y:					1 x NumX
	%	.X_FEs_Y:				cell(NumX_FEs,1)
	%		{ii}:					1 x NumFE_vals_ii
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	% Compute part with respect to X
	res.X_Y = Yval' * Xpart.X; % 1 x NumX
	
	% Compute part with respect to X_FEs
	if isfield(Xpart, 'X_FEs')
		NumFEs = size(Xpart.X_FEs,2);
		X_FEs_Y = cell(NumFEs,1);
		for ii = 1:NumFEs
			NumFE_vals_ii = Xpart.NumX_FE_vals(ii);
			X_FEs_Y{ii} = accumarray(Xpart.X_FEs(:,ii), Yval, [NumFE_vals_ii 1])'; % 1 x NumFE_vals_ii
		end
		res.X_FEs_Y = X_FEs_Y;
	else
		res.X_FEs_Y = cell(0,1);
	end
end
