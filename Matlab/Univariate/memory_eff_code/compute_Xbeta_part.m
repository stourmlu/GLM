function [res] = compute_Xbeta_part(Xpart, beta_part)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% No description available for function compute_Xbeta_part.
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Inputs:
	% Xpart:				object:
	%	.X:						dim1 x NumX
	%	.X_FEs:					dim1 x NumX_FEs
	% beta_part:			object:
	%	.beta:					NumX x 1
	%	.FE_vals:				cell(NumX_FEs, 1)
	%		{ii}:					NumFE_vals_ii x 1
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Outputs:
	% res:					dim1 x 1
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	if length(beta_part.beta) == 0
		res = zeros(size(Xpart.X,1),1);
	else
		res = Xpart.X * beta_part.beta;
	end
	if isfield(Xpart, 'X_FEs')
		for ii = 1:size(Xpart.X_FEs,2)
			FE_vals_ii = beta_part.FE_vals{ii};
			res = res + FE_vals_ii(Xpart.X_FEs(:,ii));
		end
	end
end
