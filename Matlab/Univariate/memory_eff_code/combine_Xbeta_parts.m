function [Xbeta] = combine_Xbeta_parts(dims, Xbeta_parts)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% No description available for function combine_Xbeta_parts.
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Inputs:
	% dims:					object:
	%   . dimType:				0, 1, 2 or 3
	%	.NumParts
	%	.NumObs
	%	.NumParams
	%	.mappings:				cell(NumParts,1)
	%		{ii}:					NumObs x 1 (gives integer values between 1 and dim1_ii)
	%	.dims1:					cell(NumParts,1)
	%		{ii}:					integer: gives dim1_ii that corresponds to Xparts{ii}
	%	.Xpart_2_NumX:			1 x NumParts: gives integers (dim2 of Xparts{ii}.X)
	%	.Xpart_2_NumX_FEs:		1 x NumParts: gives integers (dim2 of Xparts{ii}.X_FEs)
	%	.Xpart_2_Num_FE_vals:	cell(1,NumParts)
	%		{ii}:					1 x NumX_FEs_i  --> gives integer: number of possible values for Xparts{ii}.X_FEs(:,ff)
	%	.NumFEvals2Keep:		cell(1,NumParts) --> gives integer
	% direction:				integer (1, 2, 3, or 4)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Outputs:
	% Xbeta:					NumObs x 1
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	if dims.dimType == 0 % general
		Xbeta = combine_Xbeta_parts_general(dims, Xbeta_parts);
	end
	
	if dims.dimType == 1 % dim1
		Xbeta = combine_Xbeta_parts_dim1(dims, Xbeta_parts);
	end
	
	if dims.dimType == 2 % dim2
		Xbeta = combine_Xbeta_parts_dim2(dims, Xbeta_parts);
	end
	
	if dims.dimType == 3 % dim3
		Xbeta = combine_Xbeta_parts_dim3(dims, Xbeta_parts);
	end
end
