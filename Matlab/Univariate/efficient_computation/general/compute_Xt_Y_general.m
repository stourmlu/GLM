function [Xt_Y, Y] = compute_Xt_Y_general(dims, Xparts, Y)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% This function computes sum_{n} Y_{n} * X_{na} for all a.
	% This function returns the result, combining everything together.
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Inputs:
	% dims:					object:
	%	.mappings:				cell(NumParts,1)
	%		{ii}:					NumObs x 1 (gives integer values between 1 and dim1_ii)
	%	.dims1:					cell(NumParts,1)
	%		{ii}:					integer: gives dim1_ii that corresponds to Xparts{ii}
	%	.Xpart_2_NumX:			1 x NumParts: gives integers (dim2 of Xparts{ii}.X)
	%	.Xpart_2_NumX_FEs:		1 x NumParts: gives integers (dim2 of Xparts{ii}.X_FEs)
	%	.Xpart_2_Num_FE_vals:	cell(1,NumParts)
	%		{ii}:					1 x NumX_FEs_i  --> gives integer: number of possible values for Xparts{ii}.X_FEs(:,ff)
	%	.NumFEvals2Keep:		cell(1,NumParts) --> gives integer
	%	.NumParts
	%	.NumObs
	%	.NumParams
	% Xparts:				cell(NumParts,1)
	%	{ii}:					object
	%		.X:						dim1_ii x NumX
	%		.X_FEs:					dim1_ii x NumX_FEs
	%		.NumX_FE_vals:			1 x NumX_FEs: gives integer
	% Y:					NumObs x 1
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Outputs:
	% Xt_Y:				1 x NumParams
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
% This part does not assume that NumObs contains the full cartesian product of different dimensions in any kind of way.
% However, it is slow!
	
	% Compute each part
	NumParts = length(Xparts);
	Xt_Y_parts = cell(NumParts,1);
	for ii = 1:NumParts
		myvec = dims.mappings{ii}; % NumObs x 1 (gives integer values between 1 and dim1_ii)
		dim1_ii = dims.dims1{ii};
		Ypart_ii = accumarray(myvec, Y, [dim1_ii 1]); % dim1_ii x 1
		Xt_Y_parts{ii} = compute_Xt_Y_part(Xparts{ii}, Ypart_ii);
	end
	
	% Put all parts together
	Xt_Y = split_combine_parts(dims, 3, Xt_Y_parts); % 1 x NumParams
end
