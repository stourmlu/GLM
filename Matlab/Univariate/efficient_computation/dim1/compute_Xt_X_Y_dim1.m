function [Xt_X_Y, Y] = compute_Xt_X_Y2_dim1(dims, Xparts, Y)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% This function computes sum_{n} Y_{n} * X_{na} * X_{nb} for all (a,b).
	% This function returns the result, combining everything together.
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Inputs:
	% dims:					object:
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
	% Xt_X_Y:			NumParams x NumParams
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
% This part assumes that the data is made of 1 part.
% Full data: T x 1
% Xpart 1:   T x 1
	
	obj = compute_Xt_X_Y_samePart(Xparts{1}, Y);
	Xt_X_Y_parts = {obj};
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Xt_X_Y_parts: 			cell(NumParts, NumParts)
	%	{ii,jj}:					object:
	%		.XY_X:						NumXi x NumXj
	%		.XY_XFE:					cell(NumX_FEs_j,1):
	%			{f2}:						NumXi x NumX_FE_vals_j{f2}
	%		.XFE_YX:					cell(NumX_FEs_i,1):
	%			{f1}:						NumX_FE_vals_i{f1} x NumXj
	%		.XFE_Y_XFE:					cell(NumX_FEs_i,NumX_FEs_j):
	%			{f1,f2}:					NumX_FE_vals_i{f1} x NumX_FE_vals_j{f2}
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	% Put all parts together
	Xt_X_Y = split_combine_parts(dims, 4, Xt_X_Y_parts); % 1 x NumParams
end
