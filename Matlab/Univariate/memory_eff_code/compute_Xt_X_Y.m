function [Xt_X_Y, Y] = compute_Xt_X_Y(data, Y)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% This function computes sum_{n} Y_{n} * X_{na} * X_{nb} for all (a,b).
	% This function returns the result, combining everything together.
	%
	% I pass Y as an output argument for memory optimization reason:
	% - I need to "update" Y by reshaping it
	% - by doing so, Matlab would normally make a local copy of Y (which is memory heavy + takes time)
	% - by letting Y be an output variable (and in the calling function, "update" Y accordingly),
	%		I tell Matlab that it's okay to directly edit the same copy of Y (no need to make a local copy)
	% - for this to work, I also need to "update" Y when calling this function
	% - I only reshape Y back and forth here, so Y is identical at the start and at the end, thus it is safe.
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Inputs:
	% dims:					object:
	%   . dimType:				0, 1, 2 or 3
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
	% Xt_X_Y:			NumParams x NumParams
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	if ~isfield(data, 'Xparts')
		Xt_X_Y = data.X' * (data.X .* Y);
		return;
	end
	
	% Unbundle things
	dims   = data.dims;
	Xparts = data.Xparts;
	
	if dims.dimType == 0 % general
		[Xt_X_Y, Y] = compute_Xt_X_Y_general(dims, Xparts, Y);
	end
	
	if dims.dimType == 1 % dim1
		[Xt_X_Y, Y] = compute_Xt_X_Y_dim1(dims, Xparts, Y);
	end
	
	if dims.dimType == 2 % dim2
		[Xt_X_Y, Y] = compute_Xt_X_Y_dim2(dims, Xparts, Y);
	end
	
	if dims.dimType == 3 % dim3
		[Xt_X_Y, Y] = compute_Xt_X_Y_dim3(dims, Xparts, Y);
	end
end
