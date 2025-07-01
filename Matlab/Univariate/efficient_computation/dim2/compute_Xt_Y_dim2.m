function [Xt_Y, Y] = compute_Xt_Y_dim2(dims, Xparts, Y)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% This function computes sum_{n} Y_{n} * X_{na} for all a.
	% This function returns the result, combining everything together.
	
	% I pass Y as an output argument for memory optimization reason:
	% - I need to "update" Y by reshaping it
	% - by doing so, Matlab would normally make a local copy of Y (which is memory heavy + takes time)
	% - by letting Y be an output variable (and in the calling function, "update" Y accordingly),
	%		I tell Matlab that it's okay to directly edit the same copy of Y (no need to make a local copy)
	% - for this to work, I also need to "update" Y when calling this function
	% - I only reshape Y back and forth here, so Y is identical at the start and at the end, so it is safe.
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
	% Xt_Y:				1 x NumParams
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function assumes that the data is made of 3 parts and that the data is "complete" in the sense that we have the
% ful cartesian product.
% Full data: (T*K) x 1
% Xpart 1:   (T*K) x 1
% Xpart 2:     T   x 1
% Xpart 3:     K   x 1
	
	% Compute each part
	NumParts = length(Xparts);
	Xt_Y_parts = cell(NumParts,1);

	T = dims.dims1{2};
	K = dims.dims1{3};
	
	Y          = reshape(Y,[T K]);   % T x K
	Ysum_t     = sum(Y,2);           % T x 1
	Ysum_k     = sum(Y,1)';          % K x 1
	Y          = reshape(Y,[T*K 1]); % (T*K) x 1
	
	Xt_Y_parts{1} = compute_Xt_Y_part(Xparts{1}, Y);
	Xt_Y_parts{2} = compute_Xt_Y_part(Xparts{2}, Ysum_t);
	Xt_Y_parts{3} = compute_Xt_Y_part(Xparts{3}, Ysum_k);
	
	% Put all parts together
	Xt_Y = split_combine_parts(dims, 3, Xt_Y_parts); % 1 x NumParams
end
