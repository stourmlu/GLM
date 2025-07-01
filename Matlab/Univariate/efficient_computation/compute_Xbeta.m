function [Xbeta] = compute_Xbeta(data, beta_vec)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
	%   {ii}:					object
	%		.X:						dim1_ii x NumX
	%		.X_FEs:					dim1_ii x NumX_FEs
	% beta_parts:			cell(NumParts,1)
	%	{kk}:					object
	%		.beta:					NumX x 1
	%		.FE_vals:				cell(NumX_FEs, 1)
	%			{ii}:					NumFE_vals_ii x 1
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Outputs:
	% Xbeta:				NumObs x 1
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	if ~isfield(data, 'Xparts')
		Xbeta = data.X * beta_vec;
		return;
	end
	
	% Unbundle things
	dims   = data.dims;
	Xparts = data.Xparts;
	beta_parts = split_combine_parts(dims, 1, beta_vec);
	
	% Compute each part
	NumParts = length(Xparts);
	Xbeta_parts = cell(NumParts, 1);
	for ii = 1:NumParts
		Xbeta_parts{ii} = compute_Xbeta_part(Xparts{ii}, beta_parts{ii});
	end
	
	% Put all parts together
	Xbeta = split_combine_parts(dims, 2, Xbeta_parts); % NumObs x 1
end
