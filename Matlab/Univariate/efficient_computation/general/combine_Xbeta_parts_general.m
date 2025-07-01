function [Xbeta] = combine_Xbeta_parts_general(dims, Xbeta_parts)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% No description available for function combine_Xbeta_parts_general.
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Inputs:
	% dims:					object:
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
	
% This part does not assume that NumObs contains the full cartesian product of different dimensions in any kind of way.
% However, it can be extremely slow! Expanding Xbeta_parts{ii} using dims.mappings{ii} can be very expensive.

	NumParts = length(dims.mappings);
	NumObs = dims.NumObs;
	Xbeta = zeros(NumObs,1);
	for ii = 1:NumParts
		Xbeta = Xbeta + Xbeta_parts{ii}(dims.mappings{ii});
	end
end
