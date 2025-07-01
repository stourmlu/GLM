function [Xbeta] = combine_Xbeta_parts_dim2(dims, Xbeta_parts)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% No description available for function combine_Xbeta_parts_dim2.
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
	
% This part assumes that the data is made of 6 parts and that the data is "complete" in the sense that we have the
% ful cartesian product.
% Full data: (T*K) x 1
% Xpart 1:   (T*K) x 1
% Xpart 2:     T   x 1
% Xpart 3:     K   x 1

	T  = dims.dims1{2};
	K = dims.dims1{3};
	
	Xbeta_parts{1} = reshape(Xbeta_parts{1}, [T K]);
	Xbeta_parts{2} = reshape(Xbeta_parts{2}, [T 1]);
	Xbeta_parts{3} = reshape(Xbeta_parts{3}, [1 K]);
	
	Xbeta = Xbeta_parts{1} + Xbeta_parts{2} + Xbeta_parts{3};
	Xbeta = reshape(Xbeta, [T*K 1]);
end
