function [Xbeta] = combine_Xbeta_parts_dim3(dims, Xbeta_parts)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% No description available for function combine_Xbeta_parts_dim3.
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
% Full data: (T*K1*K2) x 1
% Xpart 1:   (T*K1*K2)  x 1
% Xpart 2:   (T*K1)    x 1
% Xpart 3:   (T*K2)    x 1
% Xpart 4:   (K1*K2)   x 1
% Xpart 5:   T         x 1
% Xpart 6:   K1        x 1
% Xpart 7:   K2        x 1

	T  = dims.dims1{5};
	K1 = dims.dims1{6};
	K2 = dims.dims1{7};
	Xbeta_parts{1} = reshape(Xbeta_parts{1}, [T K1 K2]);
	Xbeta_parts{2} = reshape(Xbeta_parts{2}, [T K1  1]);
	Xbeta_parts{3} = reshape(Xbeta_parts{3}, [T  1 K2]);
	Xbeta_parts{4} = reshape(Xbeta_parts{4}, [1 K1 K2]);
	Xbeta_parts{5} = reshape(Xbeta_parts{5}, [T  1  1]);
	Xbeta_parts{6} = reshape(Xbeta_parts{6}, [1 K1  1]);
	Xbeta_parts{7} = reshape(Xbeta_parts{7}, [1  1 K2]);
	Xbeta = Xbeta_parts{1} + Xbeta_parts{2} + Xbeta_parts{3} + Xbeta_parts{4} + Xbeta_parts{5} + Xbeta_parts{6} + Xbeta_parts{7};
	Xbeta = reshape(Xbeta, [T*K1*K2 1]);
end
