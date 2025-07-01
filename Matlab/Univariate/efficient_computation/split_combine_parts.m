function [varargout] = split_combine_parts(dims, direction, varargin)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% No description available for function split_combine_parts.
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
	% direction:			integer (1, 2, 3, or 4)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Outputs:
	% 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	NumParts = length(dims.dims1);
	NumObs = dims.NumObs;
	NumParams = dims.NumParams;
	
	% beta_vec => beta_parts
	if direction == 1
		beta_vec = varargin{1};
		beta_parts = cell(NumParts, 1);
		idx = 0;
		for ii = 1:NumParts
			NumXi = dims.Xpart_2_NumX(ii); % TO DO
			beta_parts{ii}.beta = beta_vec(idx+1:idx+NumXi);
			idx = idx + NumXi;
			
			if isfield(dims, 'Xpart_2_NumX_FEs')
				NumX_FEs_i = dims.Xpart_2_NumX_FEs(ii);
				for ff = 1:NumX_FEs_i
					Num_FE_Vals_if = dims.Xpart_2_Num_FE_vals{ii}(ff);
					FE_vals = beta_vec(idx+1:idx+Num_FE_Vals_if-1);
					idx = idx + Num_FE_Vals_if-1;
					beta_parts{ii}.FE_vals{ff} = [0;FE_vals]; % Set first one equal to zero
				end			
			end
		end
		varargout{1} = beta_parts;
	end
	
	% Xbeta_parts => Xbeta
	if direction == 2
		Xbeta_parts = varargin{1};
		varargout{1} = combine_Xbeta_parts(dims, Xbeta_parts);
	end
	
	% Xt_Y_parts => Xt_Y
	if direction == 3
		Xt_Y_parts = varargin{1};
		Xt_Y = zeros(1, NumParams);
		idx = 0;
		for ii = 1:NumParts
			NumXi = dims.Xpart_2_NumX(ii);
			Xt_Y(idx+1:idx+NumXi) = Xt_Y_parts{ii}.X_Y;
			idx = idx + NumXi;
			
			
			if isfield(dims, 'Xpart_2_NumX_FEs')
				NumX_FEs_i = dims.Xpart_2_NumX_FEs(ii);
				for ff = 1:NumX_FEs_i
					Num_FE_Vals_if = dims.Xpart_2_Num_FE_vals{ii}(ff);
					myvals = Xt_Y_parts{ii}.X_FEs_Y{ff};
					myvals = myvals(2:end); % Remove first FE values
					Xt_Y(idx+1:idx+Num_FE_Vals_if-1) = myvals;
					idx = idx + Num_FE_Vals_if-1;				
				end
			end
		end
		varargout{1} = Xt_Y;
	end
	
	% Xt_X_Y_parts => Xt_X_Y
	if direction == 4
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
		Xt_X_Y_parts = varargin{1};
		Xt_X_Y = zeros(NumParams, NumParams);
		
		idx1 = 0;
		for ii = 1:NumParts
			NumXi            = dims.Xpart_2_NumX(ii);
			if isfield(dims, 'Xpart_2_NumX_FEs')
				NumX_FEs_i      = dims.Xpart_2_NumX_FEs(ii);
				NumFEvals2Keep_i = dims.NumFEvals2Keep{ii};
			else
				NumX_FEs_i = 0;
				NumFEvals2Keep_i = 0;
			end
			idces_Xi     = idx1+1:idx1+NumXi;
			idces_XFEs_i = idx1+NumXi+1:idx1+NumXi+NumFEvals2Keep_i;
			
			idx2 = 0;
			for jj = 1:NumParts
				NumXj            = dims.Xpart_2_NumX(jj);
				if isfield(dims, 'Xpart_2_NumX_FEs')
					NumX_FEs_j      = dims.Xpart_2_NumX_FEs(jj);
					NumFEvals2Keep_j = dims.NumFEvals2Keep{jj};
				else
					NumX_FEs_j = 0;
					NumFEvals2Keep_j = 0;
				end
				idces_Xj     = idx2+1:idx2+NumXj;
				idces_XFEs_j = idx2+NumXj+1:idx2+NumXj+NumFEvals2Keep_j;
				
				XY_X_ij      = Xt_X_Y_parts{ii,jj}.XY_X;          % NumXi x NumXj
				XY_XFE_ij    = Xt_X_Y_parts{ii,jj}.XY_XFE;
				XFE_YX_ij    = Xt_X_Y_parts{ii,jj}.XFE_YX;
				XFE_Y_XFE_ij = Xt_X_Y_parts{ii,jj}.XFE_Y_XFE;
				
				% Initialize structures
				tmp1 = Xt_X_Y_parts{ii,jj}.XY_X;               % NumXi x NumXj
				tmp2 = zeros(NumXi, NumFEvals2Keep_j);            % NumXi x NumFEvals2Keep_j
				tmp3 = zeros(NumFEvals2Keep_i, NumXj);            % NumFEvals2Keep_i x NumXj
				tmp4 = zeros(NumFEvals2Keep_i, NumFEvals2Keep_j); % NumFEvals2Keep_i x NumFEvals2Keep_j
				
				% Fill up tmp2
				a2 = 0;
				for f2 = 1:NumX_FEs_j
					Num_FE_Vals_jf = dims.Xpart_2_Num_FE_vals{jj}(f2);
					tmp2(:,a2+1:a2+Num_FE_Vals_jf-1) =  XY_XFE_ij{f2}(:,2:end);
					a2 = a2 + Num_FE_Vals_jf-1;
				end
				
				% Fill up tmp3
				a1 = 0;
				for f1 = 1:NumX_FEs_i
					Num_FE_Vals_if = dims.Xpart_2_Num_FE_vals{ii}(f1);
					tmp3(a1+1:a1+Num_FE_Vals_if-1,:) =  XFE_YX_ij{f1}(2:end,:);
					a1 = a1 + Num_FE_Vals_if-1;
				end
				
				% Fill up tmp4
				a1 = 0;
				for f1 = 1:NumX_FEs_i
					Num_FE_Vals_if = dims.Xpart_2_Num_FE_vals{ii}(f1);
					a2 = 0;
					for f2 = 1:NumX_FEs_j
						Num_FE_Vals_jf = dims.Xpart_2_Num_FE_vals{jj}(f2);
						tmp4(a1+1:a1+Num_FE_Vals_if-1, a2+1:a2+Num_FE_Vals_jf-1) = XFE_Y_XFE_ij{f1,f2}(2:end, 2:end);
						a2 = a2 + Num_FE_Vals_jf-1;
					end
					a1 = a1 + Num_FE_Vals_if-1;
				end
				
				% Store everything in the right place
				Xt_X_Y(idces_Xi, idces_Xj)         = tmp1;
				Xt_X_Y(idces_Xi, idces_XFEs_j)     = tmp2;
				Xt_X_Y(idces_XFEs_i, idces_Xj)     = tmp3;
				Xt_X_Y(idces_XFEs_i, idces_XFEs_j) = tmp4;
				
				% Update idx2
				idx2 = idx2 + NumXj + NumFEvals2Keep_j;
			end
			% Update idx1
			idx1 = idx1 + NumXi + NumFEvals2Keep_i;
		end
		varargout{1} = Xt_X_Y; % NumParams x NumParams
	end
end
