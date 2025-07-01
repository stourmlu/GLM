function [Xt_X_Y, Y] = compute_Xt_X_Y2_dim2(dims, Xparts, Y)
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
	% Xt_X_Y:			NumParams x NumParams
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
% This part assumes that the data is made of 6 parts and that the data is "complete" in the sense that we have the
% ful cartesian product.
% Full data: (T*K) x 1
% Xpart 1:   (T*K) x 1
% Xpart 2:     T   x 1
% Xpart 3:     K   x 1
	
	
	NumParts = dims.NumParts;
	T = dims.dims1{2};
	K = dims.dims1{3};
	
	% Reshape Y and compute aggregate versions
	Y          = reshape(Y,[T K]);   % T x K
	Ysum_t     = sum(Y,2);           % T x 1
	Ysum_k     = sum(Y,1)';          % K x 1
	Y          = reshape(Y,[T*K 1]); % (T*K) x 1
	
	
	Xt_X_Y_parts = cell(NumParts, NumParts);
	for ii = 1:NumParts
		Xi            = Xparts{ii}.X;
		NumXi         = dims.Xpart_2_NumX(ii); % TO DO
		if isfield(Xparts{ii}, 'X_FEs')
			X_FEs_i       = Xparts{ii}.X_FEs;
			NumX_FEs_i    = dims.Xpart_2_NumX_FEs(ii); % TO DO
			Num_FE_vals_i = dims.Xpart_2_Num_FE_vals{ii};
		else
			X_FEs_i       = zeros(dims.dims1{ii}, 0);
			NumX_FEs_i    = 0;
			Num_FE_vals_i = {};
		end
		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if ii == 1
			Ypart_ii = Y;      % (T*K) x 1
		end
		if ii == 2
			Ypart_ii = Ysum_t; % T x 1
		end
		if ii == 3
			Ypart_ii = Ysum_k; % K x 1
		end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		XY_i = Ypart_ii .* Xi; % dim1_i x NumXi
		
		for jj = ii:NumParts
			%%%% Case where ii == jj
			if ii == jj
				obj = compute_Xt_X_Y_samePart(Xparts{ii}, Ypart_ii);
				
				% Store for {ii,ii}
				Xt_X_Y_parts{ii,ii} = obj;
			
			%%%% Case where ii != jj
			else
				Xj            = Xparts{jj}.X;
				NumXj         = dims.Xpart_2_NumX(jj); % TO DO
				if isfield(Xparts{ii}, 'X_FEs')
					X_FEs_j       = Xparts{jj}.X_FEs;
					NumX_FEs_j    = dims.Xpart_2_NumX_FEs(jj); % TO DO
					Num_FE_vals_j = dims.Xpart_2_Num_FE_vals{jj}; % TO DO
				else
					X_FEs_j       = zeros(dims.dims1{jj}, 0);
					NumX_FEs_j    = 0;
					Num_FE_vals_j = {};
				end
				
				% Initialize structures XY_XFE, XFE_YX and XFE_Y_XFE
				XY_XFE    = cell(NumX_FEs_j,1);
				XFE_YX    = cell(NumX_FEs_i,1);
				XFE_Y_XFE = cell(NumX_FEs_i,NumX_FEs_j);
				
				
				if ii == 1 && jj == 2 % (T*K) and T
					% Compute XY_X
					XY = sum(reshape(XY_i, [T K NumXi]), 2); % T x 1 x NumXi
					XY = reshape(XY, [T NumXi]);               % T x NumXi
					XY_X = XY' * Xj;                           % NumXi x NumXj
					
					% Compute XY_XFE
					for f2 = 1:NumX_FEs_j
						NumX_FE_vals_jf = Num_FE_vals_j(f2);
						X_FE_j = Xparts{jj}.X_FEs(:,f2); % T x 1
						for aa = 1:NumXi
							XY_XFE{f2}(aa,:) = accumarray(X_FE_j, XY(:,aa), [NumX_FE_vals_jf 1])';
						end
					end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%					
					% Compute XFE_YX
					if NumX_FEs_i > 0; error('Current code does not support FEs at tt-kk level'); end; % TO DO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					% Compute XFE_Y_XFE
					%% Would require NumX_FEs_i > 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				end
				
				if ii == 1 && jj == 3 % (T*K) and K
					% Compute XY_X
					XY = sum(reshape(XY_i, [T K NumXi]), 1); % 1 x K x NumXi
					XY = reshape(XY, [K NumXi]);               % K x NumXi
					XY_X = XY' * Xj;                           % NumXi x NumXj
					
					% Compute XY_XFE
					for f2 = 1:NumX_FEs_j
						NumX_FE_vals_jf = Num_FE_vals_j(f2);
						X_FE_j = Xparts{jj}.X_FEs(:,f2); % K x 1
						for aa = 1:NumXi
							XY_XFE{f2}(aa,:) = accumarray(X_FE_j, XY(:,aa), [NumX_FE_vals_jf 1])';
						end
					end
					
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					% Compute XFE_YX
					if NumX_FEs_i > 0; error('Current code does not support FEs at tt-kk level'); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					% Compute XFE_Y_XFE
					%% Would require NumX_FEs_i > 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				end
				
				%%% Both have one dimension, they are different
				if ii == 2 && jj == 3 % T and K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					% Compute XY_X
					XY = reshape(Y, [T K]) * Xj; % T x NumXj
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					XY_X = Xi' * XY; % NumXi x NumXj
					XY2 = Xi' * reshape(Y, [T K]); % NumXi x K
%					XY_X2 = XY2 * Xj; % NumXi x NumXj
					
					% Compute XY_XFE
					for f2 = 1:NumX_FEs_j
						NumX_FE_vals_jf = Num_FE_vals_j(f2);
						X_FE_j = Xparts{jj}.X_FEs(:,f2); % K x 1
						for aa = 1:NumXi
							XY_XFE{f2}(aa,:) = accumarray(X_FE_j, XY2(aa,:)', [NumX_FE_vals_jf 1])';
						end
					end
					
					% Compute XFE_YX
					for f1 = 1:NumX_FEs_i
						NumX_FE_vals_if = Num_FE_vals_i(f1);
						X_FE_i = Xparts{ii}.X_FEs(:,f1); % T x 1
						for aa = 1:NumXj
							XFE_YX{f1}(:,aa) = accumarray(X_FE_i, XY(:,aa), [NumX_FE_vals_if 1]);
						end
					end
					
					% Compute XFE_Y_XFE
					for f1 = 1:NumX_FEs_i
						NumX_FE_vals_if = Num_FE_vals_i(f1);
						X_FE_i = Xparts{ii}.X_FEs(:,f1); % T x 1
						X_FE_i = reshape(repmat(X_FE_i, [1 K]), [T*K 1]); % (T*K) x 1
						for f2 = 1:NumX_FEs_j
							NumX_FE_vals_jf = Num_FE_vals_j(f2);
							X_FE_j = Xparts{jj}.X_FEs(:,f2); % K x 1
							X_FE_j = reshape(repmat(X_FE_j', [T 1]), [T*K 1]); % (T*K) x 1
							XFE_Y_XFE{f1,f2} = accumarray([X_FE_i X_FE_j], Y, [NumX_FE_vals_if NumX_FE_vals_jf]);
						end
					end
				end
				
				
				%%%% Store everything
				% Store for {jj,ii}
				obj.XY_X      = XY_X;
				obj.XY_XFE    = XY_XFE;
				obj.XFE_YX    = XFE_YX;
				obj.XFE_Y_XFE = XFE_Y_XFE;
				Xt_X_Y_parts{ii,jj} = obj;
				
				% Take transpose and store for {jj,ii}
				obj2.XY_X = obj.XY_X';
				obj2.XY_XFE = obj.XFE_YX;
				for xx = 1:length(obj2.XY_XFE)
					obj2.XY_XFE{xx} = obj2.XY_XFE{xx}';
				end
				obj2.XFE_YX = obj.XY_XFE;
				for xx = 1:length(obj2.XFE_YX)
					obj2.XFE_YX{xx} = obj2.XFE_YX{xx}';
				end
				obj2.XFE_Y_XFE = obj.XFE_Y_XFE';
				for f1 = 1:size(obj2.XFE_Y_XFE,1)
				for f2 = 1:size(obj2.XFE_Y_XFE,2)
					obj2.XFE_Y_XFE{f1,f2} = obj2.XFE_Y_XFE{f1,f2}';
				end
				end
				Xt_X_Y_parts{jj,ii} = obj2;
			end
		end
	end
	
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
