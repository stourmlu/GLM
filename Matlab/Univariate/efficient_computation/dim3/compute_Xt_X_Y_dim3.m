function [Xt_X_Y, Y] = compute_Xt_X_Y2_dim3(dims, Xparts, Y)
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
% Full data: (T*K1*K2) x 1
% Xpart 1:   (T*K1*K2)  x 1
% Xpart 2:   (T*K1)    x 1
% Xpart 3:   (T*K2)    x 1
% Xpart 4:   (K1*K2)   x 1
% Xpart 5:   T         x 1
% Xpart 6:   K1        x 1
% Xpart 7:   K2        x 1
	
	
	NumParts = dims.NumParts;
	T  = dims.dims1{5};
	K1 = dims.dims1{6};
	K2 = dims.dims1{7};
	
	% Reshape Y and compute aggregate versions
	Y          = reshape(Y,[T K1 K2]);             % T x K1 x K2
	Ysum_tk1   = sum(Y,3);                         % T x K1
	Ysum_tk2   = sum(Y,2);                         % T x  1 x K2
	Ysum_k1_k2 = reshape(sum(Y,1), [K1*K2 1]);     % (K1*K2) x 1
	Ysum_t     = sum(Ysum_tk1,2);                  % T x  1
	Ysum_k1    = sum(Ysum_tk1,1)';                 % K1 x 1
	Ysum_k2    = reshape(sum(Ysum_tk2,1), [K2 1]); % K2 x 1
	Ysum_tk1   = reshape(Ysum_tk1, [T*K1 1]);      % (T*K1) x 1
	Ysum_tk2   = reshape(Ysum_tk2, [T*K2 1]);      % (T*K2) x 1
	
	Xt_X_Y_parts = cell(NumParts, NumParts);
	for ii = 1:NumParts
		Xi            = Xparts{ii}.X;
		NumXi         = dims.Xpart_2_NumX(ii); % TO DO
		if isfield(Xparts{ii}, 'X_FEs')
			X_FEs_i       = Xparts{ii}.X_FEs;
			NumX_FEs_i    = dims.Xpart_2_NumX_FEs(ii); % TO DO
			Num_FE_vals_i = dims.Xpart_2_Num_FE_vals{ii}; % TO DO
		else
			X_FEs_i       = zeros(dims.dims1{ii}, 0);
			NumX_FEs_i    = 0;
			Num_FE_vals_i = {};
		end
		
		if ii == 1
			Y  = reshape(Y,[T*K1*K2 1]);            % (T*K1*K2) x 1
			Ypart_ii = Y;
		end
		if ii == 2
			Ypart_ii = Ysum_tk1;   % (T*K1) x 1
		end
		if ii == 3
			Ypart_ii = Ysum_tk2;   % (T*K2) x 1
		end
		if ii == 4
			Ypart_ii = Ysum_k1_k2; % (K1*K2) x 1
		end
		if ii == 5
			Ypart_ii = Ysum_t;     % T x 1
		end
		if ii == 6
			Ypart_ii = Ysum_k1;    % K1 x 1
		end
		if ii == 7
			Ypart_ii = Ysum_k2;    % K2 x 1
		end
		XY_i = Ypart_ii .* Xi; % dim1_i x NumXi
		
		
		for jj = ii:NumParts
			%%%% Case where ii == jj
			if ii == jj
				obj = compute_Xt_X_Y_samePart(Xparts{ii}, Ypart_ii);
				clear Ypart_ii;
				if ii == 1
					Y  = reshape(Y,[T K1 K2]); % T x K1 x K2
				end
				
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
				
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				if ii == 1
					XY_i = reshape(XY_i, [T K1 K2 NumXi]);
				end
				if ii == 1 && jj == 2 % (T*K1*K2) and (T*K1)
					% Compute XY_X
					XY = sum(XY_i, 3); % T x K1 x 1 x NumXi
					XY = reshape(XY, [T*K1 NumXi]);              % (T*K1) x NumXi
					XY_X = XY' * Xj;                             % NumXi x NumXj
					
					if NumX_FEs_i || NumX_FEs_j  > 0 ; error('Current code does not support FEs at this level'); end; % TO DO
				end
				
				if ii == 1 && jj == 3 % (T*K1*K2) and (T*K2)
					% Compute XY_X
					XY = sum(XY_i, 2); % T x 1 x K2 x NumXi
					XY = reshape(XY, [T*K2 NumXi]);              % (T*K2) x NumXi
					XY_X = XY' * Xj;                             % NumXi x NumXj
					
					if NumX_FEs_i || NumX_FEs_j  > 0 ; error('Current code does not support FEs at this level'); end; % TO DO
				end
				
				if ii == 1 && jj == 4 % (T*K1*K2) and (K1*K2)
					% Compute XY_X
					XY = sum(XY_i, 1); % 1 x K1 x K2 x NumXi
					XY = reshape(XY, [K1*K2 NumXi]);             % (K1*K2) x NumXi
					XY_X = XY' * Xj;                             % NumXi x NumXj
					
					if NumX_FEs_i || NumX_FEs_j  > 0 ; error('Current code does not support FEs at this level'); end; % TO DO
				end
				
				if ii == 1 && jj == 5 % (T*K1*K2) and T
					% Compute XY_X
					XY = sum(XY_i, 2:3); % T x 1 x 1 x NumXi
					XY = reshape(XY, [T NumXi]);                 % T x NumXi
					XY_X = XY' * Xj;                             % NumXi x NumXj
					
					if NumX_FEs_i || NumX_FEs_j  > 0 ; error('Current code does not support FEs at this level'); end; % TO DO
				end
				
				if ii == 1 && jj == 6 % (T*K1*K2) and K1
					% Compute XY_X
					XY = sum(XY_i, [1 3]); % 1 x K1 x 1 x NumXi
					XY = reshape(XY, [K1 NumXi]);                % K1 x NumXi
					XY_X = XY' * Xj;                             % NumXi x NumXj
					
					if NumX_FEs_i || NumX_FEs_j  > 0 ; error('Current code does not support FEs at this level'); end; % TO DO
				end
				
				if ii == 1 && jj == 7 % (T*K1*K2) and K2
					% Compute XY_X
					XY = sum(XY_i, [1 2]); % 1 x 1 x K2 x NumXi
					XY = reshape(XY, [K2 NumXi]);                % K2 x NumXi
					XY_X = XY' * Xj;                             % NumXi x NumXj
					
					if NumX_FEs_i || NumX_FEs_j  > 0 ; error('Current code does not support FEs at this level'); end; % TO DO
				end
				if ii == 1
					XY_i = reshape(XY_i, [T*K1*K2 NumXi]);
				end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
				%%% One has two dims, the other has one dim that is included in the other
				if ii == 2 && jj == 5 % (T*K1) and T
					% Compute XY_X
					XY = sum(reshape(XY_i, [T K1 NumXi]), 2); % T x 1 x NumXi
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
					
					% Compute XFE_YX
					if NumX_FEs_i > 0; error('Current code does not support FEs at tt-k1 level'); end;
					
					% Compute XFE_Y_XFE
					%% Would require NumX_FEs_i > 0
				end
				
				if ii == 2 && jj == 6 % (T*K1) and K1
					% Compute XY_X
					XY = sum(reshape(XY_i, [T K1 NumXi]), 1); % 1 x K1 x NumXi
					XY = reshape(XY, [K1 NumXi]);               % K1 x NumXi
					XY_X = XY' * Xj;                           % NumXi x NumXj
					
					% Compute XY_XFE
					for f2 = 1:NumX_FEs_j
						NumX_FE_vals_jf = Num_FE_vals_j(f2);
						X_FE_j = Xparts{jj}.X_FEs(:,f2); % K1 x 1
						for aa = 1:NumXi
							XY_XFE{f2}(aa,:) = accumarray(X_FE_j, XY(:,aa), [NumX_FE_vals_jf 1])';
						end
					end
					
					% Compute XFE_YX
					if NumX_FEs_i > 0; error('Current code does not support FEs at tt-k1 level'); end;
					
					% Compute XFE_Y_XFE
					%% Would require NumX_FEs_i > 0
				end
				
				if ii == 3 && jj == 5 % (T*K2) and T
					% Compute XY_X
					XY = sum(reshape(XY_i, [T K2 NumXi]), 2); % T x 1 x NumXi
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
					
					% Compute XFE_YX
					if NumX_FEs_i > 0; error('Current code does not support FEs at tt-k2 level'); end;
					
					% Compute XFE_Y_XFE
					%% Would require NumX_FEs_i > 0
					
				end
				
				if ii == 3 && jj == 7 % (T*K2) and K2
					% Compute XY_X
					XY = sum(reshape(XY_i, [T K2 NumXi]), 1); % 1 x K2 x NumXi
					XY = reshape(XY, [K2 NumXi]);               % K2 x NumXi
					XY_X = XY' * Xj;                           % NumXi x NumXj
					
					% Compute XY_XFE
					for f2 = 1:NumX_FEs_j
						NumX_FE_vals_jf = Num_FE_vals_j(f2);
						X_FE_j = Xparts{jj}.X_FEs(:,f2); % K2 x 1
						for aa = 1:NumXi
							XY_XFE{f2}(aa,:) = accumarray(X_FE_j, XY(:,aa), [NumX_FE_vals_jf 1])';
						end
					end
					
					% Compute XFE_YX
					if NumX_FEs_i > 0; error('Current code does not support FEs at tt-k2 level'); end;
					
					% Compute XFE_Y_XFE
					%% Would require NumX_FEs_i > 0
				end

				if ii == 4 && jj == 6 % (K1*K2) and K1
					% Compute XY_X
					XY = sum(reshape(XY_i, [K1 K2 NumXi]), 2); % K1 x 1 x NumXi
					XY = reshape(XY, [K1 NumXi]);              % K1 x NumXi
					XY_X = XY' * Xj;                           % NumXi x NumXj
					
					% Compute XY_XFE
					for f2 = 1:NumX_FEs_j
						NumX_FE_vals_jf = Num_FE_vals_j(f2);
						X_FE_j = Xparts{jj}.X_FEs(:,f2); % K1 x 1
						for aa = 1:NumXi
							XY_XFE{f2}(aa,:) = accumarray(X_FE_j, XY(:,aa), [NumX_FE_vals_jf 1])';
						end
					end
					
					% Compute XFE_YX
					if NumX_FEs_i > 0; error('Current code does not support FEs at k1-k2 level'); end;
					
					% Compute XFE_Y_XFE
					%% Would require NumX_FEs_i > 0
				end
				
				if ii == 4 && jj == 7 % (K1*K2) and K2
					% Compute XY_X
					XY = sum(reshape(XY_i, [K1 K2 NumXi]), 1); % 1 x K2 x NumXi
					XY = reshape(XY, [K2 NumXi]);               % K2 x NumXi
					XY_X = XY' * Xj;                           % NumXi x NumXj
					
					% Compute XY_XFE
					for f2 = 1:NumX_FEs_j
						NumX_FE_vals_jf = Num_FE_vals_j(f2);
						X_FE_j = Xparts{jj}.X_FEs(:,f2); % K2 x 1
						for aa = 1:NumXi
							XY_XFE{f2}(aa,:) = accumarray(X_FE_j, XY(:,aa), [NumX_FE_vals_jf 1])';
						end
					end
					
					% Compute XFE_YX
					if NumX_FEs_i > 0; error('Current code does not support FEs at k1-k2 level'); end;
					
					% Compute XFE_Y_XFE
					%% Would require NumX_FEs_i > 0
				end

				
				
				%%% Both have one dimension, they are different
				if ii == 5 && jj == 6 % T and K1
					% Compute XY_X
					XY = reshape(Ysum_tk1, [T K1]) * Xj; % T x NumXj
					XY_X = Xi' * XY; % NumXi x NumXj
					XY2 = Xi' * reshape(Ysum_tk1, [T K1]); % NumXi x K1
%					XY_X2 = XY2 * Xj; % NumXi x NumXj
					
					% Compute XY_XFE
					for f2 = 1:NumX_FEs_j
						NumX_FE_vals_jf = Num_FE_vals_j(f2);
						X_FE_j = Xparts{jj}.X_FEs(:,f2); % K1 x 1
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
						X_FE_i = reshape(repmat(X_FE_i, [1 K1]), [T*K1 1]); % (T*K1) x 1
						for f2 = 1:NumX_FEs_j
							NumX_FE_vals_jf = Num_FE_vals_j(f2);
							X_FE_j = Xparts{jj}.X_FEs(:,f2); % K1 x 1
							X_FE_j = reshape(repmat(X_FE_j', [T 1]), [T*K1 1]); % (T*K1) x 1
							XFE_Y_XFE{f1,f2} = accumarray([X_FE_i X_FE_j], Ysum_tk1, [NumX_FE_vals_if NumX_FE_vals_jf]);
						end
					end
				end
				
				if ii == 5 && jj == 7 % T and K2
					% Compute XY_X
					XY = reshape(Ysum_tk2, [T K2]) * Xj; % T x NumXj
					XY_X = Xi' * XY; % NumXi x NumXj
					XY2 = Xi' * reshape(Ysum_tk2, [T K2]); % NumXi x K2
%					XY_X2 = XY2 * Xj; % NumXi x NumXj
					
					% Compute XY_XFE
					for f2 = 1:NumX_FEs_j
						NumX_FE_vals_jf = Num_FE_vals_j(f2);
						X_FE_j = Xparts{jj}.X_FEs(:,f2); % K2 x 1
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
						X_FE_i = reshape(repmat(X_FE_i, [1 K2]), [T*K2 1]); % (T*K2) x 1
						for f2 = 1:NumX_FEs_j
							NumX_FE_vals_jf = Num_FE_vals_j(f2);
							X_FE_j = Xparts{jj}.X_FEs(:,f2); % K2 x 1
							X_FE_j = reshape(repmat(X_FE_j', [T 1]), [T*K2 1]); % (T*K2) x 1
							XFE_Y_XFE{f1,f2} = accumarray([X_FE_i X_FE_j], Ysum_tk2, [NumX_FE_vals_if NumX_FE_vals_jf]);
						end
					end
				end
				
				if ii == 6 && jj == 7 % K1 and K2
					% Compute XY_X
					XY = reshape(Ysum_k1_k2, [K1 K2]) * Xj; % K1 x NumXj
					XY_X = Xi' * XY; % NumXi x NumXj
					XY2 = Xi' * reshape(Ysum_k1_k2, [K1 K2]); % NumXi x K2
%					XY_X2 = XY2 * Xj; % NumXi x NumXj
					
					% Compute XY_XFE
					for f2 = 1:NumX_FEs_j
						NumX_FE_vals_jf = Num_FE_vals_j(f2);
						X_FE_j = Xparts{jj}.X_FEs(:,f2); % K2 x 1
						for aa = 1:NumXi
							XY_XFE{f2}(aa,:) = accumarray(X_FE_j, XY2(aa,:)', [NumX_FE_vals_jf 1])';
						end
					end
					
					% Compute XFE_YX
					for f1 = 1:NumX_FEs_i
						NumX_FE_vals_if = Num_FE_vals_i(f1);
						X_FE_i = Xparts{ii}.X_FEs(:,f1); % K1 x 1
						for aa = 1:NumXj
							XFE_YX{f1}(:,aa) = accumarray(X_FE_i, XY(:,aa), [NumX_FE_vals_if 1]);
						end
					end
					
					% Compute XFE_Y_XFE
					for f1 = 1:NumX_FEs_i
						NumX_FE_vals_if = Num_FE_vals_i(f1);
						X_FE_i = Xparts{ii}.X_FEs(:,f1); % K1 x 1
						X_FE_i = reshape(repmat(X_FE_i, [1 K2]), [K1*K2 1]); % (K1*K2) x 1
						for f2 = 1:NumX_FEs_j
							NumX_FE_vals_jf = Num_FE_vals_j(f2);
							X_FE_j = Xparts{jj}.X_FEs(:,f2); % K2 x 1
							X_FE_j = reshape(repmat(X_FE_j', [K1 1]), [K1*K2 1]); % (K1*K2) x 1
							XFE_Y_XFE{f1,f2} = accumarray([X_FE_i X_FE_j], Ysum_k1_k2, [NumX_FE_vals_if NumX_FE_vals_jf]);
						end
					end
					
				end
				
				%%% One has two dims, the other has one dim that is NOT included in the other
				if ii == 4 && jj == 5 % (K1*K2) and T
					% Compute XY_X
					XY = reshape(Y,[T, K1*K2]) * Xi; % T x NumXi
					XY_X = XY' * Xj; % NumXi x NumXj
					
					% Compute XY_XFE
					for f2 = 1:NumX_FEs_j
						NumX_FE_vals_jf = Num_FE_vals_j(f2);
						X_FE_j = Xparts{jj}.X_FEs(:,f2); % T x 1
						for aa = 1:NumXi
							XY_XFE{f2}(aa,:) = accumarray(X_FE_j, XY(:,aa), [NumX_FE_vals_jf 1])';
						end
					end
					
					% Compute XFE_YX
					if NumX_FEs_i > 0; error('Current code does not support FEs at k1-k2 level'); end;
					
					% Compute XFE_Y_XFE
					%% Would require NumX_FEs_i > 0
				end
				
				if ii == 2 && jj == 7 % (T*K1) and K2
					% Compute XY_X
					XY = Xi' * reshape(Y, [T*K1 K2]); % NumXi x K2
					XY_X = XY * Xj; % NumXi x NumXj
					
					% Compute XY_XFE
					for f2 = 1:NumX_FEs_j
						NumX_FE_vals_jf = Num_FE_vals_j(f2);
						X_FE_j = Xparts{jj}.X_FEs(:,f2); % K2 x 1
						for aa = 1:NumXi
							XY_XFE{f2}(aa,:) = accumarray(X_FE_j, XY(aa,:)', [NumX_FE_vals_jf 1])';
						end
					end
					
					% Compute XFE_YX
					if NumX_FEs_i > 0; error('Current code does not support FEs at tt-k1 level'); end;
					
					% Compute XFE_Y_XFE
					%% Would require NumX_FEs_i > 0
				end
				
				if ii == 3 && jj == 6 % (T*K2) and K1
					% Compute XY_X
					XY = zeros(NumXi,K1); % NumXi x K1
					for aa = 1:NumXi
						XY(aa,:) = sum(reshape(Xi(:,aa), [T 1 K2]) .* Y, [1 3]); % 1 x K1
					end
					XY_X = XY * Xj; % NumXi x NumXj
					
					% Compute XY_XFE
					for f2 = 1:NumX_FEs_j
						NumX_FE_vals_jf = Num_FE_vals_j(f2);
						X_FE_j = Xparts{jj}.X_FEs(:,f2); % K1 x 1
						for aa = 1:NumXi
							XY_XFE{f2}(aa,:) = accumarray(X_FE_j, XY(aa,:)', [NumX_FE_vals_jf 1])';
						end
					end
					
					% Compute XFE_YX
					if NumX_FEs_i > 0; error('Current code does not support FEs at tt-k2 level'); end;
					
					% Compute XFE_Y_XFE
					%% Would require NumX_FEs_i > 0
				end
	
				%%% Both have two dimensions, with one overlapping dim
				if ii == 2 && jj == 3 % (T*K1) and (T*K2)
					% Compute XY_X
					XY_X = zeros(NumXi, NumXj);
					for aa = 1:NumXi
						XY_a = sum(reshape(Xi(:,aa), [T K1]) .* Y, 2); % T x 1 x K2
						XY_a = reshape(XY_a, [1, T*K2]); % 1 x (T*K2)
						XY_X(aa,:) = XY_a * Xj; % 1 x NumXj
					end
					
					% Compute XY_XFE
					if NumX_FEs_j > 0; error('Current code does not support FEs at tt-k2 level'); end;
					
					% Compute XFE_YX
					if NumX_FEs_i > 0; error('Current code does not support FEs at tt-k1 level'); end;
					
					% Compute XFE_Y_XFE
					%% Would require NumX_FEs_i > 0 and NumX_FEs_j > 0
				end
				
				if ii == 2 && jj == 4 % (T*K1) and (K1*K2)
					% Compute XY_X
					XY_X = zeros(NumXi, NumXj);
					for bb = 1:NumXj
						XY_b = sum(reshape(Xj(:,bb), [1 K1 K2]) .* Y, 3); % T x K1
						XY_b = reshape(XY_b, [1 T*K1]); % 1 x (T*K1)
						XY_X(:,bb) = (XY_b * Xi)'; % NumXi x 1
					end
					
					% Compute XY_XFE
					if NumX_FEs_j > 0; error('Current code does not support FEs at tt-k1 level'); end;
					
					% Compute XFE_YX
					if NumX_FEs_i > 0; error('Current code does not support FEs at k1-k2 level'); end;
					
					% Compute XFE_Y_XFE
					%% Would require NumX_FEs_i > 0 and NumX_FEs_j > 0
				end
				
				if ii == 3 && jj == 4 % (T*K2) and (K1*K2)
					% Compute XY_X
					XY_X = zeros(NumXi, NumXj);
					for bb = 1:NumXj
						XY_b = sum(reshape(Xj(:,bb), [1 K1 K2]) .* Y, 2); % T x 1 x K2
						XY_b = reshape(XY_b, [1 T*K2]); % 1 x (T*K2)
						XY_X(:,bb) = (XY_b * Xi)'; % NumXi x 1
					end
					
					% Compute XY_XFE
					if NumX_FEs_j > 0; error('Current code does not support FEs at tt-k2 level'); end;
					
					% Compute XFE_YX
					if NumX_FEs_i > 0; error('Current code does not support FEs at k1-k2 level'); end;
					
					% Compute XFE_Y_XFE
					%% Would require NumX_FEs_i > 0 and NumX_FEs_j > 0
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
	
	% Reshape Y as it was initially
	Y = reshape(Y,[T*K1*K2 1]); % (T*K1*K2) x 1
end
