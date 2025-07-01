function [Xt_X_Y, Y] = compute_Xt_X_Y2_general(dims, Xparts, Y)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% This function computes sum_{n} Y_{n} * X_{na} * X_{nb} for all (a,b).
	% This function returns the result, combining everything together.
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
	%	{ii}:					object
	%		.X:						dim1_ii x NumX
	%		.X_FEs:					dim1_ii x NumX_FEs
	%		.NumX_FE_vals:			1 x NumX_FEs: gives integer
	% Y:					NumObs x 1
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Outputs:
	% Xt_X_Y:			NumParams x NumParams
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This part does not assume that NumObs contains the full cartesian product of different dimensions in any kind of way.
% However, it can be extremely slow!
	
	NumParts = dims.NumParts;

	Xt_X_Y_parts = cell(NumParts, NumParts);
	for ii = 1:NumParts
	for jj = ii:NumParts
		if ii == jj
			myvec = dims.mappings{ii}; % NumObs x 1 (gives integer values between 1 and dim1_ii)
			dim1_ii = dims.dims1{ii};
			Ypart_ii = accumarray(myvec, Y, [dim1_ii 1]); % dim1_ii x 1
			obj = compute_Xt_X_Y_samePart(Xparts{ii}, Ypart_ii);
			Xt_X_Y_parts{ii,ii} = obj;
		else
			NumXi            = dims.Xpart_2_NumX(ii); % TO DO
			NumXj            = dims.Xpart_2_NumX(jj); % TO DO
			NumX_FEs_i      = dims.Xpart_2_NumX_FEs(ii); % TO DO
			NumX_FEs_j      = dims.Xpart_2_NumX_FEs(jj); % TO DO
			
			myvec_ii = dims.mappings{ii}; % NumObs x 1 (gives integer values between 1 and dim1_ii)
			myvec_jj = dims.mappings{jj}; % NumObs x 1 (gives integer values between 1 and dim1_jj)
			
			% Make XY_X
			XY_X = zeros(NumXi, NumXj); % NumXi x NumXj
			for aa = 1:NumXi
				Xia = Xparts{ii}.X(myvec_ii,aa);  % NumObs x 1
				Xia_Y = Xia.*Y; % NumObs x 1
				for bb = 1:NumXj
					Xjb = Xparts{jj}.X(myvec_jj,bb); % NumObs x 1
					XY_X(aa,bb) = Xia_Y'*Xjb;
				end
			end
			
			% Make XY_XFE
			XY_XFE = cell(NumX_FEs_j,1);
			for f2 = 1:NumX_FEs_j
				XY_XFE{f2} = zeros(NumXi, Xparts{jj}.NumX_FE_vals(f2));
			end
			for aa = 1:NumXi
				Xia = Xparts{ii}.X(myvec_ii,aa);  % NumObs x 1
				Xia_Y = Xia.*Y; % NumObs x 1
				for f2 = 1:NumX_FEs_j
					NumX_FE_vals_j = Xparts{jj}.NumX_FE_vals(f2);
					myFE_val = Xparts{jj}.X_FEs(myvec_jj,f2); % NumObs x 1
					XY_XFE{f2}(aa,:) = accumarray(myFE_val, Xia_Y, [NumX_FE_vals_j 1])'; % 1 x NumX_FE_vals_j
				end
			end
			
			% Make XFE_YX
			XFE_YX = cell(NumX_FEs_i,1);
			for f1 = 1:NumX_FEs_i
				XFE_YX{f1} = zeros(Xparts{ii}.NumX_FE_vals(f1), NumXj);
			end
			for bb = 1:NumXj
				Xjb = Xparts{jj}.X(myvec_jj,bb);  % NumObs x 1
				Xjb_Y = Xjb.*Y; % NumObs x 1
				for f1 = 1:NumX_FEs_i
					NumX_FE_vals_i = Xparts{ii}.NumX_FE_vals(f1);
					myFE_val = Xparts{ii}.X_FEs(myvec_ii,f1); % NumObs x 1
					XFE_YX{f1}(:,bb) = accumarray(myFE_val, Xjb_Y, [NumX_FE_vals_i 1]); %  NumX_FE_vals_i x 1
				end
			end
			
			% Make XFE_Y_XFE
			XFE_Y_XFE = cell(NumX_FEs_i, NumX_FEs_j);
			for f1 = 1:NumX_FEs_i
			for f2 = 1:NumX_FEs_j
				NumX_FE_vals_if = dims.Xpart_2_Num_FE_vals{ii}(f1);
				NumX_FE_vals_jf = dims.Xpart_2_Num_FE_vals{jj}(f2);
				myFE_val_ii = Xparts{ii}.X_FEs(myvec_ii,f1); % NumObs x 1
				myFE_val_jj = Xparts{jj}.X_FEs(myvec_jj,f2); % NumObs x 1
				XFE_Y_XFE{f1, f2} = accumarray([myFE_val_ii myFE_val_jj], Y, [NumX_FE_vals_if NumX_FE_vals_jf]); % NumX_FE_vals_if x NumX_FE_vals_jf
			end
			end
			
			% Store for {jj,ii}
			obj.XY_X      = XY_X;
			obj.XY_XFE    = XY_XFE;
			obj.XFE_YX    = XFE_YX;
			obj.XFE_Y_XFE = XFE_Y_XFE;
			Xt_X_Y_parts{ii,jj} = obj;
			
			% Take transpose for {jj,ii}
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
