function [res] = compute_Xt_X_Y_samePart(Xpart, Yval)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% This is only within the same Xpart.
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Inputs:
	% Xpart:				object:
	%	.X:						dim1 x NumX
	%	.X_FEs:					dim1 x NumX_FEs
	%	.NumX_FE_vals:			1 x NumX_FEs: gives integer
	% Yval:					dim1 x 1
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Outputs:
	% res:					object:
	%	.XY_X:					NumX x NumX
	%	.XY_XFE:				cell(NumX_FEs,1)
	%		{ii}:					NumX x NumFE_vals_ii
	%	.XFE_YX:				cell(NumX_FEs,1)
	%		{ii}:					NumFE_vals_ii x NumX
	%	.XFE_Y_XFE:				cell(NumX_FEs,NumX_FEs)
	%		{ii,jj}:				NumFE_vals_ii x NumFE_vals_jj
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	NumX     = size(Xpart.X, 2);

	if isfield(Xpart, 'X_FEs')
		NumX_FEs = size(Xpart.X_FEs, 2);
	else
		NumX_FEs = 0;
	end
	
	% Compute YX
	YX = Yval.*Xpart.X; % dim1 x NumX;
	
	% Compute part with respect to X and X
	XY_X = YX' * Xpart.X; % NumX x NumX
	
	% Compute part with respect to X and X_FEs
	XY_XFE = cell(NumX_FEs, 1);
	XFE_YX = cell(NumX_FEs, 1);
	for ii = 1:NumX_FEs
		NumX_FE_vals_ii = Xpart.NumX_FE_vals(ii);
		XY_XFE{ii} = zeros(NumX, NumX_FE_vals_ii);
		for xx = 1:NumX
			XY_XFE{ii}(xx,:) = accumarray(Xpart.X_FEs(:,ii), YX(:,xx), [NumX_FE_vals_ii 1])';
		end
		XFE_YX{ii} = XY_XFE{ii}';
	end
	
	% Compute part with respect to X_FEs and X_FEs
	XFE_Y_XFE = cell(NumX_FEs, NumX_FEs);
	for ii = 1:NumX_FEs
		% Case 1: same set of FEs (ii = jj)
		NumX_FE_vals_ii = Xpart.NumX_FE_vals(ii);
		tmp = accumarray(Xpart.X_FEs(:,ii), Yval, [NumX_FE_vals_ii 1]); % NumX_FE_vals_ii x 1
		XFE_Y_XFE{ii,ii} = diag(tmp); % NumX_FE_vals_ii x NumX_FE_vals_ii
	
		% Case 2: different sets of FEs (ii != jj => exploit symmetrys)
		for jj = ii+1:NumX_FEs
			NumX_FE_vals_jj = Xpart.NumX_FE_vals(jj);
			tmp2 = accumarray([Xpart.X_FEs(:,ii) Xpart.X_FEs(:,jj)], Yval, [NumX_FE_vals_ii NumX_FE_vals_jj]); % NumFE_vals_ii x NumFE_vals_jj
			XFE_Y_XFE{ii,jj} = tmp2;  % NumFE_vals_ii x NumFE_vals_jj
			XFE_Y_XFE{jj,ii} = tmp2'; % NumFE_vals_jj x NumFE_vals_ii
		end
	end
	
	% Output everything
	res.XY_X      = XY_X;
	res.XY_XFE    = XY_XFE;
	res.XFE_YX    = XFE_YX;
	res.XFE_Y_XFE = XFE_Y_XFE;
end
