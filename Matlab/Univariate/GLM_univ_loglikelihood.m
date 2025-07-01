function [LL, varargout] = GLM_univ_loglikelihood(modelId, data, beta_val, LL_constant, takeNegative, varargin)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% No description available for function GLM_univ_loglikelihood.
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Inputs:
	% modelId:			       integer (1 => BinomialLogit, 2 => Poisson)
	% data:			           object:
	% 	.dims:				     	object:
	%		.mappings:					cell(NumParts,1)
	%			{ii}:						NumObs x 1 (gives integer values between 1 and dim1_ii)
	%		.dims1:						cell(NumParts,1)
	%			{ii}:						integer: gives dim1_ii that corresponds to Xparts{ii}
	%		.Xpart_2_NumX:				1 x NumParts: gives integers (dim2 of Xparts{ii}.X)
	%		.Xpart_2_NumX_FEs:			1 x NumParts: gives integers (dim2 of Xparts{ii}.X_FEs)
	%		.Xpart_2_Num_FE_vals:		cell(1,NumParts)
	%			{ii}:						1 x NumX_FEs_i  --> gives integer: number of possible values for Xparts{ii}.X_FEs(:,ff)
	%		.NumFEvals2Keep:			cell(1,NumParts) --> gives integer
	%		.NumParts
	%		.NumObs
	%		.NumParams
	%	 .Xparts:				   cell(NumParts,1)
	%	 	{ii}:				       object
	%			.X:					       dim1_ii x NumX
	%			.X_FEs:					   dim1_ii x NumX_FEs
	%			.NumX_FE_vals:			   1 x NumX_FEs: gives integer
	%    .Y:                       NumObs x 1
	% beta_val:                NumBetas x 1
	% takeNegative:            boolean
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Outputs:
	% LL:            		scalar
	% grad:					NumParams x 1
	% hessian:				NumParams x NumParams
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%%% Compute V
	V = compute_Xbeta(data, beta_val);	% NumObs x 1
	
	%%% Compute log-likelihood
	if modelId == 1  % Logistic-Binomial
		NumTries = varargin{1};			% NumObs x 1
		p  = 1./(1+exp(-V));			% NumObs x 1
		LL = LL_constant + data.Y'*V + NumTries'*log(1-p);	% scalar
		clear V;
	end
	if modelId == 2  % Poisson
		logLambda  = varargin{1} + V;	% NumObs x 1
		clear V;
		lambda     = exp(logLambda);	% NumObs x 1
		LL         = LL_constant  + data.Y'*logLambda -sum(lambda,1); % scalar
		clear logLambda;
	end
	if takeNegative; LL = -LL; end;
	
	%%% Compute gradient of log-likelihood
	if nargout >= 2
		% Compute dLL_dV [NumObs x 1]
		if modelId == 1; dLL_dV = data.Y - NumTries .* p; end;		% Logistic-Binomial
		if modelId == 2; dLL_dV = data.Y - lambda; end; 	% Poisson
		grad = compute_Xt_Y(data, dLL_dV)'; % NumParams x 1
		clear dLL_dV;
		if takeNegative; grad = -grad; end;
		varargout{1} = grad;
	end
	
	%%% Compute hessian of log-likelihood
	if nargout >= 3
		% Compute d2LL_dV2 [NumObs x 1]
		if modelId == 1; d2LL_dV2 = -NumTries .* p .* (1-p); end	% Logistic-Binomial
		if modelId == 2; d2LL_dV2 = -lambda; end					% Poisson
		hessian = compute_Xt_X_Y(data, d2LL_dV2);	% NumParams x NumParams
		clear d2LL_dV2;
		if takeNegative; hessian = -hessian; end;
		varargout{2} = hessian;
	end
end
