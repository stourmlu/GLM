function [LL, varargout] = loglikelihood_MNL(Y, X, M, M_rep, nn_vec, beta_val, LL_constant, returnNegative)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% This function computes the loglikelihood, its gradient and its Hessian for a Multinomial logit model.
	%
	% M_rep should be equal to M(nn_vec). No check is done.
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Inputs:
	% Y:					dim1 x 1
	% X:					dim1 x L
	% M:					N x 1
	% M_rep:				dim1 x 1
	% nn_vec:				dim1 x 1 (integers between 1 and N)
	% beta_val:				L x 1
	% LL_constant:			scalar
	% returnNegative:		boolean
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Outputs:
	% LL:            		scalar
	% grad:					L x 1
	% hessian:				L x L
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%%% Read dimensions
	[dim1, L] = size(X);
	N = size(M,1);
	
	%%% Compute V
	V = X * beta_val; % dim1 x 1
	
	%%% Trick to avoid taking the exponential of large values
	Vmax = accumarray(nn_vec, V, [N 1], @max); % N x 1
	V = V - Vmax(nn_vec); % dim1 x 1
	
	%%% Compute logp
	tmp = log(accumarray(nn_vec, exp(V), [N 1])); % N x 1
	logp = V - tmp(nn_vec); % dim1 x 1
	
	%%% Compute LL
	LL = Y'*logp + LL_constant; % scalar
	if returnNegative
		LL = - LL;
	end
	
	if nargout <= 1; return; end;
	
	%%% Compute p
	p = exp(logp);
	
	%%% Compute grad
	grad = ((Y - M_rep .* p)'*X)'; % L x 1
	if returnNegative
		grad = - grad;
	end
	varargout{1} = grad;
	if nargout <= 2; return; end;
	
	%%% Compute A
	A = zeros(N, L); % N x L
	for ll = 1:L
		A(:,ll) = accumarray(nn_vec, p.*X(:,ll), [N 1]);
	end
	
	%%% Compute Hessian
	hess = -X'*(X.*(M_rep .* p)) + A'*(A.*M); % L x L
	if returnNegative
		hess = - hess;
	end
	varargout{2} = hess;
end
