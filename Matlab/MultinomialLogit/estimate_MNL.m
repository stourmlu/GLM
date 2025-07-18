function [beta_star, LL_star, LL_grad, FisherInfo, beta_ses] = estimate_MNL(M, X, Y, nn_vec)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% This function estimates a multinomial logit model by maximum likelihood, using the Newton-Raphson algorithm.
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Inputs:
	% M:					N x 1        (number of tries for choice set with index nn)
	% X:					dim1 x NumX  (covariates corresponding to each option, across all choice sets)
	% Y:					dim1 x 1     (number of "successful tries" for that option in that choice set, across the M_n tries)
	% nn_vec:				dim1 x 1     (maps each option to the choice set it belongs to, values between 1 and N)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Outputs:
	% beta_star:           	NumX x 1
	% LL_star:				scalar
	% grad:					NumX x 1
	% FisherInfo:			NumX x NumX
	% beta_ses:				NumX x 1
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	M_rep = M(nn_vec);
	LL_constant = sum(gammaln(M+1)) - sum(gammaln(Y+1));

	obj = @(beta_val) loglikelihood_MNL(Y, X, M, M_rep, nn_vec, beta_val, LL_constant, true);
	beta0 = zeros(size(X,2), 1);
	
	[beta_star, LL_star, LL_grad, FisherInfo] = Newton_Raphson(obj, beta0, 1e-6, 2000);	
	beta_ses = sqrt(diag(inv(FisherInfo)));
end
