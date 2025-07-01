function [beta_star, LL_star, LL_grad, FisherInfo, beta_ses] = estimate_GLM(data, modelName, varargin)
	% Read dimensions
	if isfield(data, 'Xparts')
		NumParams = data.dims.NumParams;
	else
		NumParams = size(data.X,2);
	end

	%%% Case of Binomial Logistic regression
	if strcmp(modelName, 'BinomialLogistic')
		NumTries = varargin{1};
		LL_constant = sum(gammaln(NumTries+1)) - sum(gammaln(data.Y+1)) - sum(gammaln(NumTries-data.Y+1));	
		obj = @(betaval) GLM_univ_loglikelihood(1, data, betaval, LL_constant, true, NumTries);
	end	
	
	%%% Case of Poisson regression
	if strcmp(modelName, 'Poisson')
		logLambdaOffset = varargin{1};
		LL_constant = data.Y'*logLambdaOffset - sum(gammaln(data.Y+1));
		obj = @(betaval) GLM_univ_loglikelihood(2, data, betaval, LL_constant, true, logLambdaOffset);
	end
	
	%%% Apply Newton-Raphson algorithm
	beta0 = zeros(NumParams,1);
	%%%%%%%%%%%%%%%
	[beta_star, LL_star, LL_grad, FisherInfo] = Newton_Raphson(obj, beta0, 1e-8, 2000, 1, true);
	%%%%%%%%%%%%%%%
	%%% Useful when Newton-Raphson gets stuck sometimes
%	[beta_star, LL_star, LL_grad] = BFGS(obj, beta0, 1e-8, 10, true);
%	[beta_star, LL_star, LL_grad, FisherInfo] = Newton_Raphson(obj, beta_star, 1e-8, 2000, 1, true);
	%%%%%%%%%%%%%%%
	beta_ses = sqrt(diag(inv(FisherInfo)));
end
