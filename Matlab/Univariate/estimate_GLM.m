function [beta_hat, LL_star, LL_grad, FisherInfo, beta_ses] = estimate_GLM(data, modelName, varargin)
	if length(varargin) >= 1
		extraArg = varargin{1};
	else
		extraArg = {};
	end
	if length(varargin) >= 2
		beta0 = varargin{2};
	else
		beta0 = {};
	end
	if length(varargin) >= 3
		verbose = varargin{3};
	else
		verbose = true;
	end
	
	
	% Read dimensions
	if isfield(data, 'Xparts')
		NumParams = data.dims.NumParams;
	else
		NumParams = size(data.X,2);
	end

	%%% Case of Binomial Logistic regression
	if strcmp(modelName, 'BinomialLogistic')
		NumTries = extraArg;
		LL_constant = sum(gammaln(NumTries+1)) - sum(gammaln(data.Y+1)) - sum(gammaln(NumTries-data.Y+1));	
		obj = @(betaval) GLM_univ_loglikelihood(1, data, betaval, LL_constant, true, NumTries);
	end	
	
	%%% Case of Poisson regression
	if strcmp(modelName, 'Poisson')
		logLambdaOffset = extraArg;
		LL_constant = data.Y'*logLambdaOffset - sum(gammaln(data.Y+1));
		obj = @(betaval) GLM_univ_loglikelihood(2, data, betaval, LL_constant, true, logLambdaOffset);
	end
	
	%%% Apply Newton-Raphson algorithm
	if isempty(beta0)
		beta0 = zeros(NumParams,1);	
	end
	%%%%%%%%%%%%%%%
	[beta_hat, LL_star, LL_grad, FisherInfo] = Newton_Raphson(obj, beta0, 1e-8, 2000, 1, true);
	%%%%%%%%%%%%%%%
	%%% Useful when Newton-Raphson gets stuck sometimes
%	[beta_hat, LL_star, LL_grad] = BFGS(obj, beta0, 1e-8, 10, true);
%	[beta_hat, LL_star, LL_grad, FisherInfo] = Newton_Raphson(obj, beta_hat, 1e-8, 2000, 1, true);
	%%%%%%%%%%%%%%%
	beta_ses = sqrt(diag(inv(FisherInfo)));
	
	if verbose
		mytable = table(beta_hat,beta_ses);
		if isfield(data, 'Xnames')
			mytable.Properties.RowNames = data.Xnames;
		end
		disp(mytable);
	end
end
