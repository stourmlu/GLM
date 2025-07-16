using SpecialFunctions

function compute_Xbeta(X, beta_vec)
	###########################################################################
	# 
	###########################################################################
	##### Inputs:
	# X
	###########################################################################
	##### Outputs:
	# Xbeta:				NumObs x 1
	###########################################################################
	
	return X * beta_vec
	# TO DO: handle case when X is split into different parts, by dimensions of variation
end


function compute_Xt_Y(X, Y)
	###########################################################################
	# This function computes sum_{n} Y_{n} * X_{na} for all a.
	# This function returns the result, combining everything together.
	#
	# I pass Y as an output argument for memory optimization reason:
	# - I need to "update" Y by reshaping it
	# - by doing so, Matlab would normally make a local copy of Y (which is memory heavy + takes time)
	# - by letting Y be an output variable (and in the calling function, "update" Y accordingly),
	#		I tell Matlab that it's okay to directly edit the same copy of Y (no need to make a local copy)
	# - for this to work, I also need to "update" Y when calling this function
	# - I only reshape Y back and forth here, so Y is identical at the start and at the end, so it is safe.
	###########################################################################
	##### Inputs:
	# data
	# Y:					NumObs x 1
	###########################################################################
	##### Outputs:
	# Xt_Y:				1 x NumParams
	###########################################################################
	
	return Y' * X
	# TO DO: handle case when X is split into different parts, by dimensions of variation
end



function compute_Xt_X_Y(X, Y)
	###########################################################################
	# This function computes sum_{n} Y_{n} * X_{na} * X_{nb} for all (a,b).
	# This function returns the result, combining everything together.
	#
	# I pass Y as an output argument for memory optimization reason:
	# - I need to "update" Y by reshaping it
	# - by doing so, Matlab would normally make a local copy of Y (which is memory heavy + takes time)
	# - by letting Y be an output variable (and in the calling function, "update" Y accordingly),
	#		I tell Matlab that it's okay to directly edit the same copy of Y (no need to make a local copy)
	# - for this to work, I also need to "update" Y when calling this function
	# - I only reshape Y back and forth here, so Y is identical at the start and at the end, thus it is safe.
	###########################################################################
	##### Inputs:
	# X
	# Y:					NumObs x 1
	###########################################################################
	##### Outputs:
	# Xt_X_Y:			NumParams x NumParams
	###########################################################################
	
	return X' * (X .* Y)
	# TO DO: handle case when X is split into different parts, by dimensions of variation
end



function estimate_OLS(X, Y; verbose=true, outputAll=true)
	##########################################
	##### Inputs:
	#	data
	##########################################
	##### Outputs:
	#	beta_hat:			NumX x 1
	# varargout:
	#   sigmasq_hat:	scalar
	#   Rsq:			scalar
	#   std_errors:		NumX x 1
	#   t_values:		NumX x 1
	#   p_values:		NumX x 1
	#   F_stat:			scalar
	#   residuals:		NumObs x 1
	##########################################

	g = compute_Xt_Y(X, Y)'  # NumX x 1
	
	NumObs = size(X,1)
	H = compute_Xt_X_Y(X, ones(NumObs,1)) # NumX x NumX
	
	beta_hat = H\g
	
	
	if verbose || outputAll # Calculate sigmasq_hat: estimate of sigma^2
		NumX = size(beta_hat,1)
		residuals = Y - compute_Xbeta(X, beta_hat)
		sigmasq_hat = sum(residuals.^2)/(NumObs - NumX)
	end
	
	if verbose || outputAll # Compute R^2
		SSE = sum(residuals.^2)
		SST = sum((Y .- mean(Y)).^2) # TO DO
		SSR = SST - SSE
		Rsq = 1 .- SSE/SST
	end

	if verbose || outputAll # Obtain standard errors
		std_errors = sqrt.(diag(sigmasq_hat * inv(H)))
	end

	if verbose || outputAll # Obtain t-values
		t_values = beta_hat./std_errors
	end

	if verbose || outputAll # Obtain p-values
		p_values = 2.0 .*(1 .- cdf.(TDist(NumObs-NumX), abs.(t_values)))
	end
	
	if verbose || outputAll # Obtain F statistic
		F_stat = (SSR/(NumX-1))/(SSE/(NumObs - NumX))
	end
	
	if verbose || outputAll # Output residuals
		residuals = residuals
	end
	
	if verbose
		mytable = [beta_hat,std_errors, t_values, p_values]
	end
	
	if outputAll
		return beta_hat, sigmasq_hat, Rsq, std_errors, t_values, p_values, F_stat, residuals
	else
		return beta_hat
	end
end



function GLM_univ_loglikelihood(modelId, X, Y, beta_val; LL_constant=0, takeNegative=false, extraArg=nothing, 
	returnGrad::Bool=false, returnHessian::Bool=false, returnNegative::Bool=false)
	###########################################################################
	# No description available for function GLM_univ_loglikelihood.
	###########################################################################
	##### Inputs:
	# modelId:			       integer (1 => BinomialLogit, 2 => Poisson)
	# X:					   NumObs x L
	# Y:                       NumObs x 1
	# beta_val:                NumBetas x 1
	# takeNegative:            boolean
	###########################################################################
	##### Outputs:
	# LL:            		scalar
	# grad:					NumParams x 1
	# hessian:				NumParams x NumParams
	###########################################################################
	
	### Compute V
	V = compute_Xbeta(X, beta_val)	# NumObs x 1

	
	### Compute log-likelihood
	if modelId == 1  # Logistic-Binomial
		NumTries = extraArg			# NumObs x 1
		p  = 1 ./(1 .+ exp.(-V))			# NumObs x 1
		LL = LL_constant + Y'*V + NumTries'*log.(1 .- p)	# scalar
	end
	if modelId == 2  # Poisson
		logLambda  = extraArg + V	# NumObs x 1
		lambda     = exp.(logLambda);	# NumObs x 1
		LL         = LL_constant + Y'*logLambda .- sum(lambda) # scalar
	end
	LL = LL[1]
	if takeNegative; LL = -LL; end;
	
	if !returnGrad
		return LL
	end
	
	### Compute gradient of log-likelihood
	# Compute dLL_dV [NumObs x 1]
	if modelId == 1; dLL_dV = Y .- NumTries .* p; end;		# Logistic-Binomial
	if modelId == 2; dLL_dV = Y .- lambda; end; 	# Poisson
	grad = compute_Xt_Y(X, dLL_dV)' # NumParams x 1))
	if takeNegative; grad = -grad; end;
	
	if !returnHessian
		return LL, grad
	end
	
	### Compute hessian of log-likelihood
	# Compute d2LL_dV2 [NumObs x 1]
	if modelId == 1; d2LL_dV2 = -NumTries .* p .* (1 .- p); end	# Logistic-Binomial
	if modelId == 2; d2LL_dV2 = -lambda; end					# Poisson
	hessian = compute_Xt_X_Y(X, d2LL_dV2)	# NumParams x NumParams
	if takeNegative; hessian = -hessian; end;
	
	return LL, grad, hessian
end



function estimate_GLM(X, Y, modelName; extraArg=nothing, beta0=nothing, verbose=true)
	NumParams = size(X,2)

	if beta0 == nothing
		beta0 = zeros(NumParams)	
	end

	### Case of Binomial Logistic regression
	if modelName == "BinomialLogistic"
		NumTries = extraArg
		LL_constant = sum(loggamma.(NumTries .+ 1)) - sum(loggamma.(Y .+ 1)) - sum(loggamma.(NumTries .- Y .+ 1))
		modelId = 1
	end	
	
	### Case of Poisson regression
	if modelName == "Poisson"
		logLambdaOffset = extraArg
		LL_constant = Y'*logLambdaOffset .- sum(loggamma.(Y .+ 1))
		modelId = 2
	end

	obj = (betaval; grad::Bool=false, hessian::Bool=false) -> GLM_univ_loglikelihood(modelId, X, Y, betaval, LL_constant=LL_constant, extraArg=extraArg, takeNegative=true, returnGrad=grad, returnHessian=hessian)
	
	obj(beta0)
	
	
	
	### Apply Newton-Raphson algorithm
	###############
	(beta_hat, LL_star, LL_grad, FisherInfo) = Newton_Raphson(beta0, obj, verbose=true)
	###############
	### Useful when Newton-Raphson gets stuck sometimes
#	(beta_hat, LL_star, LL_grad) = BFGS(beta0, obj, verbose=true)
#	(beta_hat, LL_star, LL_grad, FisherInfo) = Newton_Raphson(beta_hat, obj, verbose=true)
	###############
	beta_ses = sqrt.(diag(inv(FisherInfo)))
	
	if verbose
		println([beta_hat,beta_ses])
	end
	
	return beta_hat, LL_star, LL_grad, FisherInfo, beta_ses
end



