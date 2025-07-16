function accumarray!(out, ids, val, func=nothing)
	out[:] .= 0
	if func == nothing # Take the sum
    	for i in eachindex(ids, val)
			out[ids[i]] += val[i]		
		end
	elseif func == "maximum"
    	for i in eachindex(ids, val)
			id_ii = ids[i]
			val_ii = val[i]
			if out[id_ii] < val_ii
				out[id_ii] = val_ii
			end
		end
	else
		throw("Function not implemented")
    end
    return out
end

function loglikelihood_MNL(Y, X, M, nn_vec, beta_val; LL_constant=nothing, M_rep=nothing, vec_N=nothing, 
	returnGrad::Bool=false, returnHessian::Bool=false, returnNegative::Bool=false)
	###########################################################################
	# This function computes the loglikelihood, its gradient and its Hessian for a Multinomial logit model.
	#
	# M_rep should be equal to M(nn_vec). No check is done.
	###########################################################################
	##### Inputs:
	# Y:					dim1 x 1
	# X:					dim1 x L
	# M:					N x 1
	# M_rep:				dim1 x 1
	# nn_vec:				dim1 x 1 (integers between 1 and N)
	# beta_val:				L x 1
	# LL_constant:			scalar
	# Vmax:					N x 1
	# returnNegative:		boolean
	###########################################################################
	##### Outputs:
	# LL:            		scalar
	# grad:					L x 1
	# hessian:				L x L
	###########################################################################
	
	### Read dimensions
	(dim1, L) = size(X)
	N = size(M,1)
	
	### Compute V
	V = X * beta_val # dim1 x 1
	
	if M_rep == nothing
		M_rep = M[nn_vec]
	end
	if vec_N == nothing
		vec_N = zeros(N)
	end
	if LL_constant == nothing
		LL_constant = sum(loggamma.(M.+1)) - sum(loggamma.(Y.+1))
	end
	
	### Trick to avoid taking the exponential of large values
	Vmax = accumarray!(vec_N, nn_vec, V, "maximum") # N x 1
	V = V - Vmax[nn_vec] # dim1 x 1
	
	### Compute logp
	tmp = log.(accumarray!(vec_N, nn_vec, exp.(V))) # N x 1
	logp = V - tmp[nn_vec] # dim1 x 1
	
	### Compute LL
	LL = Y'*logp + LL_constant # scalar
	if returnNegative
		LL = -LL
	end
	if !returnGrad
		return LL
	end
	
	### Compute p
	p = exp.(logp)
	
	### Compute grad
	grad = ((Y - M_rep .* p)'*X)' # L x 1
	if returnNegative
		grad = - grad
	end
	if !returnHessian
		return (LL, grad)
	end
	
	### Compute A
	A = zeros(N, L) # N x L
	for ll = 1:L
		A[:,ll] = accumarray!(vec_N, nn_vec, p.*X[:,ll])
	end
	
	### Compute Hessian
	hess = -X'*(X.*(M_rep .* p)) + A'*(A.*M) # L x L
	if returnNegative
		hess = - hess
	end
	return (LL, grad, hess)
end




function estimate_MNL(M, X, Y, nn_vec)
	###########################################################################
	# This function estimates a multinomial logit model by maximum likelihood, using the Newton-Raphson algorithm.
	###########################################################################
	##### Inputs:
	# M:					N x 1        (number of tries for choice set with index nn)
	# X:					dim1 x NumX  (covariates corresponding to each option, across all choice sets)
	# Y:					dim1 x 1     (number of "successful tries" for that option in that choice set, across the M_n tries)
	# nn_vec:				dim1 x 1     (maps each option to the choice set it belongs to, values between 1 and N)
	###########################################################################
	##### Outputs:
	# beta_star:           	NumX x 1
	# LL_star:				scalar
	# grad:					NumX x 1
	# FisherInfo:			NumX x NumX
	# beta_ses:				NumX x 1
	###########################################################################
	
	M_rep = M[nn_vec]
	vec_N = zeros(N)
	LL_constant = sum(loggamma.(M.+1)) - sum(loggamma.(Y.+1))

	beta0 = zeros(size(X,2))
	obj = (beta_val; grad::Bool=false, hessian::Bool=false) -> loglikelihood_MNL(Y, X, M, nn_vec, beta_val, LL_constant=LL_constant, M_rep=M_rep, vec_N=vec_N, returnNegative=true, returnGrad=grad, returnHessian=hessian)

	(beta_star, LL_star, LL_grad, FisherInfo) = Newton_Raphson(beta0, obj, verbose=true)
	
	beta_ses = sqrt.(diag(inv(FisherInfo)))
	
	return beta_star, LL_star, LL_grad, FisherInfo, beta_ses
end
