using Random, Distributions, SpecialFunctions
include("GLM_univariate.jl")
include("../../../Optimization/Julia/Optimization.jl")

# Set dimensions
NumObs = 100000
L = 10

# Draw X and beta
X = [ones(NumObs,1) rand(Normal(0, 0.5), NumObs, L-1)]
beta_true = rand(Normal(0,1), L, 1)

### Compute V
V = compute_Xbeta(X, beta_true)


##### LINEAR MODEL
println("LINEAR MODEL")
Y = V + rand(Normal(0, 0.5), NumObs, 1)
beta_hat, sigmasq_hat, Rsq, std_errors, t_values, p_values, F_stat, residuals = estimate_OLS(X, Y, verbose=true, outputAll=true)
println([beta_hat std_errors t_values p_values])


##### POISSON
println("POISSON MODEL")
logLambdaOffset = -4 .+ 0.5 .* rand(Normal(0, 1), NumObs, 1)
logLambdaOffset = logLambdaOffset[:,1]
lambda = exp.(logLambdaOffset .+ V)
Y = [rand(Poisson(lambdaval)) for lambdaval in lambda]
Y = Y[:,1]
if maximum(Y) > 1e10
	println(maximum(Y))
	throw("Large values")
end
(beta_star, LL_star, LL_grad, FisherInfo, beta_ses) = estimate_GLM(X, Y, "Poisson", extraArg=logLambdaOffset)
println([beta_true beta_star beta_ses])

##### BINOMIAL LOGISTIC REGRESSION
println("BINOMIAL LOGISTIC MODEL")
NumTries = 1000*ones(NumObs)
p = 1 ./ (1 .+ exp.(-V))
Y = [rand(Binomial(NumTries[oo], p[oo])) for oo in 1:NumObs]

(beta_star, LL_star, LL_grad, FisherInfo, beta_ses) = estimate_GLM(X, Y, "BinomialLogistic", extraArg=NumTries)
println([beta_true beta_star beta_ses])
