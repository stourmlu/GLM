using SpecialFunctions
include("MNL.jl")
include("../../../Optimization/Julia/Optimization.jl")

using Random, Distributions

L = 15
N = 4000
J = 80
dim1 = N*J

X = rand(Normal(0, 1), dim1, L)
M = 100*ones(Int, N)

# Define choice set id
nn_vec = repeat(collect(1:N), inner=J) # This is full cartesian product

beta_true = rand(Normal(0,1), L, 1)

####################################################################################
println("GENERATING DATA...")
### Generate Y
V = X * beta_true # dim1 x 1

Vmax = zeros(N)
Vmax = accumarray!(Vmax, nn_vec, V, "maximum") # N x 1
V = V - Vmax[nn_vec] # dim1 x 1

tmp = zeros(N)
tmp = log.(accumarray!(tmp, nn_vec, exp.(V))) # N x 1
logp = V - tmp[nn_vec] # dim1 x 1

p = exp.(logp) # dim1 x 1

Y = zeros(Int,dim1)
for nn = 1:N
	idces_nn = findall(nn_vec .== nn) # Jn x 1
	p_nn = p[idces_nn] # Jn x 1
	M_nn = M[nn] # scalar
	Y[idces_nn] = rand(Multinomial(M_nn, p_nn)) # Jn x 1
end

####################################################################################
println("ESTIMATING MODEL...")

### Estimate model
(beta_star, LL_star, LL_grad, FisherInfo, beta_ses) = estimate_MNL(M, X, Y, nn_vec)

[beta_true beta_star beta_ses]
