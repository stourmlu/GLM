L = 15;
N = 4000;
J = 80;
dim1 = N*J;

X = normrnd(0,1,[dim1 L]);

M = 100*ones(N, 1);

nn_vec = repmat([1:N], J, 1); % dim1 x 1 (integers between 1 and N)
nn_vec = reshape(nn_vec, [dim1 1]); % This is full cartesian product

M_rep = M(nn_vec);
beta_true = normrnd(0,1,[L 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate Y
V = X * beta_true; % dim1 x 1
Vmax = accumarray(nn_vec, V, [N 1], @max); % N x 1
V = V - Vmax(nn_vec); % dim1 x 1
tmp = log(accumarray(nn_vec, exp(V), [N 1])); % N x 1
logp = V - tmp(nn_vec); % dim1 x 1
p = exp(logp); % dim1 x 1

Y = zeros(dim1, 1);
for nn = 1:N
	idces_nn = find(nn_vec == nn); % Jn x 1
	p_nn = p(idces_nn); % Jn x 1
	M_nn = M(nn); % scalar
	Y(idces_nn) = mnrnd(M_nn, p_nn); % Jn x 1
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Estimate model

LL_constant = sum(gammaln(M+1)) - sum(gammaln(Y+1));

addpath('../../../Optimization/Matlab');
obj = @(beta_val) loglikelihood_MNL(Y, X, M, M_rep, nn_vec, beta_val, LL_constant, true);
beta0 = zeros(L, 1);
tic
[beta_star, LL_star, LL_grad, FisherInfo] = Newton_Raphson(obj, beta0, 1e-6, 2000);
toc

% Compute std. errors
beta_ses = sqrt(diag(inv(FisherInfo)));

disp(table(beta_true, beta_star, beta_ses));
