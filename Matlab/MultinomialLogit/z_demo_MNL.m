addpath('../../../Optimization/Matlab');

L = 15;
N = 4000;
J = 80;

X = normrnd(0,1,[J*N L]);

M = 100*ones(N, 1);

nn_vec = repmat([1:N], J, 1); % This is full cartesian product
nn_vec = reshape(nn_vec, [J*N 1]); % (J*N) x 1 (integers between 1 and N)

beta_true = normrnd(0,1,[L 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate Y
V = X * beta_true; % (J*N) x 1
Vmax = accumarray(nn_vec, V, [N 1], @max); % N x 1
V = V - Vmax(nn_vec); % dim1 x 1
tmp = log(accumarray(nn_vec, exp(V), [N 1])); % N x 1
logp = V - tmp(nn_vec); % dim1 x 1
p = exp(logp); % (J*N) x 1

Y = zeros(J*N, 1);
for nn = 1:N
	idces_nn = find(nn_vec == nn); % Jn x 1
	p_nn = p(idces_nn); % Jn x 1
	M_nn = M(nn); % scalar
	Y(idces_nn) = mnrnd(M_nn, p_nn); % Jn x 1
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Estimate model
[beta_star, LL_star, LL_grad, FisherInfo, beta_ses] = estimate_MNL(M, X, Y, nn_vec);

disp(table(beta_true, beta_star, beta_ses));
