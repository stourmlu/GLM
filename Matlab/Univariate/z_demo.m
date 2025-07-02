rng(42); % For reproducibility

%%% Add paths
addpath('../../../Optimization/Matlab');
addpath('efficient_computation');
addpath('efficient_computation/dim1');
addpath('efficient_computation/dim2');
addpath('efficient_computation/dim3');
addpath('efficient_computation/general');

%%% GENERATE X
dimType = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% dim1 %%%%%
if dimType == 1
	NumObs = 1000;
	Xparts = cell(1,1);
	Xparts{1}.X            = [ones(NumObs, 1) normrnd(0,1,[NumObs, 9])]; % NumObs x NumX
	Xparts{1}.X_FEs        = zeros(NumObs, 2);                           % NumObs x NumX_FEs
	Xparts{1}.NumX_FE_vals = [3 4];                                		 % 1 x NumX_FEs
	Xparts{1}.X_FEs(:,1)   = mnrnd(1,repmat([0.4 0.3 0.3], NumObs, 1)) * [1:3]';
	Xparts{1}.X_FEs(:,2)   = mnrnd(1,repmat([0.3 0.3 0.2 0.2], NumObs, 1)) * [1:4]';
	
	data.Xparts = Xparts;
	data.dims = make_dims(data, 1);
	
	
	%%%%%
	%% For testing purposes: expand everything into one matrix X
	X = zeros(data.dims.NumObs, 0);
	for pp = 1:length(data.Xparts)
		Xpp = data.Xparts{pp}.X;
		X_FE_pp = data.Xparts{pp}.X_FEs;
		dim1 = size(Xpp,1);
		for ee = 1:size(X_FE_pp, 2)
			NumVals = Xparts{pp}.NumX_FE_vals(ee);
			Xdummies_ee = zeros(dim1, NumVals);
			for ff = 1:NumVals
				Xdummies_ee(:,ff) = double(X_FE_pp(:,ee) == ff);
			end
			Xdummies_ee = Xdummies_ee(:,2:end);
			Xpp = [Xpp Xdummies_ee];
		end
		X = [X Xpp];
	end
	data2.X = X;
	clear X Xpp X_FE_pp Xdummies_ee;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% dim2 %%%%%
if dimType == 2
	T = 100;
	K = 100;
	NumObs = T*K;
	
	Xparts = cell(3,1);
	Xparts{1}.X = normrnd(0,1,[T*K, 3]);    % (T*K) x 3
%%%
%	Xparts{1}.X_FEs        = zeros(T*K, 0); % (T*K) x NumX_FEs
%	Xparts{1}.NumX_FE_vals = [];    % 1 x NumX_FEs
	%%%
	Xparts{1}.X_FEs        = zeros(T*K, 3); % (T*K) x NumX_FEs
	Xparts{1}.NumX_FE_vals = [2 3 4];    % 1 x NumX_FEs
	Xparts{1}.X_FEs(:,1)   = mnrnd(1,repmat([0.6 0.4],         T*K, 1)) * [1:2]';
	Xparts{1}.X_FEs(:,2)   = mnrnd(1,repmat([0.4 0.3 0.3],     T*K, 1)) * [1:3]';
	Xparts{1}.X_FEs(:,3)   = mnrnd(1,repmat([0.3 0.3 0.2 0.2], T*K, 1)) * [1:4]';
%%%
	
	Xparts{2}.X            = [ones(T, 1) normrnd(0,1,[T, 3])]; % T x 4
	Xparts{2}.X_FEs        = zeros(T, 2);                      % T x NumX_FEs
	Xparts{2}.NumX_FE_vals = [2 5];                       % 1 x NumX_FEs
	Xparts{2}.X_FEs(:,1)   = mnrnd(1,repmat([0.5 0.5],             T, 1)) * [1:2]';
	Xparts{2}.X_FEs(:,2)   = mnrnd(1,repmat([0.2 0.2 0.2 0.2 0.2], T, 1)) * [1:5]';
	
	Xparts{3}.X = normrnd(0,1,[K, 4]);    % K x 4
	Xparts{3}.X_FEs        = zeros(K, 1); % K x NumX_FEs
	Xparts{3}.NumX_FE_vals = [10];  % 1 x NumX_FEs
	Xparts{3}.X_FEs(:,1)   = mnrnd(1,repmat([0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1], K, 1)) * [1:10]';
	
	data.Xparts = Xparts;
	data.dims = make_dims(data, 2);
	
	%%%%%
	%% For testing purposes: expand everything into one matrix X
	X = zeros(data.dims.NumObs, 0);
	for pp = 1:length(data.Xparts)
		Xpp = data.Xparts{pp}.X;
		X_FE_pp = data.Xparts{pp}.X_FEs;
		dim1 = size(Xpp,1);
		for ee = 1:size(X_FE_pp, 2)
			NumVals = Xparts{pp}.NumX_FE_vals(ee);
			Xdummies_ee = zeros(dim1, NumVals);
			for ff = 1:NumVals
				Xdummies_ee(:,ff) = double(X_FE_pp(:,ee) == ff);
			end
			Xdummies_ee = Xdummies_ee(:,2:end);
			Xpp = [Xpp Xdummies_ee];
		end
		% Resize/repeat Xpp if necessary (depending on part pp)
		if pp == 1
			1;
		else if pp == 2
			Xpp = reshape(Xpp, [T 1 size(Xpp,2)]); Xpp = repmat(Xpp, 1, K, 1); Xpp = reshape(Xpp,[T*K, size(Xpp,3)]);
		else if pp == 3
			Xpp = reshape(Xpp, [1 K size(Xpp,2)]); Xpp = repmat(Xpp, T, 1, 1); Xpp = reshape(Xpp,[T*K, size(Xpp,3)]);
		end; end; end;
		X = [X Xpp];
	end
	data2.X = X;
	clear X Xpp X_FE_pp Xdummies_ee;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% dim3 %%%%%
if dimType == 3
	T = 40;
	K1 = 45;
	K2 = 50;
	NumObs = T*K1*K2;
	
	Xparts = cell(7,1);
	Xparts{1}.X            = [ones(T*K1*K2, 1) normrnd(0,1,[T*K1*K2, 3])]; % (T*K1*K2) x 3
	Xparts{1}.X_FEs        = zeros(T*K1*K2, 2);                            % (T*K1*K2) x NumX_FEs
	Xparts{1}.NumX_FE_vals = [2 3];                                        % 1 x NumX_FEs
	Xparts{1}.X_FEs(:,1)   = mnrnd(1,repmat([0.6 0.4],     T*K1*K2, 1)) * [1:2]';
	Xparts{1}.X_FEs(:,2)   = mnrnd(1,repmat([0.5 0.3 0.2], T*K1*K2, 1)) * [1:3]';
	
	Xparts{2}.X            = [normrnd(0,1,[T*K1, 4])];                     % (T*K1)    x 4
	Xparts{2}.X_FEs        = zeros(T*K1, 3);                               % (T*K1) x NumX_FEs
	Xparts{2}.NumX_FE_vals = [2 3 4];                                        % 1 x NumX_FEs
	Xparts{2}.X_FEs(:,1)   = mnrnd(1,repmat([0.6 0.4],     T*K1, 1)) * [1:2]';
	Xparts{2}.X_FEs(:,2)   = mnrnd(1,repmat([0.5 0.3 0.2], T*K1, 1)) * [1:3]';
	Xparts{2}.X_FEs(:,3)   = mnrnd(1,repmat([0.25 0.25 0.25 0.25], T*K1, 1)) * [1:4]';
	
	Xparts{3}.X            = normrnd(0,1,[T*K2, 5]);                       % (T*K2)    x 4
	Xparts{3}.X_FEs        = zeros(T*K2, 2);                               % (T*K2) x NumX_FEs
	Xparts{3}.NumX_FE_vals = [5 2];                                        % 1 x NumX_FEs
	Xparts{3}.X_FEs(:,1)   = mnrnd(1,repmat([0.2 0.2 0.2 0.2 0.2],     T*K2, 1)) * [1:5]';
	Xparts{3}.X_FEs(:,2)   = mnrnd(1,repmat([0.5 0.5], T*K2, 1)) * [1:2]';
	
	Xparts{4}.X            = normrnd(0,1,[K1*K2, 6]);                      % (K1*K2)   x 4
	Xparts{4}.X_FEs        = zeros(K1*K2, 2);                               % (K1*K2) x NumX_FEs
	Xparts{4}.NumX_FE_vals = [4 5];                                         % 1 x NumX_FEs
	Xparts{4}.X_FEs(:,1)   = mnrnd(1,repmat([0.25 0.25 0.25 0.25],     K1*K2, 1)) * [1:4]';
	Xparts{4}.X_FEs(:,2)   = mnrnd(1,repmat([0.2 0.2 0.2 0.2 0.2], K1*K2, 1)) * [1:5]';
	
	Xparts{5}.X            = normrnd(0,1,[T, 7]);                          % T         x 4
	Xparts{5}.X_FEs        = zeros(T, 1);                                  % T x NumX_FEs
	Xparts{5}.NumX_FE_vals = [2];                                        % 1 x NumX_FEs
	Xparts{5}.X_FEs(:,1)   = mnrnd(1,repmat([0.6 0.4],     T, 1)) * [1:2]';
	
	Xparts{6}.X            = normrnd(0,1,[K1, 8]);                         % K1        x 4
	Xparts{6}.X_FEs        = zeros(K1, 1);                                  % K1 x NumX_FEs
	Xparts{6}.NumX_FE_vals = [3];                                         % 1 x NumX_FEs
	Xparts{6}.X_FEs(:,1)   = mnrnd(1,repmat([0.5 0.3 0.2], K1, 1)) * [1:3]';
	
	
	Xparts{7}.X            = normrnd(0,1,[K2, 9]);                         % K2        x 4
	Xparts{7}.X_FEs        = zeros(K2, 2);                                  % K2 x NumX_FEs
	Xparts{7}.NumX_FE_vals = [4 4];                                         % 1 x NumX_FEs
	Xparts{7}.X_FEs(:,1)   = mnrnd(1,repmat([0.25 0.25 0.25 0.25],     K2, 1)) * [1:4]';
	Xparts{7}.X_FEs(:,2)   = mnrnd(1,repmat([0.5 0.2 0.2 0.1], K2, 1)) * [1:4]';
	
	data.Xparts = Xparts;
	data.dims = make_dims(data, 3);
	
	%%%%%%
	%% For testing purposes: expand everything into one matrix X
	X = zeros(data.dims.NumObs, 0);
	for pp = 1:length(data.Xparts)
		Xpp = data.Xparts{pp}.X;
		X_FE_pp = data.Xparts{pp}.X_FEs;
		dim1 = size(Xpp,1);
		for ee = 1:size(X_FE_pp, 2)
			NumVals = Xparts{pp}.NumX_FE_vals(ee);
			Xdummies_ee = zeros(dim1, NumVals);
			for ff = 1:NumVals
				Xdummies_ee(:,ff) = double(X_FE_pp(:,ee) == ff);
			end
			Xdummies_ee = Xdummies_ee(:,2:end);
			Xpp = [Xpp Xdummies_ee];
		end
		% Resize/repeat Xpp if necessary (depending on part pp)
		if pp == 1
			1;
		else if pp == 2
			Xpp = reshape(Xpp, [T K1 1 size(Xpp,2)]); Xpp = repmat(Xpp, 1, 1, K2, 1); Xpp = reshape(Xpp,[T*K1*K2, size(Xpp,4)]);
		else if pp == 3
			Xpp = reshape(Xpp, [T 1 K2 size(Xpp,2)]); Xpp = repmat(Xpp, 1, K1, 1, 1); Xpp = reshape(Xpp,[T*K1*K2, size(Xpp,4)]);
		else if pp == 4
			Xpp = reshape(Xpp, [1 K1 K2 size(Xpp,2)]); Xpp = repmat(Xpp, T, 1, 1, 1); Xpp = reshape(Xpp,[T*K1*K2,  size(Xpp,4)]);
		else if pp == 5
			Xpp = reshape(Xpp, [T 1 1 size(Xpp,2)]); Xpp = repmat(Xpp, 1, K1, K2, 1); Xpp = reshape(Xpp,[T*K1*K2, size(Xpp,4)]);
		else if pp == 6
			Xpp = reshape(Xpp, [1 K1 1 size(Xpp,2)]); Xpp = repmat(Xpp, T, 1, K2, 1); Xpp = reshape(Xpp,[T*K1*K2, size(Xpp,4)]);
		else if pp == 7
			Xpp = reshape(Xpp, [1 1 K2 size(Xpp,2)]); Xpp = repmat(Xpp, T, K1, 1, 1); Xpp = reshape(Xpp,[T*K1*K2, size(Xpp,4)]);
		end;end;end;end;end;end;end;
		X = [X Xpp];
	end
	data2.X = X;
	clear X Xpp X_FE_pp Xdummies_ee;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% general %%%%%
if dimType == 0
	NumObs = 10000;
	
	Xparts = cell(3,1);
	NumObs_pp = 200;
	Xparts{1}.X            = [ones(NumObs_pp, 1) normrnd(0,1,[NumObs_pp, 3])]; % NumObs_pp x 4
	Xparts{1}.X_FEs        = zeros(NumObs_pp, 2);                            % NumObs_pp x NumX_FEs
	Xparts{1}.NumX_FE_vals = [2 3];                                        % 1 x NumX_FEs
	Xparts{1}.X_FEs(:,1)   = mnrnd(1,repmat([0.6 0.4],     NumObs_pp, 1)) * [1:2]';
	Xparts{1}.X_FEs(:,2)   = mnrnd(1,repmat([0.5 0.3 0.2], NumObs_pp, 1)) * [1:3]';
	Xparts{1}.mapping = mnrnd(1, repmat(1/NumObs_pp*ones(1,NumObs_pp), NumObs, 1)) * [1:NumObs_pp]';
	
	
	
	NumObs_pp = 300;
	Xparts{2}.X            = [normrnd(0,1,[NumObs_pp, 4])];                     % NumObs_pp    x 4
	Xparts{2}.X_FEs        = zeros(NumObs_pp, 3);                               % NumObs_pp x NumX_FEs
	Xparts{2}.NumX_FE_vals = [2 3 4];                                        % 1 x NumX_FEs
	Xparts{2}.X_FEs(:,1)   = mnrnd(1,repmat([0.6 0.4],     NumObs_pp, 1)) * [1:2]';
	Xparts{2}.X_FEs(:,2)   = mnrnd(1,repmat([0.5 0.3 0.2], NumObs_pp, 1)) * [1:3]';
	Xparts{2}.X_FEs(:,3)   = mnrnd(1,repmat([0.25 0.25 0.25 0.25], NumObs_pp, 1)) * [1:4]';
	Xparts{2}.mapping = mnrnd(1, repmat(1/NumObs_pp*ones(1,NumObs_pp), NumObs, 1)) * [1:NumObs_pp]';
	
	NumObs_pp = 500;
	Xparts{3}.X            = normrnd(0,1,[NumObs_pp, 5]);                       % NumObs_pp    x 4
	Xparts{3}.X_FEs        = zeros(NumObs_pp, 2);                               % NumObs_pp x NumX_FEs
	Xparts{3}.NumX_FE_vals = [5 2];                                        % 1 x NumX_FEs
	Xparts{3}.X_FEs(:,1)   = mnrnd(1,repmat([0.2 0.2 0.2 0.2 0.2],     NumObs_pp, 1)) * [1:5]';
	Xparts{3}.X_FEs(:,2)   = mnrnd(1,repmat([0.5 0.5], NumObs_pp, 1)) * [1:2]';
	Xparts{3}.mapping = mnrnd(1, repmat(1/NumObs_pp*ones(1,NumObs_pp), NumObs, 1)) * [1:NumObs_pp]';
	
	data.Xparts = Xparts;
	data.dims = make_dims(data, 0);
	
	%%%%%%
	%% For testing purposes: expand everything into one matrix X
	X = zeros(data.dims.NumObs, 0);
	for pp = 1:length(data.Xparts)
		Xpp = data.Xparts{pp}.X;
		X_FE_pp = data.Xparts{pp}.X_FEs;
		dim1 = size(Xpp,1);
		for ee = 1:size(X_FE_pp, 2)
			NumVals = Xparts{pp}.NumX_FE_vals(ee);
			Xdummies_ee = zeros(dim1, NumVals);
			for ff = 1:NumVals
				Xdummies_ee(:,ff) = double(X_FE_pp(:,ee) == ff);
			end
			Xdummies_ee = Xdummies_ee(:,2:end);
			Xpp = [Xpp Xdummies_ee];
		end
		% Resize/repeat Xpp if necessary (depending on part pp)
		Xpp = Xpp(data.Xparts{pp}.mapping,:);
		
		X = [X Xpp];
	end
	data2.X = X;
	clear X Xpp X_FE_pp Xdummies_ee;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Generate beta
beta_true = [-7; normrnd(0,1,[data.dims.NumParams-1 1])];

%%% Compute V
V = compute_Xbeta(data, beta_true);


%%%%% LINEAR MODEL
%data.Y = V + normrnd(0,0.5,[NumObs,1]);
%tic
%[beta_star] = estimate_OLS(data, true);
%toc
%
%data2.Y = data.Y;
%tic
%[beta_star2] = estimate_OLS(data2);
%toc
%table(beta_true, beta_star, beta_star2)

%%%%% POISSON
%logLambdaOffset = -8 + 0.5*normrnd(0,1, [NumObs,1]);
%lambda = exp(logLambdaOffset + V);
%Y = poissrnd(lambda);
%if max(Y) > 1e10
%disp(max(Y));
%error('Large values')
%end
%data.Y = Y;
%[beta_star, LL_star, LL_grad, FisherInfo, beta_ses] = estimate_GLM(data, 'Poisson', logLambdaOffset);


%%%%% BINOMIAL LOGISTIC REGRESSION
NumTries = 1000*ones(NumObs,1);
p = 1./(1+exp(-V));
Y = mnrnd(NumTries, [p 1-p]); Y = Y(:,1);
data.Y = Y;
tic
[beta_star, LL_star, LL_grad, FisherInfo, beta_ses] = estimate_GLM(data, 'BinomialLogistic', NumTries);
toc

data2.Y = Y;
tic
[beta_star2, LL_star2, LL_grad2, FisherInfo2, beta_ses2] = estimate_GLM(data2, 'BinomialLogistic', NumTries);
toc

%%%% DISPLAY RESULTS
table(beta_true, beta_star, beta_ses)
table(beta_true, beta_star, beta_star2)
max(abs(beta_true - beta_star))
max(abs(beta_star - beta_star2))
