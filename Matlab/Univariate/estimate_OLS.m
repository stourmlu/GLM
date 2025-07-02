function [beta_hat, varargout] = estimate_OLS(data, varargin)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Inputs:
	%	data
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%% Outputs:
	%	beta_hat:			NumX x 1
	% varargout:
	%   {1}=sigmasq_hat:	scalar
	%   {2}=Rsq:			scalar
	%   {3}=std_errors:		NumX x 1
	%   {4}=t_values:		NumX x 1
	%   {5}=p_values:		NumX x 1
	%   {6}=F_stat:			scalar
	%   {7}=residuals:		NumObs x 1
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	if length(varargin) >= 1
		verbose = varargin{1};
	else
		verbose = true;
	end
	
	g = compute_Xt_Y(data, data.Y)';  % NumX x 1
	
	if isfield(data, 'Xparts')
		NumObs = data.dims.NumObs;
	else
		NumObs = size(data.X,1);
	end
	H = compute_Xt_X_Y(data, ones(NumObs,1)); % NumX x NumX
	
	beta_hat = H\g;
	
	
	if verbose || nargout >= 2 % Calculate sigmasq_hat: estimate of sigma^2
		NumX = size(beta_hat,1);
		residuals = data.Y - compute_Xbeta(data, beta_hat);
		sigmasq_hat = sum(residuals.^2)/(NumObs - NumX);
		varargout{1} = sigmasq_hat;
	end
	
	if verbose || nargout >= 3 % Compute R^2
		SSE = sum(residuals.^2);
		SST = sum((data.Y - mean(data.Y)).^2); % TO DO
		SSR = SST - SSE;
		Rsq = 1 - SSE/SST;
		varargout{2} = Rsq;
	end

	if verbose || nargout >= 4 % Obtain standard errors
		std_errors = sqrt(diag(sigmasq_hat * inv(H)));
		varargout{3} = std_errors;
	end

	if verbose || nargout >= 5 % Obtain t-values
		t_values = beta_hat./std_errors;
		varargout{4} = t_values;
	end

	if verbose || nargout >= 6 % Obtain p-values
		p_values = 2*(1 - tcdf(abs(t_values), NumObs-NumX));
		varargout{5} = p_values;
	end
	
	if verbose || nargout >= 7 % Obtain F statistic
		F_stat = (SSR/(NumX-1))/(SSE/(NumObs - NumX));
		varargout{6} = F_stat;
	end
	
	if verbose || nargout >= 8 % Output residuals
		varargout{7} = residuals;
	end
	
	if verbose
		mytable = table(beta_hat,std_errors, t_values, p_values);
		if isfield(data, 'Xnames')
			mytable.Properties.RowNames = data.Xnames;
		end
		disp(mytable);
	end
end
