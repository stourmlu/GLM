function [Xnames] = make_Xnames(data)
	
	Xnames = cell(data.dims.NumParams, 1);
	idx = 0;
	for pp = 1:length(data.Xparts)
		myNumParams = size(data.Xparts{pp}.X, 2);
		if isfield(data.Xparts{pp}, 'Xnames')
			myParamNames = data.Xparts{pp}.Xnames;
		else
			myParamNames = strcat(sprintf('X%d', pp), sprintfc('_%d', 1:myNumParams));
		end
		Xnames(idx+1:idx+myNumParams) = myParamNames;
		idx = idx+myNumParams;
		
		for ee = 1:size(data.Xparts{pp}.X_FEs, 2)
			myNumParams = data.Xparts{pp}.NumX_FE_vals(ee) - 1; % (1st value has no coef in dummy coding)
			
			if isfield(data.Xparts{pp}, 'Xnames')
				myParamNames = data.Xparts{pp}.Xnames;
			else
				myParamNames = strcat(sprintf('X%d_FE_%d', pp, ee), sprintfc('_value_%d', 2:myNumParams+1));
			end
			
			Xnames(idx+1:idx+myNumParams) = myParamNames;
			idx = idx + myNumParams;
		end	
	end
end
