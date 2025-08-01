function [dims] = make_dims(data, dimType)
	
	
	
	
	% Make dims object (specific)
	dims.dimType = dimType;
	
	% Make dims object (general)
	dims.NumParts = length(data.Xparts);
	dims.dims1               = cell(1,dims.NumParts);
	dims.Xpart_2_NumX        = zeros(1,dims.NumParts);
	dims.Xpart_2_NumX_FEs    = zeros(1,dims.NumParts);
	dims.Xpart_2_Num_FE_vals = cell(1,dims.NumParts);
	dims.NumFEvals2Keep      = cell(1,dims.NumParts);
	NumParams                = 0;
	for pp = 1:dims.NumParts
		dims.dims1{pp}               = size(data.Xparts{pp}.X,1);
		dims.Xpart_2_NumX(pp)        = size(data.Xparts{pp}.X,2);
		dims.Xpart_2_NumX_FEs(pp)    = size(data.Xparts{pp}.X_FEs,2);
		dims.Xpart_2_Num_FE_vals{pp} = data.Xparts{pp}.NumX_FE_vals;
		dims.NumFEvals2Keep{pp}      = sum(dims.Xpart_2_Num_FE_vals{pp}) - dims.Xpart_2_NumX_FEs(pp);
		NumParams                    = NumParams + dims.Xpart_2_NumX(pp) + dims.NumFEvals2Keep{pp};
	end
	if dimType == 1 || dimType == 2 || dimType == 3
		dims.NumObs = dims.dims1{1};
	end
	if dimType == 0
		dims.NumObs = size(data.Xparts{1}.mapping, 1);
	end
	dims.NumParams = NumParams;
	
	if dimType == 0
		dims.mappings = cell(dims.NumParts,1);
		for pp = 1:dims.NumParts
			dims.mappings{pp} = data.Xparts{pp}.mapping;
		end
	end
end