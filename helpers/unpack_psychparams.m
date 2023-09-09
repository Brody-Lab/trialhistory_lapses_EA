function[prm_struct] = unpack_psychparams(params, varargin)


if isempty(varargin)
    variable = 'beta';
else
    variable = varargin{1};
end


if params.trust == 1
    prm_struct.left_lapse = params.(variable)(1);
    prm_struct.right_lapse = 1 - params.(variable)(2) - params.(variable)(1);
    prm_struct.sens = params.(variable)(3);
    prm_struct.bias = params.(variable)(4);
else
    prm_struct.left_lapse = nan;
    prm_struct.right_lapse = nan;
    prm_struct.sens = nan;
    prm_struct.bias = nan;
end
    
end
