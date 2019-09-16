function out       	= dpm_get_empty_out(model,inp,par,grd,options)
%GET_EMPTY_OUT   Gets an empty result struct
%   OUT   = GET_EMPTY_OUT(MODEL,OPTIONS)
%   Gets an empty output struct from model.
%
%   See also dpm_forward_sim
%
%   Author(s): Olle L. Sundström, 29.11.2006
%   Copyright 2006- Olle L. Sundström
if ~exist('options')
    options = 'nan';
end

[X C I out] = feval(model,inp,par);
out.X = X;
out.C = C;
out.I = I;
out = dpm_setallfield(out,options);

function S         	= dpm_setallfield(S1,options)

try
    if strcmp(options,'nan')
        value = nan;
    elseif strcmp(options,'zero')
        value = 0;
    elseif strcmp(options,'inf')
        value = inf;
    end

    S = S1;
    names = fieldnames(S1);
    for i=1:length(names)
        if isstruct(S1.(names{i}))
            S.(names{i}) = dpm_setallfield(S1.(names{i}), options);
        elseif iscell(S1.(names{i}))
            for j=1:numel(S1.(names{i}))
                S.(names{i}){j} = value;
            end
        elseif isnumeric(S1.(names{i})) || islogical(S1.(names{i}))
            S.(names{i}) = value;
        else
            S.(names{i}) = S1.(names{i});
        end
    end
catch
    error('mergestruct: S1 and S2 have different structures.')
end
