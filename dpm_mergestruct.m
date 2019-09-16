function S        	= dpm_mergestruct(S1,S2)
% WARNING: might cause memory problems!
try
    S = S1;
    names = fieldnames(S1);
    for i=1:length(names)
        if isstruct(S1.(names{i}))
            S.(names{i}) = dpm_mergestruct(S1.(names{i}), S2.(names{i}));
        elseif iscell(S1.(names{i}))
            for j=1:numel(S1.(names{i}))
                S.(names{i}){j} = [S1.(names{i}){j} S2.(names{i}){j}];
            end
        elseif isnumeric(S1.(names{i})) || islogical(S1.(names{i}))
            S.(names{i}) = [S1.(names{i}) S2.(names{i})];
        else
            S.(names{i}) = S2.(names{i});
        end
    end
catch
    error('mergestruct: S1 and S2 have different structures.')
end
