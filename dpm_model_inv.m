function [X C I out]= dpm_model_inv(inp,par)
inpo = inp;
iterations = 0;
dSOC = inf;

while max(abs(reshape(dSOC,1,numel(dSOC)))) > par.options.Tol && iterations < par.options.Iter
%     [out.X out.C out.I] = eval([par.model '(inp,par);']);
    [X C I] = feval(par.model,inp,par);
    dSOC   = X{1} - inpo.X{1};
    inp.X{1} = inp.X{1} - dSOC;
    iterations = iterations+1;
end
% out.I = bitor(out.I,bitor(out.X{1}<lim.Xn{1}.lo, out.X{1}>lim.Xn{1}.hi));
X = inp.X;
C{2} = C{1};
C{1} = (X{1}-inpo.X{1});
out = [];
