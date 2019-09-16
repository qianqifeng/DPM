function grd        = input_check_grd(grd,T)%检测grd的状态和输入变量，并扩展为T+1个

% grd = 
%     Xn: {[1x1 struct]}
%     XN: {[1x1 struct]}
%     X0: {[250]}
%     Nx: {[1x1001 double]}
%     Nu: {[1x1001 double]}
%     Un: {[1x1 struct]}
if T < 1
        error('DPM:Internal','prb.N must be greater than 0')
end
for i=1:length(grd.Nx)
    if grd.Nx{i} < 1
        error('DPM:Internal','grd.Nx{.} must be equal or greater than 1')
    end
    if length(grd.Xn{i}.lo)==1
        grd.Xn{i}.lo = repmat(grd.Xn{i}.lo,1,T+1);
    elseif length(grd.Xn{i}.lo)~=T+1
        error('DPM:Internal','grd.Xn{.}.lo must be a scalar OR have the same length as the problem')
    end
    if length(grd.Xn{i}.hi)==1
        grd.Xn{i}.hi = repmat(grd.Xn{i}.hi,1,T+1);
    elseif length(grd.Xn{i}.hi)~=T+1
        error('DPM:Internal','grd.Xn{.}.hi must be a scalar OR have the same length as the problem')
    end
    if length(grd.Nx{i})==1
        grd.Nx{i} = repmat(grd.Nx{i},1,T+1);
    elseif length(grd.Nx{i})~=T+1
        error('DPM:Internal','grd.Nx{.} must be a scalar OR have the same length as the problem')
    end
end%检测grd的状态变量，并扩展为N+1个
for i=1:length(grd.Nu)
    if grd.Nu{i} < 1
        error('DPM:Internal','grd.Nu{.} must be equal or greater than 1')
    end
    if length(grd.Un{i}.lo)==1
        grd.Un{i}.lo = repmat(grd.Un{i}.lo,1,T);
    elseif length(grd.Un{i}.lo)~=T
        error('DPM:Internal','grd.Un{.}.lo must be a scalar OR have the same length as the problem')
    end
    if length(grd.Un{i}.hi)==1
        grd.Un{i}.hi = repmat(grd.Un{i}.hi,1,T);
    elseif length(grd.Un{i}.hi)~=T
        error('DPM:Internal','grd.Un{.}.hi must be a scalar OR have the same length as the problem')
    end
    if length(grd.Nu{i})==1
        grd.Nu{i} = repmat(grd.Nu{i},1,T);
    elseif length(grd.Nu{i})~=T
        error('DPM:Internal','grd.Nu{.} must be a scalar OR have the same length as the problem')
    end
end%检测grd的输入变量，并扩展为N+1个
