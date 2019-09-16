function inp       	= dpm_get_empty_inp(grd,dis,options)
if strcmp(options,'nan')
    value = nan;
elseif strcmp(options,'zero')
    value = 0;
elseif strcmp(options,'inf')
    value = inf;
end

for i=1:length(grd.Nx)
    inp.X{i} = value;
end
for i=1:length(grd.Nu)
    inp.U{i} = value;
end
for i=1:length(dis.W)
    inp.W{i} = value;
end
inp.Ts = dis.Ts;
