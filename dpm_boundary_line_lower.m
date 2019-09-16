function LineLower  = dpm_boundary_line_lower(model,par,grd,dis,options)

optb               = options;
parb               = par;
grdb               = grd;
grdb.Xn_true{1}.lo = grd.Xn{1}.lo;
grdb.Xn_true{1}.hi = grd.Xn{1}.hi;
grdb.Nx{1}         = ones(1,dis.N+1);
grdb.Xn{1}.lo      = grd.X0{1}.*ones(1,dis.N+1);
grdb.Xn{1}.hi      = grd.X0{1}.*ones(1,dis.N+1);
optb.UseLine       = 0;
optb.SaveMap       = 1;
parb.model         = model;
parb.options.Iter  = options.Iter;
parb.options.Tol   = options.Tol;
optb.CalcLine      = 1;

if isfield(options,'gN')
    warning('options.gN can only be used w/ 1 state boundary')
    optb.gN{1} = grd.XN{1}.lo;%.*sub_x_ones;
    optb.gN{2} = dpm_interpn(grd.X{1},options.gN{1},grdb.X{1});%.*sub_x_ones);
else
    optb.gN{1} = grd.XN{1}.lo;%.*sub_x_ones;
    optb.gN{2} = zeros(size(grd.XN{1}.lo));%zeros(size(sub_x_ones));
end

optb.MyInf        = options.MyInf;
optb.Minimize     = 1; % PEL:2012 - this used to be =options.Minimize; 
optb.Warnings     = 'off';
dynb              = dpm(@dpm_model_inv,parb,grdb,dis,optb);
% convert cellarray to vectors
vsize = size(dynb.Jo);
for j=1:vsize(1)
    V_new{j} = nan(1,length(dynb.Jo{j,1}),vsize(2));
    for i=1:vsize(2)
        V_new{j}(1,:,i) = [dynb.Jo{j,i}];
    end
end
usize = size(dynb.Uo);
for j=1:usize(1)
    O_new{j} = nan(1,length(dynb.Uo{j,1}),usize(2));
    for i=1:usize(2)
%         if i<usize(2)
            O_new{j}(1,:,i) = [dynb.Uo{j,i}];
%         end
    end
end
dynb.Jo = V_new;
dynb.Uo = O_new;
% convert infeasible points to minimum state boundary
for i=1:dis.N
    dynb.Jo{1}(1,dynb.Jo{1}(1,:,i)>grd.Xn{1}.hi(i),i) = grd.Xn{1}.lo(i);
    dynb.Jo{1}(1,isnan(dynb.Jo{1}(1,:,i)),i)          = grd.Xn{1}.lo(i).*ones(size(dynb.Jo{1}(1,isnan(dynb.Jo{1}(1,:,i)),i)));
end
dynb.Jo{2}(isnan(dynb.Jo{2}) | dynb.Jo{2}>=optb.MyInf) = options.MyInf;
for i=1:length(dynb.Uo)
    dynb.Uo{i}(isnan(dynb.Uo{i})) = 1;
end
% insert lower boundary into original problem definition
LineLower.Xo       = dynb.Jo{1};
LineLower.Uo       = dynb.Uo;
LineLower.Jo       = dynb.Jo{2}; % add final cost term
