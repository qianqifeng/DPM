function dyn       	= dpm_backward(model,par,grd,dis,options)
% DP_BACKWARD   Computes the optimal input matrix and the optimal cost-to-go.
%
%    Examples:
%
%       Computes the optimal input matrix and the optimal cost-to-go
%          dyn = dpm_backward(model,par,grd,dis,options)
%
%   -THIS IS THE FIRST CORE OF THE DYNAMIC PROGRAMMING-
%
%   Author(s): Olle L. Sundstr鰉, 29.11.2006
%   Copyright 2006- Olle L. Sundstr鰉

% Calculate boundaries
%xs = 0;

% Initialize optimal input, and cost maps
%in = 0;
x_sze = ones(size(grd.Nx));
for i=1:length(grd.Nx)
    x_sze(i) = grd.Nx{i}(end);
end
u_sze = ones(size(grd.Nu));
for i=1:length(grd.Nu)
    u_sze(i) = grd.Nu{i}(end);
end

% If UseLine ______________________________________________________________
if options.UseLine
    dyn.B.lo  = dpm_boundary_line_lower(model,par,grd,dis,options);
    dyn.B.hi  = dpm_boundary_line_upper(model,par,grd,dis,options);
end        


% generate current state grid
for i=1:length(grd.Nx)
    if options.UseLine && i==1 && ~options.FixedGrid
        if length(grd.Nx)==1
            current_grd.X{i} = linspace(dyn.B.lo.Xo(end),dyn.B.hi.Xo(end),grd.Nx{i}(end))';
        else
            current_grd.X{i} = linspace(min(dyn.B.lo.Xo(:,i)),max(dyn.B.hi.Xo(:,i)),grd.Nx{i}(i))';
        end
    else
        current_grd.X{i} = linspace(grd.Xn{i}.lo(end),grd.Xn{i}.hi(end),grd.Nx{i}(end))';
        if options.CalcLine
            grd.Xn{1}.lo = grd.Xn_true{1}.lo;
            grd.Xn{1}.hi = grd.Xn_true{1}.hi;
        end
    end
end



% generate current input grid
for i=1:length(grd.Nu)
    current_grd.U{i} = linspace(grd.Un{i}.lo(end),grd.Un{i}.hi(end),grd.Nu{i}(end))';
end

inp0  = dpm_get_empty_inp(grd,dis,'zero');
out0  = dpm_get_empty_out(model,inp0,par,grd,'zero');

% Initialize the outputs dyn.Uo and dyn.Jo
for i=1:length(out0.C)
    % if dyn.Jo is specified in options
    if isfield(options,'gN') && length(options.gN)>= i
        if dpm_sizecmp(options.gN{i},zeros(get_size(current_grd)))
            dyn.Jo{i} = options.gN{i};
        else
            error('DPM:Internal',['options.gN{' num2str(i) '} has incorrect dimesions']);
        end
    % if dyn.Jo is NOT specified in options
    else
        dyn.Jo{i} = zeros(get_size(current_grd));
    end
end

% if not using boundary line set cost of states outside feasible region
% to myinf
if ~options.UseLine && ~options.CalcLine && ~options.UseLevelSet
    for i=1:length(grd.Nx)
        eval(['dyn.Jo{1}(' repmat(':,',1,i-1) 'current_grd.X{i} > grd.XN{i}.hi' repmat(',:',1,length(grd.Nx)-i) ') = options.MyInf;'])
        eval(['dyn.Jo{1}(' repmat(':,',1,i-1) 'current_grd.X{i} < grd.XN{i}.lo' repmat(',:',1,length(grd.Nx)-i) ') = options.MyInf;'])
    end
end


if options.UseLevelSet
    %set up level set function  
    dyn.Jo{end+1} = -inf*zeros(get_size(current_grd));
    %initialize by V_N = h(x_N)
    if length(grd.Nx)>1
        code_fin_cst = '[';
        for i=1:length(grd.Xn)
            code_fin_cst = [code_fin_cst, 'x{', num2str(i),'} '];
        end
        code_fin_cst = [code_fin_cst, '] = ndgrid('];
        for i=1:length(grd.Xn)
            code_fin_cst = [code_fin_cst, '(current_grd.X{', num2str(i),'}),'];
        end
        code_fin_cst = code_fin_cst(1:end-1);
        code_fin_cst = [code_fin_cst, ');'];
    else
        code_fin_cst = 'x{1} = current_grd.X{1};';
    end
    eval(code_fin_cst);
    for i=1:length(x)
        dyn.Jo{end} = max(dyn.Jo{end},max(grd.XN{i}.lo-x{i},x{i}-grd.XN{i}.hi)); 
    end
    if length(grd.Nx)>1
        code_x_grd = '[';
        for i=1:length(grd.Xn)
            code_x_grd = [code_x_grd, 'x', num2str(i),' '];
        end
        code_x_grd = [code_x_grd, '] = ndgrid('];
        for i=1:length(grd.Xn)
            code_x_grd = [code_x_grd, '(1:x_sze(', num2str(i),'))'','];
        end
        code_x_grd = code_x_grd(1:end-1);
        code_x_grd = [code_x_grd, ');'];
    else
        code_x_grd = 'x1 = (1:x_sze(1))'';';
    end
end

% initialization if the entire cost-to-go map should be saved 初始化
if options.SaveMap
    V_map = cell(length(dyn.Jo),dis.N+1);
    for i=1:length(dyn.Jo)
        for j=1:dis.N+1
            x_sze_n = [];
            for k=1:length(grd.Nx)
                x_sze_n = [x_sze_n grd.Nx{k}(j)];
            end
            if length(x_sze_n)==1
                x_sze_n = [x_sze_n 1];
            end
            V_map{i,j}          = nan(x_sze_n);
        end
        V_map{i,dis.N+1} = dyn.Jo{i};
    end
end

% initialize the optimal input map(s)初始化
dyn.Uo = cell(length(grd.Nu),dis.N);
for i=1:length(grd.Nu)
    for j=1:dis.N
        x_sze_n = [];
        for k=1:length(grd.Nx)
            x_sze_n = [x_sze_n grd.Nx{k}(j)];
        end
        if length(x_sze_n)==1
            x_sze_n = [x_sze_n 1];
        end
        dyn.Uo{i,j} = nan.*ones(x_sze_n);
    end
end


% GENERATE CODE FOR GENERATING GRID
code_generate_grid = '[';
for i=1:length(grd.Nx)
    code_generate_grid = [code_generate_grid 'inp.X{' num2str(i) '} '];
end
for i=1:length(grd.Nu)
    code_generate_grid = [code_generate_grid 'inp.U{' num2str(i) '} '];
end
code_generate_grid = [code_generate_grid '] = ndgrid('];
for i=1:length(grd.Nx)
    code_generate_grid = [code_generate_grid 'current_grd.X{' num2str(i) '},'];
end
for i=1:length(grd.Nu)
    code_generate_grid = [code_generate_grid 'current_grd.U{' num2str(i) '},'];
end
code_generate_grid = code_generate_grid(1:end-1);
code_generate_grid = [code_generate_grid ');'];


% GENERATE CODE FOR COST-TO-GO INTERPOLATION
eval(['xsize = [' dpm_code('length(current_grd.X{#}) ',1:length(grd.Nx)) '];']);
for k=1:length(dyn.Jo)
    code_cost_to_go_interp{k} = ['cost_to_go{' num2str(k) '} = dpm_interpn('];
    for i=fliplr(find(xsize>1))
        code_cost_to_go_interp{k} = [code_cost_to_go_interp{k} 'previous_grd.X{' num2str(i) '},'];
    end
    code_cost_to_go_interp{k} = [code_cost_to_go_interp{k} 'dyn.Jo{' num2str(k) '},'];
    for i=fliplr(find(xsize>1))
        if options.CalcLine
            code_cost_to_go_interp{k} = [code_cost_to_go_interp{k} 'inp.X{' num2str(i) '},'];
        else
            code_cost_to_go_interp{k} = [code_cost_to_go_interp{k} 'out.X{' num2str(i) '},'];
        end
    end
    code_cost_to_go_interp{k} = code_cost_to_go_interp{k}(1:end-1);
    code_cost_to_go_interp{k} = [code_cost_to_go_interp{k} ');'];
end
% GENERATE CODE FOR CONVERTING IND TO SUB
code_ind2str = '[';
for i=1:length(grd.Nu)
    code_ind2str = [code_ind2str 'uo' num2str(i) ' '];
end
code_ind2str = [code_ind2str '] = ind2sub(u_sze,ui);'];


% DYNAMIC PROGRAMMING
% Display progres bar of Dynamic Programming backward iteration
if ~isdeployed && strcmp(options.Waitbar,'on')
    if ~options.CalcLine
        h = waitbar(1,'DP running backwards. Please wait...');
    else
        h = waitbar(1,'DP calculating boundary line. Please wait...');
    end
    set(h,'name','DPM:Waitbar');
end

if ~isdeployed && strcmp(options.Verbose,'on')
    if ~options.CalcLine
        fprintf('%s','DP running backwards:     %%');
    else
        fprintf('%s','DP calculating boundary line:     %%');
	end
end


if options.UseLine
    if length(grd.Nx)==1
        isfeas = dyn.B.lo.Xo(:,1) > grd.X0{1} | dyn.B.hi.Xo(:,1) < grd.X0{1};
    elseif length(grd.Nx)==2
        isfeas = dpm_interpn(current_grd.X{2},dyn.B.lo.Xo(:,:,1),grd.X0{2}) > grd.X0{1} | dpm_interpn(current_grd.X{2},dyn.B.hi.Xo(:,:,1),grd.X0{2}) < grd.X0{1};
    end
    if sum(reshape(isfeas,1,numel(isfeas)))~=0
        warning('DPM:Backward','Initial value not feasible!')
    end
end

% flag for warnings (only warned one time)
iswarned = 0;

% start of dp-backward loop
% for n=dis.N-1:-1:1
n = dis.N+1;
while n > 1
    n = n-1;
    if n == 1
        debug = 1;
    end
    previous_grd = current_grd;
    x_sze = nan(1,length(grd.Nx));
    u_sze = nan(1,length(grd.Nu));
    for i=1:length(grd.Nx)
        if options.UseLine && i==1 && ~options.FixedGrid
            if length(grd.Nx)==1
                current_grd.X{i} = linspace(dyn.B.lo.Xo(n),dyn.B.hi.Xo(n),grd.Nx{i}(n))';
            else
                current_grd.X{i} = linspace(min(min(dyn.B.lo.Xo(1,:,n:n+1))),max(min(dyn.B.hi.Xo(1,:,n:n+1))),grd.Nx{i}(n))';
            end
        elseif ~options.CalcLine || i~=1
            current_grd.X{i} = linspace(grd.Xn{i}.lo(n),grd.Xn{i}.hi(n),grd.Nx{i}(n))';
        end
        x_sze(i) = grd.Nx{i}(n);
    end
    for i=1:length(grd.Nu)
        current_grd.U{i} = linspace(grd.Un{i}.lo(n),grd.Un{i}.hi(n),grd.Nu{i}(n))';
        u_sze(i) = grd.Nu{i}(n);
    end

    % generate input and state grid
    eval(code_generate_grid);
    % inp.X{1} inp.U{1}   
    if options.CalcLine
        % if calculating boundary line initialize the first input state
        eval(['inp.X{1} = repmat(dyn.Jo{1},[ones(1,length(grd.Nx)) ' dpm_code('length(current_grd.U{#}) ',1:length(grd.Nu)) ']);']);
    end

    % Call model function _________________________________________________
    % generate disturbance
    for w = 1:length(dis.W)
        inp.W{w} = dis.W{w}(n);
    end
    inp.Ts   = dis.Ts;

try
    % call model function
    if isfield(options,'Signals') && ~isempty(options.Signals)
        [out.X out.C out.I signals] = feval(model,inp,par);
    else
        [out.X out.C out.I] = feval(model,inp,par);
    end
    
    if options.UseLevelSet
        %Calculate Level set function:
        if n==dis.N %check analytically
            Vt = -inf*ones(size(out.X{1}));
            for i=1:length(grd.Xn)
                Vt = max(Vt,max(grd.XN{i}.lo-out.X{i},out.X{i}-grd.XN{i}.hi));
            end
        else %check by interpolation
            eval(code_cost_to_go_interp{end});
            Vt = cost_to_go{end};
        end
        Vt(out.I==1) = options.MyInf;
        if length(grd.Nu)>1
            Vt = reshape(Vt,[x_sze prod(u_sze)]);
        end
        [dyn.Jo{end} ub] = min(Vt,[],length(x_sze)+1);

    end
    
    % determine the arc-cost
    for i=1:length(grd.Nx)
        out.I = bitor(out.I,out.X{i}>grd.Xn{i}.hi(n+1));
        out.I = bitor(out.I,out.X{i}<grd.Xn{i}.lo(n+1));
    end
    J  = (out.I==0).*out.C{1} + out.I.* options.MyInf;
    % _____________________________________________________________________



    % Calculate cost for entire grid
    
        if options.UseLine
                if length(grd.Nx)==1
                    cost_to_go{1} = dpm_interpf1mb(previous_grd.X{1}, dyn.Jo{1}, out.X{1}, [dyn.B.lo.Xo(n+1) dyn.B.hi.Xo(n+1)],[dyn.B.lo.Jo(n+1) dyn.B.hi.Jo(n+1)],options.MyInf);
                else
                    cost_to_go{1} = dpm_interpf2mb(previous_grd.X{1},previous_grd.X{2}, dyn.Jo{1}, out.X{1},out.X{2}, [dyn.B.lo.Xo(:,:,n+1);dyn.B.hi.Xo(:,:,n+1)],[dyn.B.lo.Jo(:,:,n+1);dyn.B.hi.Jo(:,:,n+1)],options.MyInf);
                end
        else
            if length(dyn.Jo{1})==1
                cost_to_go{1} = dyn.Jo{1};
            else
                eval(code_cost_to_go_interp{1});
            end
        end
        % total cost = arc-cost + cost-to-go!
        Jt = J + cost_to_go{1};

        if options.Minimize
            Jt(Jt>options.MyInf) = options.MyInf;
            if options.CalcLine
                Jt(Jt<grd.Xn{1}.lo(n)) = options.MyInf;
            end
        else
            Jt(Jt<options.MyInf) = options.MyInf;
            if options.CalcLine
                Jt(Jt>grd.Xn{1}.hi(n)) = options.MyInf;
            end
        end

catch
    err = lasterror;
        if exist('out','var')
            if sum(reshape(isnan(out.I),1,numel(out.I))) > 0
                error('DPM:Internal','Make sure the model does not output NaN in the variable I')
            end
            if sum(reshape(isnan(out.C{1}),1,numel(out.C{1}))) > 0
                error('DPM:Internal','Make sure the model does not output NaN in the variable C')
            end
            if sum(reshape(isnan(out.X{1}),1,numel(out.X{1}))) > 0
                error('DPM:Internal','Make sure the model does not output NaN in the variable X')
            end
            err.message = [err.message ' Error in dpm_backward at n=' num2str(n)];
        end
        rethrow(err);
    end
    if length(grd.Nu)>1
        Jt = reshape(Jt,[x_sze prod(u_sze)]);
    end
    % minimize the cost-to-go
    if options.Minimize
        [Q ui] = min(Jt,[],length(x_sze)+1);
    else
        [Q ui] = max(Jt,[],length(x_sze)+1);
    end
    
    if options.UseLevelSet
        %handle infeasible states:
        Q_inf = J + cost_to_go{1};
        Q_inf = reshape(Q_inf,[x_sze prod(u_sze)]);
        eval(code_x_grd);
        switch length(grd.Nx)
            case 1
                Qinf = Q_inf(sub2ind([x_sze prod(u_sze)],x1,ub));
            case 2
                Qinf = Q_inf(sub2ind([x_sze prod(u_sze)],x1,x2,ub));
            case 3
                Qinf = Q_inf(sub2ind([x_sze prod(u_sze)],x1,x2,x3,ub));
            case 4
                Qinf = Q_inf(sub2ind([x_sze prod(u_sze)],x1,x2,x3,x4,ub));
        end
        Q(min(Vt,[],length(x_sze)+1)>0) = Qinf(min(Vt,[],length(x_sze)+1)>0);
        ui(min(Vt,[],length(x_sze)+1)>0) = ub(min(Vt,[],length(x_sze)+1)>0);
    end

    if ~options.CalcLine && isfield(options,'Signals') && ~isempty(options.Signals)
        for i=1:length(options.Signals)
            try
                eval(['dyn.' options.Signals{i} '{n} = signals.' options.Signals{i} '(sub2ind([x_sze u_sze],' dpm_code('(1:x_sze(#))''',1:length(x_sze)) ',ui));']);
            catch
                warning('DPM:Backward','options.signals element is not found in model output.')
            end
        end
    end
    
   
    if sum(reshape(Q,1,numel(Q))==options.MyInf) == numel(Q)
        if options.UseLine % if using line then all points can be infeasible as long as eventually one point is inside boundary lines.
            if dyn.B.hi.Jo(1,1,n)==options.MyInf && dyn.B.lo.Jo(1,1,n)==options.MyInf || sum(reshape(out.X{1}>dyn.B.lo.Xo(n+1) & out.X{1}<dyn.B.hi.Xo(n+1),1,numel(out.X{1})))>0 && sum(reshape(inp.X{1}(out.X{1}>dyn.B.lo.Xo(n+1) & out.X{1}<dyn.B.hi.Xo(n+1))>dyn.B.lo.Xo(n) & inp.X{1}(out.X{1}>dyn.B.lo.Xo(n+1) & out.X{1}<dyn.B.hi.Xo(n+1))<dyn.B.hi.Xo(n),1,numel(inp.X{1}(out.X{1}>dyn.B.lo.Xo(n+1) & out.X{1}<dyn.B.hi.Xo(n+1)))))>0
                if isfield(options,'DebugMode') && options.DebugMode
                    fprintf('DPM:Model function error \n \t Entering model function at the instance where the error occured \n\t Check if the entire grid generates infeasible solutions.\n')
                    if isa(model, 'function_handle')
                        eval(['dbstop in ' func2str(model) ' at 1']);
                    else
                        eval(['dbstop in ' model ' at 1']);
                    end
                    n = n+1;
                    continue
                end
                warning('DPM:Backward','No feasible solution Q(i,j,..) = Inf   for all i,j,...')
                break;  % 在计算上边界的第495个时进入此处
            end
        else
            if ~options.CalcLine && isfield(options,'DebugMode') && options.DebugMode
                if isa(model, 'function_handle')
                    eval(['dbstop in ' func2str(model) ' at 1']);
                else
                    eval(['dbstop in ' model ' at 1']);
                end
                n = n+1;
                continue
            end
            warning('DPM:Backward','No feasible solution Q(i,j,..) = Inf   for all i,j,...')
            break;
        end            
    end    
    

    % Update optimal cost dyn.Jo with the minimum cost-to-go Q
    if options.UseLine
        if length(grd.Nx)==1
            below  = current_grd.X{1}<dyn.B.lo.Xo(n);
            above  = current_grd.X{1}>dyn.B.hi.Xo(n);
            if (sum(below)~=0|sum(above)~=0)
                deg = 1;
            end
            inside = current_grd.X{1}>=dyn.B.lo.Xo(n) & current_grd.X{1}<=dyn.B.hi.Xo(n);
        else
            below  = current_grd.X{1}<min(dyn.B.lo.Xo(:,:,n));
            above  = current_grd.X{1}>max(dyn.B.hi.Xo(:,:,n));
            inside = current_grd.X{1}>=min(dyn.B.lo.Xo(:,:,n)) & current_grd.X{1}<=max(dyn.B.hi.Xo(:,:,n));
        end
        dyn.Jo{1} = nan(size(Q));

        eval(['dyn.Jo{1}(:' repmat(',:',1,length(grd.Nx)-1) ')     = Q;']);
        % Single State: if all points are infeasible and some points are between
        % boundaries: use interpolation between boundary-data
        if dyn.B.lo.Jo(n)<options.MyInf && dyn.B.hi.Jo(n)<options.MyInf && length(grd.Nx)==1 && sum(reshape(Q,1,numel(Q))==options.MyInf) == numel(Q) && sum(reshape(current_grd.X{1}>dyn.B.lo.Xo(n) & current_grd.X{1}<dyn.B.hi.Xo(n),1,numel(current_grd.X{1})))>0
            dyn.Jo{1}(current_grd.X{1}>dyn.B.lo.Xo(n) & current_grd.X{1}<dyn.B.hi.Xo(n)) = dpm_interpn([dyn.B.lo.Xo(n) dyn.B.hi.Xo(n)],[dyn.B.lo.Jo(n) dyn.B.hi.Jo(n)],  current_grd.X{1}(current_grd.X{1}>dyn.B.lo.Xo(n) & current_grd.X{1}<dyn.B.hi.Xo(n)));
        end

        eval(['dyn.Jo{1}(above' repmat(',:',1,length(grd.Nx)-1) ') = options.MyInf;']);
        eval(['dyn.Jo{1}(below' repmat(',:',1,length(grd.Nx)-1) ') = options.MyInf;']);
    else
        dyn.Jo{1}       = Q;
    end

    eval(code_ind2str);

    for i=2:length(out.C)
        Xi{1}=[];
        for j=2:length(grd.Nx)
            Xi{j} = reshape(out.X{j},[prod(x_sze) prod(u_sze)]);
        end
        ind2 = dpm_sub2ind([prod(x_sze) prod(u_sze)],1:prod(x_sze),ui);

        Ci  = reshape(out.C{i},[prod(x_sze) prod(u_sze)]);
        if length(grd.Nx)>1
            eval(['cost_to_go{2} = dpm_interpn(' dpm_code('current_grd.X{#},',2:length(grd.Nx)) 'dyn.Jo{i}' dpm_code(',out.X{#}',2:length(grd.Nx)) ');']);
            dyn.Jo{i} = (reshape(cost_to_go{2}(ind2),x_sze)~=options.MyInf & reshape(out.I(ind2),x_sze)==0).*(reshape(Ci(ind2),x_sze) + reshape(cost_to_go{2}(ind2),x_sze)) + (reshape(cost_to_go{2}(ind2),x_sze)==options.MyInf | reshape(out.I(ind2),x_sze)~=0).*options.MyInf;
        else
            dyn.Jo{i} = (dyn.Jo{i}~=options.MyInf & out.I(ind2)==0).*(Ci(ind2) + dyn.Jo{i}) + (dyn.Jo{i}==options.MyInf | out.I(ind2)~=0).*options.MyInf;
        end

        if options.Minimize
            dyn.Jo{i}(dyn.Jo{i}>options.MyInf) = options.MyInf;
        else
            dyn.Jo{i}(dyn.Jo{i}<options.MyInf) = options.MyInf;
        end  

    end
    % Store the optimal input that minimized the cost Jt
    if options.UseLine
        for i=1:length(grd.Nu)
            if length(grd.Nx) > 1
                eval(['dyn.Uo{i,n}(:' repmat(',:',1,length(grd.Nx)-1) ') = current_grd.U{i}(uo' num2str(i) ');']);

                [ind col] = dpm_sub2indr([grd.Nx{1}(n) grd.Nx{2}(n)],ones(1,grd.Nx{2}(n)),dpm_findl(current_grd.X{1},dyn.B.lo.Xo(:,:,n)),2);               
                [s1 s2] = ind2sub([grd.Nx{1} grd.Nx{2}],ind);
                dyn.Uo{i,n}(dpm_sub2ind(size(dyn.Uo{i,n}),s1,s2)) = dyn.B.lo.Uo{i}(dpm_sub2ind(size(dyn.B.lo.Uo{i}),ones(size(col)),col,n.*ones(size(col))));

                [ind col] = dpm_sub2indr([grd.Nx{1}(n) grd.Nx{2}(n)],dpm_findu(current_grd.X{1},dyn.B.hi.Xo(:,:,n)),grd.Nx{1}(n).*ones(1,grd.Nx{2}(n)),2);               
                [s1 s2] = ind2sub([grd.Nx{1}(n) grd.Nx{2}(n)],ind);
                dyn.Uo{i,n}(dpm_sub2ind(size(dyn.Uo{i,n}),s1,s2)) = dyn.B.hi.Uo{i}(dpm_sub2ind(size(dyn.B.lo.Uo{i}),ones(size(col)),col,n.*ones(size(col))));
            else


                    % if cost-to-go is infeasible and some grid points are
                    % still inside boundaries
                    eval(['dyn.Uo{i,n}(:) = current_grd.U{i}(uo' num2str(i) ');']);
                    if dyn.B.lo.Jo(n)<options.MyInf && dyn.B.hi.Jo(n)<options.MyInf && sum(reshape(Q,1,numel(Q))==options.MyInf) == numel(Q) && sum(reshape(current_grd.X{1}>dyn.B.lo.Xo(n) & current_grd.X{1}<dyn.B.hi.Xo(n),1,numel(current_grd.X{1})))>0
                        dyn.Uo{i,n}(current_grd.X{1}>dyn.B.lo.Xo(n) & current_grd.X{1}<dyn.B.hi.Xo(n)) = dpm_interpn([dyn.B.lo.Xo(n) dyn.B.hi.Xo(n)],[dyn.B.lo.Uo{i}(n) dyn.B.hi.Uo{i}(n)],  current_grd.X{1}(current_grd.X{1}>dyn.B.lo.Xo(n) & current_grd.X{1}<dyn.B.hi.Xo(n)));
                    end
                    dyn.Uo{i,n}(below) = dyn.B.lo.Uo{i}(n);
                    dyn.Uo{i,n}(above) = dyn.B.hi.Uo{i}(n);
            end
        end
    else
        for i=1:length(grd.Nu)
            eval(['dyn.Uo{i,n}(:' repmat(',:',1,length(grd.Nx)-1) ') = current_grd.U{i}(uo' num2str(i) ');']);
        end
    end
    
    % store cost-to-go if savemap options is set
    if options.SaveMap
        for i=1:length(dyn.Jo)
            eval(['V_map{i,n}(:' repmat(',:',1,length(grd.Nx)-1) ') = dyn.Jo{i};']);
        end
    end

    % Update progres bar for the dynamic programming backward
    if ~isdeployed && strcmp(options.Waitbar,'on')
        waitbar(n/dis.N,h);
    end
    if ~isdeployed && strcmp(options.Verbose,'on') && mod(n-1,floor(dis.N/100))==0 && round(100*n/dis.N)<100
        fprintf('%s%2d %%',ones(1,4)*8,round(100*n/dis.N));
    end
    if n<651
        dubug = 1;
    end
end % finish time loop 

% Clear dyn if in low mem mode
if ~options.SaveMap
    dyn.Jo = [];
else
    dyn.Jo = V_map;
end

% Close progres bar
if ~isdeployed && strcmp(options.Waitbar,'on')
    waitbar(1,h)
    close(h)
end

if ~isdeployed && strcmp(options.Verbose,'on')
    fprintf('%s Done!\n',ones(1,5)*8)
end
