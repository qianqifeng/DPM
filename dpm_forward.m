function out       	= dpm_forward(dyn,model,par,grd,dis,options)
% DPM_FORWARD   Simulates the model using the optimal input
%
%    Examples:
%
%       Simulates the model using the optimal input
%          out = dpm_forward(dyn,model,par,grd,dis,options)
%
%   -THIS IS THE SECOND CORE OF THE DYNAMIC PROGRAMMING-
%
%   Author(s): Olle L. Sundström, 29.11.2006
%   Copyright 2006- Olle L. Sundström


if ~isfield(options,'InputType') || length(options.InputType)~= length(grd.Nu)
    options.InputType = repmat('c',1,length(grd.Nu));
end


code_grid_states = '';
for i=length(grd.Nx):-1:1
    code_grid_states = [code_grid_states 'current_grd.X{' num2str(i) '},'];
end
code_input_states = '';
for i=length(grd.Nx):-1:1
    code_input_states = [code_input_states ',inp.X{' num2str(i) '}'];
end
code_states_nearest = 'ixm1';
for i=2:length(grd.Nx)
    code_states_nearest = [code_states_nearest ',ixm' num2str(i)];
end

inp  = dpm_get_empty_inp(grd,dis,'zero');
outn = dpm_get_empty_out(model,inp,par,grd,'zero');

% initialize state
for i=1:length(grd.Nx)
    inp.X{i}  = grd.X0{i};
    outn.X{i} = grd.X0{i};
end


% Display progres bar of Dynamic Programming forward iteration
if ~isdeployed && strcmp(options.Waitbar,'on')
    h = waitbar(0,'DP running forwards. Please wait...');
    set(h,'name','DPM:Waitbar');
end

if ~isdeployed && strcmp(options.Verbose,'on')
    fprintf('DP running forwards:     0 %%');
end

% backward compability of dis.N0
if ~isfield(dis,'T0')
    dis.N0 = 1;
end

if ~options.UseUmap
    % GENERATE CODE FOR GENERATING GRID
    if length(grd.Nu)>1
        code_generate_grid = '[';
        for i=1:length(grd.Nu)
            code_generate_grid = [code_generate_grid 'inpt.U{' num2str(i) '} '];
        end
        code_generate_grid = [code_generate_grid '] = ndgrid('];
        for i=1:length(grd.Nu)
            code_generate_grid = [code_generate_grid 'current_grd.U{' num2str(i) '},'];
        end
        code_generate_grid = code_generate_grid(1:end-1);
        code_generate_grid = [code_generate_grid ');'];
    else
        code_generate_grid = 'inpt.U{1} = current_grd.U{1};';
    end
    
    % GENERATE CODE FOR COST-TO-GO INTERPOLATION
    Jsze = size(dyn.Jo);
    for k=1:Jsze(1)
        code_cost_to_go_interp{k} = ['cost_to_go{' num2str(k) '} = dpm_interpn('];
        for i=length(grd.Nx):-1:1
            code_cost_to_go_interp{k} = [code_cost_to_go_interp{k} 'next_grd.X{' num2str(i) '},'];
        end
        code_cost_to_go_interp{k} = [code_cost_to_go_interp{k} 'dyn.Jo{' num2str(k) ',n+1},'];
        for i=length(grd.Nx):-1:1
            code_cost_to_go_interp{k} = [code_cost_to_go_interp{k} 'X{' num2str(i) '},'];
        end
        code_cost_to_go_interp{k} = code_cost_to_go_interp{k}(1:end-1);
        code_cost_to_go_interp{k} = [code_cost_to_go_interp{k} ');'];
    end
end

iswarned = 0;
% Forward sim
for n=dis.N0:dis.N
    for i=1:length(grd.Nx)
        if options.UseLine && length(grd.Nx)==1 && i==1 && ~options.FixedGrid
            current_grd.X{i} = linspace(dyn.B.lo.Xo(n),dyn.B.hi.Xo(n),grd.Nx{i}(n))';
        elseif options.UseLine && ~options.FixedGrid
            current_grd.X{i} = linspace(min(dyn.B.lo.Xo(:,n)),max(dyn.B.hi.Xo(:,n)),grd.Nx{i}(n))';
        else
            current_grd.X{i} = linspace(grd.Xn{i}.lo(n),grd.Xn{i}.hi(n),grd.Nx{i}(n))';
        end
    end
    for i=1:length(grd.Nu)
        current_grd.U{i} = linspace(grd.Un{i}.lo(n),grd.Un{i}.hi(n),grd.Nu{i}(n))';
    end
    
    %make inp struct
    for w = 1:length(dis.W)
        inp.W{w} = dis.W{w}(n);
    end
    inp.Ts   = dis.Ts;
    
    % if Umap is not used, try all possible u-candidates:
    if ~options.UseUmap
        for i=1:length(grd.Nx)
            next_grd.X{i} = linspace(grd.Xn{i}.lo(n+1),grd.Xn{i}.hi(n+1),grd.Nx{i}(n+1))';
        end
        
        % set up grid of all possible input candidates
        % make inp struct
        for w = 1:length(dis.W)
            inpt.W{w} = dis.W{w}(n);
        end
        inpt.Ts   = dis.Ts;
        eval(code_generate_grid)
        for i=1:length(grd.Nx)
            inpt.X{i} = inp.X{i}*ones(size(inpt.U{1}));
        end
        
        % input to system
        [X C I] = feval(model,inpt,par);
        % take care of bounds
        for i=1:length(grd.Nx)
            I = bitor(I,X{i}>grd.Xn{i}.hi(n+1));
            I = bitor(I,X{i}<grd.Xn{i}.lo(n+1));
        end
        % if line was caluclated take care of bounds again
        if options.UseLine
            I = bitor(I,X{i}>dyn.B.hi.Xo(n+1));
            I = bitor(I,X{i}<dyn.B.lo.Xo(n+1));
        end
        
        %arc cost
        J  = (I==0).*C{1} + I.* options.MyInf;
        
        if options.UseLevelSet
            % get information about feasiblity:
            eval(code_cost_to_go_interp{end});
            J(cost_to_go{end}>0) = options.MyInf;
        end
        
        % minimize total cost
        % Calculate cost for entire grid
        if length(dyn.Jo{1})==1
            cost_to_go{1} = dyn.Jo{1};
        else
            %interpolate from cost to go map
            eval(code_cost_to_go_interp{1});
        end
        Jt = J + cost_to_go{1};
        
        %if no valid input signal is found:
        if options.UseLevelSet && ~any(reshape(cost_to_go{end},numel(cost_to_go{end}),1)<=0 & reshape(I,numel(I),1)==0)
            Jt = cost_to_go{2};
            Jt(I~=0) = options.MyInf;
        end
        
        if options.Minimize
            Jt(Jt>options.MyInf) = options.MyInf;
        else
            Jt(Jt<options.MyInf) = options.MyInf;
        end
        
        if length(grd.Nu)>1
            Jt = reshape(Jt,[1,numel(Jt)]);
        end
        
        % minimize the cost-to-go
        if options.Minimize
            [Q ui] = min(Jt);
        else
            [Q ui] = max(Jt);
        end
        
        % use input that minimizes total cost
        for i=1:length(inpt.U)
            if length(grd.Nu)>1
                inpt.U{i} = reshape(inpt.U{i},[1,numel(inpt.U{i})]);
            end
            inp.U{i} = inpt.U{i}(ui(1));
        end
        
    else
        if options.UseLine
            for i=1:length(grd.Nu)
                xi = cell(1,length(grd.Nx));
                for j=length(grd.Nx):-1:1
                    if options.InputType(i) == 'd' % Discrete
                        if j==1
                            xistr = dpm_code('xi{#},',2:length(current_grd.X));
                            eval(['x1vec = [dyn.B.lo.Xo(1,' xistr 'n); current_grd.X{j}(current_grd.X{j}>dyn.B.lo.Xo(1,' xistr 'n) & current_grd.X{j}<dyn.B.hi.Xo(1,' xistr 'n)); dyn.B.hi.Xo(1,' xistr 'n)];']);
                            [temp xi{j}] = min(abs(x1vec-inp.X{j}));
                            Xin{j} = x1vec(xi{j});
                        else
                            xi{j} = round((inp.X{j}-current_grd.X{j}(1))/(current_grd.X{j}(2)-current_grd.X{j}(1))) + 1;
                            xi{j} = max(xi{j},1);
                            xi{j} = min(xi{j},length(current_grd.X{j}));
                            Xin{j} = current_grd.X{j}(xi{j});
                        end
                    else % Continuous
                        Xin{j} = inp.X{j};
                    end
                end
                
                if length(grd.Nx) > 1
                    inp.U{i}   = dpm_interpf2sbh(current_grd.X{1},current_grd.X{2}, dyn.Uo{i,n}(:,:), Xin{1},Xin{2}, [dyn.B.lo.Xo(:,:,n); dyn.B.hi.Xo(:,:,n)]);
                else
                    inp.U{i}   = dpm_interpf1sbh(current_grd.X{1}, dyn.Uo{i,n}, Xin{1}, [dyn.B.lo.Xo(n) dyn.B.hi.Xo(n)]);
                end       %  current_grd.X{i} = linspace(dyn.B.lo.Xo(n),dyn.B.hi.Xo(n),grd.Nx{i}(n))';
            end
        elseif ~isempty(grd.Nu)
            for i=1:length(grd.Nu)
                if options.InputType(i) == 'c' % Continuous
                    eval(['inp.U{i}   = dpm_interpn(' code_grid_states 'dyn.Uo{i,n}(:' repmat(',:',1,length(grd.Nx)-1) ')' code_input_states ');']);
                else % Discrete
                    for j=1:length(grd.Nx)
                        eval(['ixm' num2str(j) ' = round((inp.X{j}-current_grd.X{j}(1))/(current_grd.X{j}(2)-current_grd.X{j}(1)) + 1);']);
                    end
                    eval(['inp.U{i}    = dyn.Uo{i,n}(' code_states_nearest ');']);
                end
            end
        end
    end
    
    %call model with optimal input:
    [X C I outn] = feval(model,inp,par);
    outn.X = inp.X;
    outn.C = C;
    outn.I = I;
    if outn.I~=0 && ~iswarned
        warning('DPM:Forward','Infeasible Solution!')
        iswarned = 1;
    end
    inp.X = X;
    

    if n > dis.N0
        out = dpm_mergestruct(out,outn);
    else
        out = outn;
    end

    % Update progres bar for the Dynamic Programming
    if ~isdeployed && strcmp(options.Waitbar,'on')
        waitbar((n-dis.N0)/(dis.N-dis.N0),h);
    end   
    
    if ~isdeployed && strcmp(options.Verbose,'on') && mod(n-1,floor(dis.N/100))==0 && round(100*n/dis.N)<100
        fprintf('%s%2d %%',ones(1,4)*8,round(100*n/dis.N));
    end

end


for i=1:length(outn.X)
    out.X{i} = [out.X{i} inp.X{i}];
end


% Close progres bar
if ~isdeployed && strcmp(options.Waitbar,'on')
    waitbar(1,h)
    close(h)
end

if ~isdeployed && strcmp(options.Verbose,'on')
    fprintf('%s Done!\n',ones(1,5)*8);
end
