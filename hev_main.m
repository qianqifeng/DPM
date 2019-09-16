% load driving cycle
load JN1015

% create grid
clear grd
grd.Nx{1}    = 61;  % 点的个数 
grd.Xn{1}.hi = 0.7; % 上下边界
grd.Xn{1}.lo = 0.4;

grd.Nu{1}    = 21; %控制点的个数
grd.Un{1}.hi = 1;   % 控制输入的上下边界
grd.Un{1}.lo = -1;	% Att: Lower bound may vary with engine size.

% set initial state
grd.X0{1} = 0.55;   % 初始值

% final state constraints
grd.XN{1}.hi = 0.56; %终值的上下边界
grd.XN{1}.lo = 0.55;

% define problem
clear prb
prb.W{1} = speed_vector; % (661 elements)
prb.W{2} = acceleration_vector; % (661 elements)
prb.W{3} = gearnumber_vector; % (661 elements)
prb.Ts = 1;
prb.N  = 660*1/prb.Ts + 1;
 
% set options
options = dpm();
options.MyInf = 1000;
options.BoundaryMethod = 'none'; % also possible: 'none' or 'LevelSet';
% options.Waitbar     = 'on';
if strcmp(options.BoundaryMethod,'Line') 
    %these options are only needed if 'Line' is used
    options.Iter = 5;
    options.Tol = 1e-8;
    options.FixedGrid = 0;
end
[res dyn] = dpm(@hev,[],grd,prb,options);
