function y         	= dpm_interpf2mb(xx1,xx2,YY,A1,A2,xlim,ylim,myInf)
%MY_INTERPF1M Computes the 1D interpolation for the given set A
%   using the function YY(xx). USES EXTRAPOLATION
%
%   Y = MY_INTERPF1M(XX,YY,A)
%
%   XX    = Axis vector for the value matrix YY
%   YY    = Value matrix
%   A     = Set of XX-values to be interpolated
%
%   Y     = Set of interpolated values (same structure as A)
%
%   Assumes equally spaced XX
%
%   Author(s): Olle L. Sundström, 29.11.2006
%   Copyright 2006- Olle L. Sundström


XX1   = repmat(xx1',length(xx2),1);
XLIMl = repmat(xlim(1,:)',1,length(xx1));
XLIMu = repmat(xlim(2,:)',1,length(xx1));

% find grid point just inside lower boundary
[r c] = find(XX1>XLIMl);
[ru in]=unique(r,'first');
Iinl = c(in);
% Iinl = findl(XX1>XLIMl,'first');

% find grid point just inside upper boundary
[r c] = find(XX1<XLIMu);
[ru in]=unique(r,'last');
Iinu = c(in);
% Iinu = findl(XX1<XLIMu,'last');

% Iinl = find(repmat(xx1',7,1)>repmat(xlim(:,1),1,31),1,'first');
% Iinu = find(xx<xlim(2),1,'last');
xliml = xlim(1,:);
xlimu = xlim(2,:);
yliml = ylim(1,:);
ylimu = ylim(2,:);
xx1l  = xx1(Iinl);
xx1u  = xx1(Iinu);

% find interpolation points between lower boundary and closest grid point
Ibel = find(A1 >=  dpm_interpn(xx2,xliml,A2) & A1 <  dpm_interpn(xx2,xx1l,A2));
% find interpolation points between upper boundary and closest grid point
Ibeu = find(A1 <=  dpm_interpn(xx2,xlimu,A2) & A1 >  dpm_interpn(xx2,xx1u,A2));
% find interpolation points outside boundary
Iout = find(A1 <   dpm_interpn(xx2,xliml,A2) | A1 >  dpm_interpn(xx2,xlimu,A2));
% find interpolation points inside boundary
Iin  = find(A1 >=  dpm_interpn(xx2,xx1l,A2)  & A1 <= dpm_interpn(xx2,xx1u,A2));
% find interpolation points between lower and upper boundary
Iinx = find(A1 >=  dpm_interpn(xx2,xliml,A2) & A1 <= dpm_interpn(xx2,xlimu,A2));

% initialize output
y = nan(size(A1));

% interpolate as usual with interior points
if ~isempty(Iin)
    y(Iin) = dpm_interpn(xx2,xx1,YY,A2(Iin),A1(Iin));
end
% set outside points to inf
if ~isempty(Iout)
    y(Iout)= myInf;
end
% if there are grid points between boundaries
% if ~isempty(find(xx1<xlimu & xx1>xliml,1))
if ~isempty(min(Iinu - Iinl)>0)
    % interpolate points between lower boundary and closest feasible grid point
    if ~isempty(Ibel)
        Xl = dpm_interpn(xx2,xliml,A2(Ibel));
        Xu = dpm_interpn(xx2,xx1l,A2(Ibel));
        Yl = dpm_interpn(xx2,yliml,Xl);
        Yu = dpm_interpn(xx2,YY(dpm_sub2ind(size(YY),Iinl,(1:length(xx2))')),Xu);
        y(Ibel) = Yl + (A1(Ibel) - Xl)./(Xu-Xl).*(Yu-Yl);
        %         y(Ibel) = dpm_interpn([xliml xx1l],[yliml yy(Iinl)],A(Ibel));
    end
    % interpolate points between upper boundary and closest feasible grid point
    if ~isempty(Ibeu)
        Xu = dpm_interpn(xx2,xlimu,A2(Ibeu));
        Xl = dpm_interpn(xx2,xx1u,A2(Ibeu));
        Yu = dpm_interpn(xx2,ylimu,Xu);
        Yl = dpm_interpn(xx2,YY(dpm_sub2ind(size(YY),Iinu,(1:length(xx2))')),Xl);
        y(Ibeu) = Yl + (A1(Ibeu) - Xl)./(Xu-Xl).*(Yu-Yl);
        %         y(Ibeu) = dpm_interpn([xx1u xlimu'],[yy(Iinu) ylim(2)],A2(Ibeu),A1(Ibeu));
    end
else
    % if there are no grid points between boundaries
    y(Iinx) = dpm_interpn(xlim,ylim,A(Iinx));
end
