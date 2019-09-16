function y         	= dpm_interpf1mb(xx,yy,A,xlim,ylim,myInf)
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

% find grid point just inside lower boundary
Iinl = find(xx>xlim(1),1,'first');
% find grid point just inside upper boundary
Iinu = find(xx<xlim(2),1,'last');

% find interpolation points between lower boundary and closest grid point
Ibel = find(A >=  xlim(1) & A < xx(Iinl));
% find interpolation points between upper boundary and closest grid point
Ibeu = find(A <=  xlim(2) & A > xx(Iinu));
% find interpolation points outside boundary
Iout = find(A <   xlim(1) | A > xlim(2));
% find interpolation points inside boundary
Iin  = find(A >= xx(Iinl) & A <=xx(Iinu));
% find interpolation points between lower and upper boundary
Iinx = find(A >=  xlim(1) & A <=xlim(2));

% initialize output
y = zeros(size(A));

% interpolate as usual with interior points
if ~isempty(Iin)
    y(Iin) = dpm_interpn(xx',yy',A(Iin));
end
% set outside points to inf
if ~isempty(Iout)
    y(Iout)= myInf;
end
% if there are grid points between boundaries
if ~isempty(find(xx<xlim(2) & xx>xlim(1),1))
    % interpolate points between lower boundary and closest feasible grid point
    if ~isempty(Ibel)
        y(Ibel)= dpm_interpn([xlim(1) xx(Iinl)],[ylim(1) yy(Iinl)],A(Ibel));
    end
    % interpolate points between upper boundary and closest feasible grid point
    if ~isempty(Ibeu)
        y(Ibeu)= dpm_interpn([xx(Iinu) xlim(2)],[yy(Iinu) ylim(2)],A(Ibeu));
    end
else
    % if there are no grid points between boundaries
    y(Iinx)= dpm_interpn(xlim,ylim,A(Iinx));
end
