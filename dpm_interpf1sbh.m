function y         	= dpm_interpf1sbh(xx,yy,a,lim)
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
xlu = find(xx >  lim(1),1,'first');
xll = find(xx <= lim(1),1,'last');
xuu = find(xx >= lim(2),1,'first');
xul = find(xx <  lim(2),1,'last');

% if a is between lower limit and regular grid
if a <= lim(1)
    y = yy(xll);
    % if a is outside upper limit
elseif a >= lim(2)
    y = yy(xuu);
    % if a is inside limits and within regular grid
elseif a < xx(xlu) && a > lim(1)
%     % if close to engine off
%     if yy(xlu) == 1 || yy(xll) == 1
%         [tmp ind] = min(abs([xx(xlu) lim(1)]-a));
%         ytmp = [yy(xlu) yy(xll)];
%         y = ytmp(ind);
%     else
        dy  = yy(xlu)-yy(xll);
        dx  = xx(xlu)-lim(1);
        y     = (a-lim(1))*dy/dx + yy(xll);
%     end
% if a is between upper limit and regular grid
elseif  a < lim(2) && a > xx(xul)
%     % if close to engine off
%     if yy(xuu) == 1 || yy(xul) == 1
%         [tmp ind] = min(abs([lim(2) xx(xul)]-a));
%         ytmp = [yy(xuu) yy(xul)];
%         y = ytmp(ind);
%     else
        dy  = yy(xuu)-yy(xul);
        dx  = lim(2)-xx(xul);
        y     = (a-xx(xul))*dy/dx + yy(xul);
%     end
    % if a is outside lower limit
else
%     % if close to engine off
%     il = find(xx<a,1,'last');
%     iu = find(xx>a,1,'first');
%     if yy(il) == 1 || yy(iu)==1
%         ix     = round(dpm_interpn([xx(il) xx(iu)],[1 2],a));
%         yy = [yy(il) yy(iu)];
%         y = yy(ix);
%     else
        y     = dpm_interpn(xx,yy,a);
%     end
end
