function y         	= dpm_interpf2sbh(xx1,xx2,YY,a1,a2,lim)
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
lim2(1) = dpm_interpn(xx2,lim(1,:),a2);
lim2(2) = dpm_interpn(xx2,lim(2,:),a2);
yy      = dpm_interpn(xx2,xx1,YY,a2.*ones(size(xx1)),xx1);
y       = dpm_interpf1sbh(xx1,yy,a1,lim2);
