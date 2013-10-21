function yi = interpol(x,Y,xi)
% interpolation function
% This is a wrapper function for more advanced, like the ones in curve
% fitting toolbox interpolation functions. This is to pretain the (possible
% future) compatibility with Python
%
% Created:
% 11 Oct Antti Ripatti
%
% (c) Antti Ripatti
%
yi = interp1(x,Y,xi,'spline');