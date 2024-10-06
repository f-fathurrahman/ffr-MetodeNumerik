function [dxdt]=expfun(t,x)
% PURPOSE: Compute numerically the solution 
%  for the differential equation y'=y 
% USAGE: Call from appAp419driverode.m 
%  or  appAp419_420driverode.m
% REVISION HISTORY: 22-May-2014 H.-G. Matuttis
%

dxdt=-x;
return