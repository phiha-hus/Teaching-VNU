function [u]=beonestep(t,u,y,h,f,fn,varargin)
% BEONESTEP  one step of the backward Euler method
u = u + h*f(t,y,varargin{:});
return
