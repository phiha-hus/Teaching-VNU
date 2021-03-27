function [t,u]=newmark(odefun,tspan,y0,Nh,param,...
               varargin)
%NEWMARK Solves second order differential equations
%  using the Newmark method
%  [T,Y]=NEWMARK(ODEFUN,TSPAN,Y0,NH,PARAM) with TSPAN =
%  [T0 TF] integrates the system of differential
%  equations y''=f(t,y,y') from time T0 to TF with
%  initial conditions Y0=(y(t0),y'(t0)) using the
%  Newmark method on an equispaced grid of NH steps.
%  PARAM holds parameters zeta and theta
%  Function ODEFUN(T,Y) must return a vector, whose
%  elements hold the evaluation of f(t,y), of the
%  same dimension of Y.
%  Each row in the solution array Y corresponds to a
%  time returned in the  column vector T.
tt=linspace(tspan(1),tspan(2),Nh+1);
y=y0(:); u=y.';
global glob_h glob_t glob_y glob_odefun;
global glob_zeta glob_theta glob_varargin glob_fn;
glob_h=(tspan(2)-tspan(1))/Nh;
glob_y=y; glob_odefun=odefun;
glob_zeta = param(1); glob_theta = param(2);
glob_varargin=varargin;
if ( exist('OCTAVE_VERSION') )
o_ver=OCTAVE_VERSION;
version=str2num([o_ver(1),o_ver(3),o_ver(5)]);
end

if ( ~exist( 'OCTAVE_VERSION' )  | version >= 320 )
 options=optimset;
 options.Display='off';
 options.TolFun=1.e-12;
 options.MaxFunEvals=10000;
end
glob_fn =odefun(tt(1),glob_y,varargin{:});
for glob_t=tt(2:end)
if ( exist( 'OCTAVE_VERSION' ) & version < 320 )
  w = fsolve('newmarkfun', glob_y );
else
  w = fsolve(@(w) newmarkfun(w),glob_y,options);
end
  glob_fn =odefun(glob_t,w,varargin{:});
  u = [u; w.']; glob_y = w;
end
t=tt';
clear glob_h glob_t glob_y glob_odefun;
clear glob_zeta glob_theta glob_varargin glob_fn;
end

function z=newmarkfun(w)
 global glob_h glob_t glob_y glob_odefun;
 global glob_zeta glob_theta glob_varargin glob_fn;
 fn1=glob_odefun(glob_t,w,glob_varargin{:});
 z(1)=w(1) - glob_y(1) -glob_h*glob_y(2)-...
   glob_h^2*(glob_zeta*fn1+(0.5-glob_zeta)*glob_fn);
 z(2)=w(2) - glob_y(2) -...
   glob_h*((1-glob_theta)*glob_fn+glob_theta*fn1);
end
