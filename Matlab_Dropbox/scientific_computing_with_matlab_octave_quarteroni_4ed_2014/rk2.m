function [tt,u]=rk2(odefun,tspan,y0,Nh,varargin)
tt=linspace(tspan(1),tspan(2),Nh+1);
h=(tspan(2)-tspan(1))/Nh;  hh=h*0.5;
u=y0;
for t=tt(1:end-1)
  y = u(end,:);
  k1=odefun(t,y,varargin{:});
  t1 = t + h; y = y + h*k1;
  k2=odefun(t1,y,varargin{:});
  u = [u; u(end,:) + hh*(k1+k2)];
end
tt=tt';
