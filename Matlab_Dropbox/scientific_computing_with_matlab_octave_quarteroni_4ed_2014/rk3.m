function [tt,u]=rk3(odefun,tspan,y0,Nh,varargin);
tt=linspace(tspan(1),tspan(2),Nh+1);
h=(tspan(2)-tspan(1))/Nh; hh=h*0.5; h2=2*h;
u=y0; h6=h/6;
for t=tt(1:end-1)
  y = u(end,:);
  k1=odefun(t,y,varargin{:});
  t1 = t + hh; y1 = y + hh* k1;
  k2=odefun(t1,y1,varargin{:});
  t1 = t + h; y1 = y + h*(2*k2-k1);
  k3=odefun(t1,y1,varargin{:});
  u = [u; u(end,:) + h6*(k1+4*k2+k3)];
end
tt=tt';
