% Shiming Duan, Aug. 18, 2010
% Evaluate int(e^(-s*(t+h))*g(t-h),0,h) at s = Sk
% Convert preshape function g(t) into G(s) at s = Sk
function G = preshape_conv(g0,s,h)

if sum(isnan(g0))>0
    var_g =findsym(g0);
    eval(['syms ' var_g]);
else
    syms t
    var_g = 't';
end

G = [];
for i = 1:length(g0)
    G(i) = vpa(eval(['feval(@int,expm(-s*(' var_g '+h))*g0(i),0,h)']));
end

G = reshape(G,[],1);
