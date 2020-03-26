function P=gg_200128_numint_geom_tdep_2(Kini,b,g,M,N,t,marg)
% M = nmax;
if marg
    M=1;
end

l = 0:(M-1);
k = 0:(N-1);
u_ = exp(-2*pi*1i*l/M)-1;
v_ = exp(-2*pi*1i*k/N)-1;
[U,V] = ndgrid(u_,v_); u = U(:); v = V(:);

fun = @(x) (b*exp(-g*x).*(u+g*v*x))./(1-b*exp(-g*x).*(u+g*v*x));
I___ = exp(Kini*integral(fun,0,t,'ArrayValued',true));

I___ = reshape(I___',[M,N])';
P=real(ifft2(I___))';
return

function I = numint(KINI,b,g,u,v,t)
f = @(x) (b*exp(-g*x).*(u+g*v*x))./(1-b*exp(-g*x).*(u+g*v*x));
INTVAL = integral(f,0,t);
I=exp(KINI*INTVAL);
return
