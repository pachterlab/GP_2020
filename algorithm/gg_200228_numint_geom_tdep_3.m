function P=gg_200228_numint_geom_tdep_3(Kini,b,g,M,N,t,marg)
if strcmp(marg,'mature')
    M=1;
elseif strcmp(marg,'nascent')
    N=1;
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

