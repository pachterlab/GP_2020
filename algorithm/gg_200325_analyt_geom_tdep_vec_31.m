%compute the solution to bursty transcription/splicing/degradation system
%using Taylor and Laurent approximations of provided order, on an MxN grid,
%at time T. 
function P=gg_200325_analyt_geom_tdep_vec_31(Kini,b,g,M,N,T,...
    marg,N_approx_taylor,N_approx_laurent)
%check if a marginal distribution is desired
if strcmp(marg,'mature')
    M=1;
elseif strcmp(marg,'nascent')
    N=1;
end

%define grid
l = 0:(M-1);
k = 0:(N-1);
th1 = -2*pi*l/M;
th2 = -2*pi*k/N;
u_ = exp(-1i*th1)-1;
v_ = exp(-1i*th2)-1;
[TH1,TH2] = ndgrid(th1,th2); TH1 = TH1(:); TH2 = TH2(:);
[U,V] = ndgrid(u_,v_); U = U(:); V = V(:);
[I,J] = ndgrid(1:M,1:N); I = I(:); J = J(:);

%define threshold to switch between approximations
alpha=((1+sqrt(3))/2)/b;

%identify where function |U| hits this threshold
[roots,~,~]=root_ID(g,TH1,TH2,alpha,I,J);
num_intervals = sum(roots<T,2)+1;
interval_ID_taylor = mod(sum(isfinite(roots),2),2)==0;
roots = [zeros(M*N,1) roots];

%initialize variables to store integral values
INTVAL = zeros(M*N,1);
ICOMP = zeros(M*N,1);

%iterate over successive intervals and compute approximations
for in_ = 1:(max(num_intervals)-1)
    %identify and compute Taylor approximation
    active = in_ <= (num_intervals-1);
    active_taylor = interval_ID_taylor & active;
    ICOMP(active_taylor) = ...
        analytint_taylor(...
        b,g,U(active_taylor),V(active_taylor), ...
        I(active_taylor),J(active_taylor),...
        roots(active_taylor,in_), roots(active_taylor,in_+1),  ...
        N_approx_taylor);
    
    %identify and compute Laurent approximation
    active_laurent = (~interval_ID_taylor) & active;
    ICOMP(active_laurent) = ...
        analytint_laurent(...
        b,g,U(active_laurent),V(active_laurent), ...
        I(active_laurent),J(active_laurent),...
        roots(active_laurent,in_), roots(active_laurent,in_+1),  ...
        N_approx_laurent);
    INTVAL(active) = INTVAL(active)+ICOMP(active);
    interval_ID_taylor = ~interval_ID_taylor;
end


%identify and compute Taylor approximation after last root
final_taylor = mod(sum(isfinite(roots),2)-sum(roots<T,2),2)==0;
final_root = roots;
final_root(final_root>T)=NaN;
final_root = max(final_root,[],2);
final_laurent = ~final_taylor;
ICOMP(final_taylor) = ...
    analytint_taylor(...
    b,g,U(final_taylor),V(final_taylor), ...
    I(final_taylor),J(final_taylor),...
    final_root(final_taylor), T*ones(sum(final_taylor),1),  ...
    N_approx_taylor);

%identify and compute Laurent approximation after last root
ICOMP(final_laurent) = ...
    analytint_laurent(...
    b,g,U(final_laurent),V(final_laurent), ...
    I(final_laurent),J(final_laurent),...
    final_root(final_laurent), T*ones(sum(final_laurent),1),  ...
    N_approx_laurent);
INTVAL = INTVAL + ICOMP;

%compute characteristic function
I___ = exp(Kini*INTVAL);
I___ = reshape(I___,[M,N])';
%convert back to discrete probability domain using IDFT
P=real(ifft2(I___))';
return

%compute all roots and extrema of |U|-alpha
function [s,r,s0]=root_ID(g,th1,th2,alpha,i,j)
n = length(th1);

a2=alpha.^2;
nr_tol=1e-4;
bisec_tol = nr_tol;

u2=2*(1-cos(th1));
v2=2*(1-cos(th2));
uvvu=2*(1+cos(th2-th1)-cos(th1)-cos(th2));

s=NaN(n,3);
r=NaN(n,3);

exact=true(n,1);

%%%%%
% filter for i=j=1
filt_0 = (i==1 & j==1);

%%%%%
% filter for u=0
filt_1 = (i==1 & j>1); %filt_1_ind = find(filt_1); 
filt_1_rt = filt_1 & (v2/exp(2) > a2);
s(filt_1_rt,1) = -1/g * lambertw(0,-alpha./(sqrt(v2(filt_1_rt))));
s(filt_1_rt,2) = -1/g * lambertw(-1,-alpha./(sqrt(v2(filt_1_rt))));

%%%%%
% filter for v=0
filt_2 = (j==1 & i>1);
filt_2_rt = filt_2 & (u2 > a2);
s(filt_2_rt,1) = -log(alpha./(sqrt(u2(filt_2_rt))))/g;


%%%%%
% define subset that requires numerical root finding

filt_inexact = ~(filt_0 | filt_1 | filt_2);
exact(filt_inexact) = false;
inexact_ind = find(filt_inexact);

%%%%%
% calculate polynomial parameters for the derivative function
B = g * (v2(filt_inexact) - uvvu(filt_inexact));
A = -g^2 * v2(filt_inexact);
C = uvvu(filt_inexact)/2 - u2(filt_inexact);
discrim = B.^2 - 4.*A.*C;


% 

r1 = (-B+sqrt(discrim))./(2*A); r1(discrim<0 | r1<0) = NaN;
r2 = (-B-sqrt(discrim))./(2*A); r2(discrim<0 | r2<0) = NaN;
r_extrema = [zeros(size(r1)), r1, r2]; r_extrema=sort(r_extrema,2);

%%%%%%%%%%%%
%%%%%%%%%%%
val_func = @(s,v2,u2,uvvu) exp(-2*g*s).*(g^2*v2.*s.^2 + g*uvvu.*s + u2);
val_0 = u2(filt_inexact);
val_at_extrema = [val_0, val_func(...
    r_extrema(:,2:3),v2(filt_inexact),u2(filt_inexact),uvvu(filt_inexact))];

val_at_extrema = val_at_extrema-a2;
s_guess = NaN(size(r_extrema));
for i = 1:2
    filt = all(~isnan(val_at_extrema(:,i:i+1)),2) & ...
        (sign(val_at_extrema(:,i)) ~= sign(val_at_extrema(:,i+1)));
    s_guess(filt,i) = r_extrema(filt,i) + ...
        (r_extrema(filt,i+1) - r_extrema(filt,i)) .* ...
        val_at_extrema(filt,i)./(val_at_extrema(filt,i)-val_at_extrema(filt,i+1));
end
last_extremum_location = sum(~isnan(r_extrema),2);

for i = 1:3
    filt_loc = last_extremum_location==i;
    root_filt_ind = find(filt_loc & (val_at_extrema(:,i) > 0));
    s_guess(root_filt_ind,3) = ...
        r_extrema(root_filt_ind,i) ...
        - log(a2./(a2+val_at_extrema(root_filt_ind,i)))/(2*g);
end

r(filt_inexact,:) = r_extrema;
s(filt_inexact,:) = s_guess;

%%%%%
% save initial guess
s0 = s(filt_inexact,:);
if sum(~exact)>0
    P1 = [g^2*v2(filt_inexact), g*uvvu(filt_inexact), u2(filt_inexact)];
    P2 = [A,B,C];
    [s_nr,fail]=nr(s0,P1, P2, g, a2, nr_tol);
%     fprintf('%i out of %i NR searches failed.\n',sum(fail),sum(~exact))
    s(inexact_ind(~fail),:)=s_nr(~fail,:);
    
    
    errf = @(x,P1) (P1(:,1).*x.^2 + P1(:,2).*x + P1(:,3)).*exp(-2*g*x) - a2;

    if sum(fail)>0
        filt_1 = fail & ~isnan(s0(:,1)); n = sum(filt_1);
        if n>0
            if any(isnan(r(inexact_ind(filt_1),2))) 
                error('invalid bisection 1!')
            end
            s(inexact_ind(filt_1),1) = bisec(r(inexact_ind(filt_1),1), r(inexact_ind(filt_1),2),...
                P1(filt_1,:),g,a2,bisec_tol);
            if max(errf(s(inexact_ind(filt_1),1),P1(filt_1,:))) > bisec_tol || ...
                any(~isfinite(errf(s(inexact_ind(filt_1),1),P1(filt_1,:))))
                error('bisec 1 false positive!');
            end
        end
        
        filt_2 = fail & ~isnan(s0(:,2));  n=sum(filt_2);
        if n>0
            if any(isnan(r(inexact_ind(filt_2),2))) || any(isnan(r(inexact_ind(filt_2),3)))
                error('invalid bisection 2!')
            end
            s(inexact_ind(filt_2),2) = bisec( ...
                r(inexact_ind(filt_2),2), ...
                r(inexact_ind(filt_2),3),...
                P1(filt_2,:) ,g,a2,bisec_tol);
            if max(errf(s(inexact_ind(filt_2),2),P1(filt_2,:))) > bisec_tol|| ...
                any(~isfinite(errf(s(inexact_ind(filt_2),2),P1(filt_2,:))))
                error('bisec 2 false positive!');
            end
        end
        
        filt_3 = fail & ~isnan(s0(:,3)); n=sum(filt_3);
        if n>0
            maxr = max(r(inexact_ind(filt_3),:),[],2);
            s(inexact_ind(filt_3),3) = bisec( ...
                maxr, ...
                maxr+10/g,...
                P1(filt_3,:),g,a2,bisec_tol);
            if any(isnan(maxr)) 
                error('invalid bisection 3!')
            end
            if max(errf(s(inexact_ind(filt_3),3),P1(filt_3,:))) > bisec_tol || ...
                any(~isfinite(errf(s(inexact_ind(filt_3),3),P1(filt_3,:))))
                error('bisec 3 false positive!');
            end
        end
    end
    if any(any(errf(s(inexact_ind,:),P1)> bisec_tol))  
        error('false positive!!');
    end

end
s(s==0)=NaN;
s=sort(s,2);
return

%use Newton-Raphson method to find roots of |U|-alpha
function [x,fail]=nr(x,P1,P2,g,a2,tol)
n=size(x,1);
x0=x;
err=inf*ones(n,3);
i=1;
imax=20;
relax=1;
fail=false(n,1);

while max(abs(err))>tol
    err = (P1(:,1).*x.^2 + P1(:,2).*x + P1(:,3)).*exp(-2*g*x) - a2;
    x = x - relax * err ./ (2*g*(P2(:,1).*x.^2 + P2(:,2).*x + P2(:,3)).*exp(-2*g*x));
    i=i+1;
    if i>imax 
        fail(any(abs(err)>tol,2) | ...
            (sum(isfinite(x0),2) ~= sum(isfinite(x),2)) | ...
            sum(x<0,2) > 0) = true;
        break
    end
end

fail(any(abs(err)>tol,2) | ...
    (sum(isfinite(x0),2) ~= sum(isfinite(x),2)) | ...
    sum(x<0,2) > 0) = true; 
return

%use bisection method to find roots of |U|-alpha
function c=bisec(a,b,P1,g,a2,tol)
EC=inf;
errf = @(x) (P1(:,1).*x.^2 + P1(:,2).*x + P1(:,3)).*exp(-2*g*x) - a2;
if any(sign(errf(a))==sign(errf(b)))
    error('Wrong bracketing!')
end
while any(abs(EC)>tol)
    c = (a+b)/2;
    EC = errf(c);
    
    filt = sign(EC)==sign(errf(a));
    a(filt) = c(filt);
    b(~filt) = c(~filt);
end
return

%compute weights (omega) for Taylor approximation power series
function [lnw,w]=weight_calc_taylor(b,N)
lnw=NaN(1,N);
% w=lnw;
w=NaN;
k = 1:N;
% lnw=k.*log(b./k) + ...
%         log(1-exp(-(N+2)*log(2) + gammaln(N+2) - gammaln(k+1) - gammaln(N+2-k));
for i = k
    lnw(i) =  log(mhygfx(1,N+2,-i+N+2,1/2));
end
lnw =k.*log(b./k) + ...
        log(1-exp(-(N+2)*log(2) + gammaln(N+2) - gammaln(k+1) - gammaln(N+2-k) ...
        + lnw));
return

%compute integral of Taylor approximation
function INTVAL = analytint_taylor(b,g,u,v,I,J,t_1,t_2,N_approx)
n=length(u);
XX = NaN(n,N_approx);

lnw = weight_calc_taylor(b,N_approx);
k=1:N_approx;

%%%%%
% filter for i = j = 1
filt_0 = I==1 & J==1;

%%%%%
% filter for j = 1 -> v = 0
filt_1 = J==1 & I>1; n_filt_1 = sum(filt_1);
if n_filt_1>0
    XX(filt_1,:) = k.*log(u(filt_1)) + k.*log(k) - log(k*g) ...
                 + log(exp(-g*k.*t_2(filt_1))-exp(-g*k.*t_1(filt_1)));
end

%%%%%
% filter for i = 1 -> u = 0
filt_2 = I==1 & J>1; n_filt_2=sum(filt_2);
if n_filt_2>0
    XX(filt_2,:) = k.*log(v(filt_2)) - log(k*g) + gammaln(k+1) + ...
            log(gammainc(real(g*k.*t_2(filt_2)),repmat(k+1,[n_filt_2,1]),'upper')-...
            gammainc(real(g*k.*t_1(filt_2)),repmat(k+1,[n_filt_2,1]),'upper'));
end

%%%%%
% filter for i,j > 1
filt_3 = ~(filt_0 | filt_1 | filt_2); n_filt_3 = sum(filt_3);
if n_filt_3>0
    MM = exp_int(N_approx,u(filt_3),v(filt_3),g,t_2(filt_3)) ...
        - exp_int(N_approx,u(filt_3),v(filt_3),g,t_1(filt_3));
    XX(filt_3,:) = gammaln(k+1) +k.*log(v(filt_3))-log(g*k)+ log(MM);
end

INTVAL = sum(-exp(XX+lnw),2);
INTVAL(filt_0) = 0;
return

%compute integral of Laurent approximation
function INTVAL = analytint_laurent(b,g,u,v,I,J,t_1,t_2,N_approx)

n=length(u);
XX = NaN(n,N_approx);

k=1:N_approx;
lnw = -k*log(b);

%%%%%
% filter for i = j = 1
filt_0 = I==1 & J==1;

%%%%%
% filter for j = 1 -> v = 0
filt_1 = J==1 & I>1;
XX(filt_1,:) = -k.*log(u(filt_1)) - log(k*g) ...
               + log(exp(g*k.*t_2(filt_1))-exp(g*k.*t_1(filt_1)));

%%%%%
% filter for i,j > 1
filt_2 = ~(filt_0 | filt_1 ); n=sum(filt_2);
if n>0
    MM = exp_int_gen(N_approx,u(filt_2),v(filt_2),g,t_2(filt_2)) ...
        - exp_int_gen(N_approx,u(filt_2),v(filt_2),g,t_1(filt_2));
    XX(filt_2,:) = -gammaln(k) - log(g*v(filt_2)) + log(MM);
end

INTVAL = sum(-exp(XX+lnw),2) - (t_2-t_1);
return

%compute Taylor exponential integral e_n(z)
function MM=exp_int(N,u,v,g,t)
n=length(u);
MM=zeros(n,N);
if t==inf
    return
end
k=1:N;

RR = log(u./v+g*t);
QQ = g*t.*k; 
for j = 0:N
    k_ind = max(1,j):N;
    MM(:,k_ind) = MM(:,k_ind) + ...
        exp(j*RR + j*log(k(k_ind))-gammaln(j+1)-QQ(:,k_ind));
end
return

%compute exponential integrak E_N(z)
function MM=exp_int_gen(N,u,v,g,t)
n=length(u);
k=1:N;

z=(k.*u./v+k*g.*t);
MM=zeros(n,N);

RR = log(u+g*v.*t);
QQ = g*t.*k;

for j = 0:(N-2)
    k_ind = (j+2):N;
    MM(:,k_ind) = MM(:,k_ind) - exp(gammaln(k(k_ind)-j-1)+QQ(:,k_ind) + j.*log(k(k_ind)./v) + (j+1-k(k_ind)) .* RR);
end
MM = MM+exp(-k.*u./v + (k-1).*log(k./v)) .* (-me1z_gg200131_comb_4(-z));

return
