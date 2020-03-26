function P=gg_200130_analyt_geom_tdep_vec_26(Kini,b,g,M,N,T,marg,N_approx_taylor,N_approx_laurent)
if marg
    M=1;
end

l = 0:(M-1);
k = 0:(N-1);
th1 = -2*pi*l/M;
th2 = -2*pi*k/N;
u_ = exp(-1i*th1)-1;
v_ = exp(-1i*th2)-1;

[TH1,TH2] = ndgrid(th1,th2); TH1 = TH1(:); TH2 = TH2(:);
[U,V] = ndgrid(u_,v_); U = U(:); V = V(:);
[I,J] = ndgrid(1:M,1:N); I = I(:); J = J(:);


alpha=((1+sqrt(3))/2)/b;
roots=root_ID(g,TH1,TH2,alpha,I,J);

num_intervals = sum(roots<T,2)+1;
interval_ID_taylor = mod(sum(isfinite(roots),2),2)==0;
roots = [zeros(M*N,1) roots];

INTVAL = zeros(M*N,1);
ICOMP = zeros(M*N,1);

for in_ = 1:(max(num_intervals)-1)
    active = in_ <= (num_intervals-1);
    active_taylor = interval_ID_taylor & active;
    ICOMP(active_taylor) = ...
        analytint_taylor(...
        b,g,U(active_taylor),V(active_taylor), ...
        I(active_taylor),J(active_taylor),...
        roots(active_taylor,in_), roots(active_taylor,in_+1),  ...
        N_approx_taylor);
    
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

ICOMP(final_laurent) = ...
    analytint_laurent(...
    b,g,U(final_laurent),V(final_laurent), ...
    I(final_laurent),J(final_laurent),...
    final_root(final_laurent), T*ones(sum(final_laurent),1),  ...
    N_approx_laurent);
INTVAL = INTVAL + ICOMP;

I___=exp(Kini*INTVAL);

I___ = reshape(I___,[M,N])';
P=real(ifft2(I___))';
return

function s=root_ID(g,th1,th2,alpha,i,j)
n = length(th1);

a2=alpha.^2;
nr_tol=1e-4;
bisec_tol = nr_tol;

u2=2*(1-cos(th1));
v2=2*(1-cos(th2));
uvvu=2*(1+cos(th2-th1)-cos(th1)-cos(th2));

s=NaN(n,3);
r=NaN(n,2);

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
% calculate polynomial parameters for subset
B = g * (v2(filt_inexact) - uvvu(filt_inexact));
A = -g^2 * v2(filt_inexact);
C = uvvu(filt_inexact)/2 - u2(filt_inexact);
discrim = B.^2 - 4.*A.*C;
val_0 = u2(filt_inexact);

%%%%%
% filter cases with 1 positive root
filt_r1 = C>0; filt_r1_ind =find(filt_r1);
r1 = (-B(filt_r1)-sqrt(discrim(filt_r1)))./(2*A(filt_r1));
val_r1 = exp(-2*g*r1).*(g^2*v2(inexact_ind(filt_r1)).*r1.^2 ...
    + g*uvvu(inexact_ind(filt_r1)).*r1+u2(inexact_ind(filt_r1)));
%intersection right of the root
filt_r1_s3 = val_r1>a2; %filt_r1_s3_ind = find(filt_r1_s3);
s(inexact_ind(filt_r1_ind(filt_r1_s3)),3) = r1(filt_r1_s3)-log(a2./val_r1(filt_r1_s3))/(2*g);
r(inexact_ind(filt_r1_ind(filt_r1_s3)),1) = r1(filt_r1_s3);
%intersection left of the root
filt_r1_s1 = filt_r1_s3 & (val_0(filt_r1) < a2); 
s(inexact_ind(filt_r1_ind(filt_r1_s1)),1) = ...
    r1(filt_r1_s1) .* (a2-val_0(filt_r1_ind(filt_r1_s1))) ./ ...
                      (val_r1(filt_r1_s1) - val_0(filt_r1_ind(filt_r1_s1)));


%%%%%
% filter cases with no positive roots
filt_r0 = ~filt_r1 & ((discrim<0) | (B./A>0)); filt_r0_ind = find(filt_r0);
filt_r0_s3 = val_0(filt_r0)>a2; 
s(inexact_ind(filt_r0_ind(filt_r0_s3)),3) = ...
    -log(a2./val_0(filt_r0_ind(filt_r0_s3)))/(2*g);
r(inexact_ind(filt_r0_ind(filt_r0_s3)),2) = 0;


%%%%%
% filter cases with two positive roots
filt_r2 = ~(filt_r1 | filt_r0); filt_r2_ind = find(filt_r2);
r1 = (-B(filt_r2) + sqrt(discrim(filt_r2))) ./ (2*A(filt_r2));
r2 = (-B(filt_r2) - sqrt(discrim(filt_r2))) ./ (2*A(filt_r2));
% if any(isnan(r1)) || any(isnan(r2))
%     error('')
% end
val_r1 = exp(-2*g*r1) .*...
    (g^2*v2(inexact_ind(filt_r2)) .* r1.^2 ...
    + g*uvvu(inexact_ind(filt_r2)) .* r1 ...
    + u2(inexact_ind(filt_r2)));
val_r2 = exp(-2*g*r2) .*...
    (g^2*v2(inexact_ind(filt_r2)) .* r2.^2 ...
    + g*uvvu(inexact_ind(filt_r2)) .* r2 ...
    + u2(inexact_ind(filt_r2)));
%intersection left of left root
filt_r2_s1 = (val_0(filt_r2) > a2) & (val_r1<a2);
s(inexact_ind(filt_r2_ind(filt_r2_s1)),1) = ...
    r1(filt_r2_s1) .* (a2 - val_0(filt_r2_ind(filt_r2_s1))) ./ ...
              (val_r1(filt_r2_s1) - val_0(filt_r2_ind(filt_r2_s1)));
r(inexact_ind(filt_r2_ind(filt_r2_s1)),1) = r1(filt_r2_s1);
%intersection right of right root
filt_r2_s3 = val_r2 > a2;
s(inexact_ind(filt_r2_ind(filt_r2_s3)),3) = ...
    r2(filt_r2_s3)-log(a2./val_r2(filt_r2_s3))/(2*g);
r(inexact_ind(filt_r2_ind(filt_r2_s3)),2) = r2(filt_r2_s3);
%intersection between roots
filt_r2_s2 = filt_r2_s3 & (val_r1<a2);
s(inexact_ind(filt_r2_ind(filt_r2_s2)),2) = ...
    r1(filt_r2_s2) + (r2(filt_r2_s2)-r1(filt_r2_s2))...
              .* (a2 - val_r1(filt_r2_s2)) ./ ...
              (val_r2(filt_r2_s2) - val_r1(filt_r2_s2));
r(inexact_ind(filt_r2_ind(filt_r2_s2)),1) = r1(filt_r2_s2);

% find(~isnan(s(:,2)))
% disp('');
          
%%%%%
% save initial guess
s0 = s(filt_inexact,:);
if sum(~exact)>0
    P1 = [g^2*v2(filt_inexact), g*uvvu(filt_inexact), u2(filt_inexact)];
    P2 = [A,B,C];
    [s_nr,fail]=nr(s0,P1, P2, g, a2, nr_tol);
%     fprintf('%i NR searches failed.\n',sum(fail))
    s(inexact_ind(~fail),:)=s_nr(~fail,:);
    
    
    errf = @(x,P1) (P1(:,1).*x.^2 + P1(:,2).*x + P1(:,3)).*exp(-2*g*x) - a2;

    if sum(fail)>0
        filt_1 = fail & ~isnan(s0(:,1)); n = sum(filt_1);
        if n>0
            s(inexact_ind(filt_1),1) = bisec(zeros(n,1), r(inexact_ind(filt_1),1),...
                P1(filt_1,:) ,g,a2,bisec_tol);
            if any(isnan(r(inexact_ind(filt_1),1))) 
                error('invalid bisection 1!')
            end
            if max(errf(s(inexact_ind(filt_1),1),P1(filt_1,:))) > bisec_tol || ...
                any(~isfinite(errf(s(inexact_ind(filt_1),1),P1(filt_1,:))))
                error('bisec 1 false positive!');
            end
        end
        
        filt_2 = fail & ~isnan(s0(:,2));  n=sum(filt_2);
        if n>0
            if any(isnan(r(inexact_ind(filt_2),1))) || any(isnan(r(inexact_ind(filt_2),2)))
                error('invalid bisection 2!')
            end
            s(inexact_ind(filt_2),2) = bisec( ...
                r(inexact_ind(filt_2),1), ...
                r(inexact_ind(filt_2),2),...
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
                error('bisec 2 false positive!');
            end
        end
%     s(inexact_ind(filt_1),1)
    end
    if any(any(errf(s(inexact_ind,:),P1)> bisec_tol))  
        error('false positive!!');
    end

end
s(s==0)=NaN;
s=sort(s,2);
return
function c=bisec(a,b,P1,g,a2,tol)
EC=inf;
errf = @(x) (P1(:,1).*x.^2 + P1(:,2).*x + P1(:,3)).*exp(-2*g*x) - a2;
% a'
% b'
% errf(a)
% errf(b)
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
% errf(c)
return
function [x,fail]=nr(x,P1,P2,g,a2,tol)
n=size(x,1);
x0=x;
err=inf*ones(n,3);
i=1;
imax=20;
n_restarts = 0;
relax=1;
max_restarts=3;
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

function [lnw,w]=weight_calc_taylor(b,N)
lnw=NaN(1,N);
% w=lnw;
w=NaN;
k = 1:N;
for i = k
    lnw(i) =  log(mhygfx(1,N+2,-i+N+2,1/2));
end
lnw =k.*log(b./k) + ...
        log(1-exp(-(N+2)*log(2) + gammaln(N+2) - gammaln(k+1) - gammaln(N+2-k) ...
        + lnw));
return
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
            log(gammainc(g*k.*t_2(filt_2),repmat(k+1,[n_filt_2,1]),'upper')-...
            gammainc(g*k.*t_1(filt_2),repmat(k+1,[n_filt_2,1]),'upper'));
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
function MM=exp_int(N,u,v,g,t)
n=length(u);
MM=zeros(n,N);
%for Taylor
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
