function ce1=me1z_gg200131_comb_4(z)
%       ====================================================
%       Purpose: Compute complex exponential integral E1(z)
%       Input :  z   --- Argument of E1(z)
%       Output:  CE1 --- E1(z)
%       ====================================================
SH = size(z);
z = reshape(z,[prod(SH),1]);
Rez = real(z);
a0=abs(z);

%%%%%
% initialize soln array
ce1 = NaN(size(z));

%%%%%
% set singular value
filt = a0==0;
ce1(filt) = inf;

act = ~filt;

%%%%%
% filter 4: Pade around "finger"
filt = act & ((Rez/17+0.3824).^2 + (imag(z)/13).^2 > 1);n=sum(filt);
%filt = act & (abs(imag(z))>12 | (Rez>0 & a0>12)); n=sum(filt);
if n>0
    z_filt = z(filt);
    ce1(filt) = e1_pade_n_b(z_filt,6);
end
% size(filt)
% size(act)
act = act & ~filt;


%%%%%
% filter 2: pade

filt = act & (((Rez+10)/15).^2 + (imag(z)/9.5).^2 > 1); n=sum(filt);
if n>0
    z_filt = z(filt);
    ce1(filt) = e1_pade_n_c(z_filt,10);
end

act = act & ~filt;

%%%%%
% filter : pade 
filt = act & (((Rez+0.65)/4.05).^2 + (imag(z)/4).^2 <1 ) ; n=sum(filt);
if n>0 
%     disp(n)
    z_filt = z(filt);
    ce1(filt) = e1_pade_n(z_filt,10);
end

act = act & ~filt;

%%%%%
% filter : pade 
filt = act & (((Rez+0.65)/4.05).^2 + (imag(z)/4).^2 <1 ) ; n=sum(filt);
if n>0 
%     disp(n)
    z_filt = z(filt);
    ce1(filt) = e1_pade_n(z_filt,10);
end

act = act & ~filt;

%%%%%
% filter : cheb ei gt8
filt = act & (Rez<-8 & (abs(imag(z)) < (-Rez-8)*0.5294)) ; n=sum(filt);
if n>0 
%     disp(n)
    z_filt = z(filt);
    ce1(filt) = -ei_cheb_gt8(-z_filt,20)-sign(angle(z_filt))*1i*pi;
end

act = act & ~filt;
%%%%%
% filter : cheb ei lt8
filt = act & (((Rez+4.5)/4.5).^2 + (imag(z)/2.3).^2 <1 ) ; n=sum(filt);
if n>0 
%     disp(n)
    z_filt = z(filt);
    ce1(filt) = -ei_cheb_lt8(-z_filt,20)-sign(angle(z_filt))*1i*pi;
end

act = act & ~filt;

%%%%%
% filter 1: interior series
filt = act ; n=sum(filt);
if n>0
    z_filt = z(filt);
    ce1(filt) = e1_series(z_filt,55);
end

ce1 = reshape(ce1,SH);
return

