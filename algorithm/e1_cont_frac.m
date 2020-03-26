function ct = e1_cont_frac(z,ord)
n=length(z);
% el=0.5772156649015328;
ct0 = complex(zeros(n,1));
for k = ord:-1:1
    ct0 = k./(1+k./(z+ct0));
end
ct = 1./(z+ct0);
ct=exp(-z) .* ct - pi*1i*(imag(z)==0 & real(z)<0);
return