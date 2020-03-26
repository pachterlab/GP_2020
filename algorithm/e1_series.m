%Jin & Zhang
function CE1 = e1_series(z,ord)
n=length(z);
el=0.5772156649015328;
CE1=complex(ones(n,1));
cr=CE1;
k_ = (1:ord) ./ ((1:ord)+1).^2;
for kV=1:int32(ord)
    cr = -cr .* z * k_(kV) ;
    CE1 = CE1 + cr;
end
CE1 = -el - log(z) + z.*CE1;
return