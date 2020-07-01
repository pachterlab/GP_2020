% 55th order series approximation (Zhang and Jin 1996)
function CE1 = e1_series_2(z,ord)
n=length(z);
el=0.5772156649015328;
CE1=complex(ones(n,1));
% tol = 1e-10;
cr=CE1;
k_ = 1:ord;
% k_ = (1:ord) ./ ((1:ord)+1).^2;
% CR = (-z).^k_ ./(k_+1).^2 ./ gamma(k_+1);
% CE1 = 1+ sum(CR,2);
% act = true(size(z));
% zz=z;r
for kV=1:ord
    cr = -cr .* z * kV/(kV+1)^2;
    CE1 = CE1 + cr;
%     zz(act) = ;
%     if all(~act)
%         break;
%     end
end
CE1 = -el - log(z) + z.*CE1;
return