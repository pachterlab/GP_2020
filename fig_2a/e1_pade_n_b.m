% 6th order Padé approximation (Luke 1975, p.110-112 , table 4.3)
function ei_ = e1_pade_n_b(z,n)
el=0.5772156649015328;
switch n
    case 0
        p1 = [1];
        p2 = [1];
    case 1
        p1 = [1,1];
        p2 = [1,2];
    case 2
        p1 = [1,5,2];
        p2 = [1,6,6];
    case 3
        p1 = [1,11,26,6];
        p2 = [1,12,36,24];
    case 4
        p1 = [1,19,102,154,24];
        p2 = [1,20,120,240,120];
    case 5
        p1 = [1,29,272,954,1044,120];
        p2 = [1,30,300,1200,1800,720];
    case 6
        p1 = [1,41,590,3648,9432,8028,720];
        p2 = [1,42,630,4200,12600,15120,5040];
end
if length(p1)~=length(p2)
    error('incorrect input');
end
ei_=polyval(p1,z)./polyval(p2,z)./z .*exp(-z) ;
return