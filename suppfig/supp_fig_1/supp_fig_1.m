function supp_fig_1
clc;clear;close all;

%%%%%%
% Initialize the complex plane points to sample
%width of the square in the complex plane, centered on zero
D = 70;
%number of points per side of grid
nnn = 1000;
R = linspace(-D/2,D/2,nnn);
[R_,I_] = ndgrid(R,R); 
z = R_(:) + I_(:)*1i;

purp = [239,171,245]/255;


%calculate the MATLAB built-in benchmark results
tic;
tru = expint(z);
t_matlab = toc;

%define functions to evaluate
func = {  ...
    @(z)e1_pade_n_b(z,6),...
    @(z)e1_pade_n_c(z,10),...
    @(z)e1_pade_n(z,10),...
    @(z)-ei_cheb_lt8(-z,20)-sign(angle(z))*1i*pi,...
    @(z)-ei_cheb_gt8(-z,20)-sign(angle(z))*1i*pi,...
    @(z)e1_series_2(z,55),...
    @(z)me1z_gg200131_comb_4(z)};
n=8;

p1=[227,172,82]/255;
p2=[252,222,164]/255;
p3=[90,180,172]/255;
map = NaN(n*2,3);
for i = 1:3
    map(1:n,i) = linspace(p1(i),p2(i),n);
    map((n+1):end,i) = linspace(p2(i),p3(i),n);
end

name = {'pade1','pade2','pade3','cheb1','cheb2','taylor','Combined'};


%define axes and elliptical regions to draw
axBIG = subplot(2,5,[4,5,9,10]);hold on;
cent = [-6.5,-10,-.65,-4.5];
rad_1 = [17,15,4.05,4.5];
rad_2 = [13,9.5,4,2.3];

%begin to evaluate approximations
for i = 7:-1:1
    
    tic
    est=func{i}(z);
    t_ = toc;
    
    if i>3,JJ=i+2;else JJ=i;end
    if i<7
        ax = subplot(2,5,JJ);
        axes(ax)
    else
        axes(axBIG);
    end
    
    %calculate error value, clamping at 1e-18
    val = log10(abs(est-tru)./abs(tru)); 
    val(val<-18)=-18;
    colormap(map);
    
    %draw behavior of approximation error throughout the complex plane
    [~,HH]=contourf(R_,I_,reshape(val,size(R_)),30); hold on; 
    set(HH,'LineColor','none')

    %center colormap at 1e-8
    caxis([-18,2]);
    axis(D/2*[-1 1 -1 1]);
    title(name{i},'FontWeight','Normal');
    
    %print error and timing diagonstics
    fprintf('%i: %s maxerr: %.3f, Tref/T = %.3f.\n',i,name{i},...
        max(log10(abs(est-tru)./abs(tru))),t_matlab/t_);
    set(gca,'xcolor','none','ycolor','none','FontSize',12,'color','none')

    %plot the approximation bounds
    if i < 5
        th_ = linspace(0,2*pi,100);
        x_ = rad_1(i)*cos(th_)+cent(i);
        y_ = rad_2(i)*sin(th_);
        axes(ax)
        plot(x_,y_,'-','Color',purp,'LineWidth',2);
        axes(axBIG)
        plot(x_,y_,'-','Color',purp,'LineWidth',2);hold on;
    elseif i==5
        xl=xlim;
        x_ = linspace(xl(1),-8,50); x_ = [x_ fliplr(x_)];
        y_ = linspace(0.5294*(-8-xl(1)),0,50); y_ = [y_ -fliplr(y_)];
        axes(ax)
        plot(x_,y_,'-','Color',purp,'LineWidth',2);
        axes(axBIG)
        plot(x_,y_,'-','Color',purp,'LineWidth',2);hold on;
    end
end

set(gcf,'Units','Normalized','Position',[0.2188 0.3630 0.4885 0.3083])
return