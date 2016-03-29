%% roterande skålen
clc; clf; clear
eps=1e-6;
T=[eps 0.7875 0.6851 0.6337 0.6048 0.7079 0.9691 0.8268 0.7488 0.7452 0.8202 1.0129 0.8815 0.8082 0.6781 0.5795 0.5791 0.7750 0.6822 0.5796 0.6567];
h=0*4e-3+[0 75 57 50 43 63 104 67 62 63 76 129 89 72 52 37 38 75 56 41 53]*1e-3;
hold on
plot(T.^2/(pi^2*8), h,'.')
g=h./T.^2*8*pi^2
medel=mean(g)
standerror=std(g)/sqrt(length(T))
f = @(k) sum((h-k*T.^2/(pi^2*8)).^2);
g=fminbnd(f,9,10)

%% korrektion
clc; clf; clear
hold on
rho=5e-3;

r1=7.5e-2;
r2=6.3e-2;

T1=[0.7875 0.6851 0.6337 0.6048 0.7079 0.9691  0.7488 0.7452 0.8202 1.0129 0.8815 0.8082 0.6781 0.5795 0.5791];
%T1=[0.7875 0.6851 0.6337 0.6048 0.7079 0.9691  0.7488 0.7452 0.8202 1.0129 0.8815 0.8082 0.6781];
T2=[0.7750 0.6822 0.5796 0.6567];
T=[T1, T2];

w1=2*pi./T1;
w2=2*pi./T2;

h1=[75 57 50 43 63 104  62 63 76 129 89 72 52 37 38]*1e-3;
%h1=[75 57 50 43 63 104  62 63 76 129 89 72 52]*1e-3;
h2=[75 56 41 53]*1e-3;

h=[h1,h2];
r=[r1*ones(1,length(h1)), r2*ones(1,length(h2))];
W=[w1, w2];

%f1=@(v) cot(v)+r/(2*(r-rho))*tan(v/2)-h./(r-rho);
alfa=zeros(1,length(h1));
alfa2=zeros(1,length(h1));
for i=1:length(h)
    f=@(v) cot(v)+r(i)/(2*(r(i)-rho))*tan(v/2)-h(i)./(r(i)-rho);
    f2=@(v) cot(v)+r(i)/(2*r(i))*tan(v/2)-h(i)./r(i);
    
    alfa(i)=fsolve(f,atan(h(i)/r(i)));clc
    alfa2(i)=fsolve(f2,atan(h(i)/r(i)));clc
end


delta=rho*cot(alfa);
max(delta)*1e3
min(delta)*1e3

H=h+delta;

omega_kvadrat=T.^2/(pi^2*8);

g=H./omega_kvadrat;
g1=mean(g);
fprintf('Medelvärdet av g är: \t %0.3f\n', g1)
fprintf('Standfelet är: \t\t %0.3f \n', std(g)/sqrt(length(T)))

f = @(k) sum((H-k*T.^2/(pi^2*8)).^2);
g2=fminbnd(f,9,10);
fprintf('g är med minsta kvadrat: %0.3f \n\n', g2)

[g3, std_g]=lscov(omega_kvadrat.', H.');

x=linspace(0,.016);
hold on
%plot(x,g1*x);
err1=delta(1:length(T1));
err2=delta((1+length(T1)):end);

%errorbar(omega_kvadrat(1:length(T1)) , H(1:length(T1))*1e3, err1*1e3, '.', 'markersize',16)
%errorbar(omega_kvadrat((1+length(T1)):end) , H((1+length(T1)):end)*1e3, err2*1e3, '*r', 'markersize',8)

plot(omega_kvadrat(1:length(T1)) , H(1:length(T1))*1e3, '.', 'markersize',16)
plot(omega_kvadrat((1+length(T1)):end) , H((1+length(T1)):end)*1e3, '*r', 'markersize',8)

p=plot(x,x*g2*1e3, 'k');
axis([0, 0.016, 0, 160])
%grid on

xtext='$1/(2\omega^{2})\:/ [\mathrm{s}^{2}]$';
ytext='$h\:/[\mathrm{mm}]$';

xlabel(xtext, 'Interpreter', 'Latex', 'FontSize', 20, 'Color', 'k');
ylabel(ytext, 'Interpreter', 'Latex', 'FontSize', 20, 'Color', 'k');
set(gca,'FontSize',15, 'xtick',(0:4:16)*1e-3);


l=legend('$r=7{,}4$\,cm', '$r=6{,}3$\,cm', '$g/(2\omega^2)$: $g=9.78$\,m\,s$^{-2}$');
set(l,'interpreter', 'latex', 'fontsize', 18, 'location', 'southeast')































