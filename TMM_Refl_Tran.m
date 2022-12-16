syms r t k0 k1 k2 k3 k4 k5 k6 k7 k8 k9 d1 d2 d3 d4 d5 d6 d7 d8 d9;

% parameter definition
lamda=(100:0.01:400)*1e-9;
c=2.99792458e8;
nu=c./lamda;
omiga=2*pi.*nu;

ng=1.5395;
n1=2.8;
n2=1.43;
npr=1.3078;

kg=2*pi*nu*ng/c;
k1=2*pi*nu*n1/c;
k2=2*pi*nu*n2/c;
k3=k1;
k4=k2;
k5=k1;
k6=k2;
k7=k1;
k8=k2;
k9=k1;
kpr=2*pi*nu*npr/c;

dg=2e-9;
d1=17.2e-9;
d2=33.7e-9;
d3=d1;
d4=d2;
d5=d1;
d6=d2;
d7=d1;
d8=d2;
d9=d1;
dpr=100e-9; % PR²ãµÄ²ÎÊý

r=zeros(1,30001);
t=zeros(1,30001);

for i=1:30001
    
% transfer matrix mothod
% Ma=[cos(k0(1,i).*d0) sin(k0(1,i).*d0)./k0(1,i); -k0(1,i).*sin(k0(1,i).*d0) cos(k0(1,i).*d0)];
M0=[cos(kg(1,i).*dg) sin(kg(1,i).*dg)./kg(1,i); -kg(1,i).*sin(kg(1,i).*dg) cos(kg(1,i).*dg)];
M1=[cos(k1(1,i).*d1) sin(k1(1,i).*d1)./k1(1,i); -k1(1,i).*sin(k1(1,i).*d1) cos(k1(1,i).*d1)];
M2=[cos(k2(1,i).*d2) sin(k2(1,i).*d2)./k2(1,i); -k2(1,i).*sin(k2(1,i).*d2) cos(k2(1,i).*d2)];
M3=[cos(k3(1,i).*d3) sin(k3(1,i).*d3)./k3(1,i); -k3(1,i).*sin(k3(1,i).*d3) cos(k3(1,i).*d3)];
M4=[cos(k4(1,i).*d4) sin(k4(1,i).*d4)./k4(1,i); -k4(1,i).*sin(k4(1,i).*d4) cos(k4(1,i).*d4)];
M5=[cos(k5(1,i).*d5) sin(k5(1,i).*d5)./k5(1,i); -k5(1,i).*sin(k5(1,i).*d5) cos(k5(1,i).*d5)];
M6=[cos(k6(1,i).*d6) sin(k6(1,i).*d6)./k6(1,i); -k6(1,i).*sin(k6(1,i).*d6) cos(k6(1,i).*d6)];
M7=[cos(k7(1,i).*d7) sin(k7(1,i).*d7)./k7(1,i); -k7(1,i).*sin(k7(1,i).*d7) cos(k7(1,i).*d7)];
M8=[cos(k8(1,i).*d8) sin(k8(1,i).*d8)./k8(1,i); -k8(1,i).*sin(k8(1,i).*d8) cos(k8(1,i).*d8)];
M9=[cos(k9(1,i).*d9) sin(k9(1,i).*d9)./k9(1,i); -k9(1,i).*sin(k9(1,i).*d9) cos(k9(1,i).*d9)];
M10=[cos(kpr(1,i).*dpr) sin(kpr(1,i).*dpr)./kpr(1,i); -kpr(1,i).*sin(kpr(1,i).*dpr) cos(kpr(1,i).*dpr)];
M11=[cos(kg(1,i).*dg) sin(kg(1,i).*dg)./kg(1,i); -kg(1,i).*sin(kg(1,i).*dg) cos(kg(1,i).*dg)];
% M=M0*M1*M2*M3*M4*M5*M6*M7*M8*M9*M10*M11;
M=M11*M10*M9*M8*M7*M6*M5*M4*M3*M2*M1*M0;

% r(1,i)=1-((2.*(M(2,1)+1i*k0(1,i).*M(2,2)))./(M(2,1)+1i*k0(1,i).*M(2,2)+1i*k0(1,i).*(M(1,1)+1i*k0(1,1).*M(1,2))));
r(1,i)=(-1i*kg(1,i)*M(1,1)+kg(1,i)^2*M(1,2)+M(2,1)+1i*kg(1,i)*M(2,2))/(1i*kg(1,i)*M(1,1)+kg(1,i)^2*M(1,2)-M(2,1)+1i*kg(1,i)*M(2,2));
% t(1,i) = (1+r(1,i))./(M(1,1)+1i*k0(1,i)*M(1,2));
t(1,i)=M(1,1)*(1+r(1,i))+M(1,2)*1i*kg(1,i)*(1-r(1,i));
end

R=abs(r);
T=abs(t);
plot(lamda*1e9, R, 'r', 'linewidth', 1.5);
hold on
plot(lamda*1e9, T, 'b', 'linewidth', 1.5);
hold off

set(gca, 'fontname','Times New Roman','fontsize', 12);
yyaxis left
xlabel('Wavelength(nm)');
ylabel('Reflectivity(a.u.)');

yyaxis right
ylabel('Transmittance(a.u.)');