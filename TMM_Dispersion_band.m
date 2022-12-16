clear all;
clc;

s=2.5;

lamda0=193e-9;
c=2.99792458e8;
omiga0=2*pi*c/lamda0;
afreq=(0:0.01:4)*omiga0;

n0=1;
n1=2.8;
n2=1.43;

d1=lamda0/(4*n1);
d2=lamda0/(4*n2);

TEkd=[];
TMkd=[];

for i=1:401
    for j=1:401
        
        omiga=afreq(1,j);
        k0=omiga*n0/c;
        kx=0.01*(i-1)*k0;
%         sst0=kx/k0;
%         sst1=n0*sst0/n1;
%         cst1=sqrt(1-sst1^2);
%         sst2=n0*sst0/n2;
%         cst2=sqrt(1-sst2^2);
        k1=sqrt((n1*omiga/c)^2-kx^2);
        k2=sqrt((n2*omiga/c)^2-kx^2);
%         k1=omiga*n1.*cst1/(c);
%         k2=omiga*n2.*cst2/(c);
        TEkd(i,j)=cos(k1*d1)*cos(k2*d2)-0.5*(k2/k1+k1/k2)*sin(k1*d1)*sin(k2*d2);      
        TMkd(i,j)=cos(k1*d1)*cos(k2*d2)-0.5*(n2^2*k1/(n1^2*k2)+n1^2*k2/(n2^2*k1))*sin(k1*d1)*sin(k2*d2);        
        
    end
end

[x,y]=meshgrid(afreq/omiga0, afreq/omiga0);
L1=abs(TEkd')<=1;
L2=abs(TMkd')<=1;

yte=y(L1);
ytm=y(L2);
xte = x(L1); 
xtm = x(L2);

% Tekd=abs(TEkd);
% mesh(Tekd)
scatter(xte, yte,s, 'b','filled');
hold on
% figure
scatter(-xtm, ytm,s, 'b','filled');
hold off
set(gca, 'fontname','Times New Roman','fontsize', 12);
title('Dispersion Band'); 
xlabel('{\itk_x}/{\itk}_0'); 
ylabel('\omega/\omega_0');
grid on
box on



