%%
clear;clc;close all;fclose all;
%%
%定义常值
q = 1e8; %集度 10N/m
l = 1;  %长度 1m
E = 200e9; %杨氏模量 Pa
I = 1e1; %一次矩
n = 5; %阶数

vMaxComp = -5/384*q*l^4/(E*I);
syms xVariable vVariable;
xValue = 0:l/100:l;
vValue = ones(1,length(xValue));
vMax = zeros(1,n);

figure(1);
for i = 1 : n
    if (i == 1)
        vVariable = -4*q*l^4/(E*I*pi^5) * 1/i^5*sin(i*pi*xVariable/l);
    else
        vVariable = vVariable - 4*q*l^4/(E*I*pi^5) * 1/i^5*sin(i*pi*xVariable/l);
    end
   vValue = double(subs(vVariable,xVariable,xValue));
   axes1 = subplot(2,2,1:2); plot(xValue,vValue,'LineWidth',1.5); grid on; hold on;
   legend1 = legend(axes1,'show');
   set(legend1,'Position',[0.129880951597578 0.583478675333497 0.123928572211946 0.176190480050586]);
   vMax(i) = subs(vVariable,xVariable,l/2);
end

subplot(2,2,1:2);legend('n=1','n=2','n=3','n=4','n=5');
xlabel('Length-m');ylabel('u_y/m');title('u_y - Length');

vMaxError = (vMax - vMaxComp * ones(1,length(vMax)) )/ vMaxComp;

subplot(2,2,3);scatter(1:n,vMax,35,'MarkerEdgeColor',[0.0 0.0 0.0],'MarkerFaceColor',[0 .5 .5]);grid on;box on;
xlabel('Order-1');ylabel('u_y_ _m_a_x/m');title('u_y_ _m_a_x - Order');
subplot(2,2,4);plot(1:n,vMaxError * 100,'LineWidth',1.5);hold on; scatter(1:n,vMaxError * 100,'LineWidth',1.0); grid on; box on;
xlabel('Order-1');ylabel('Relative Error - %');title('Relative Error - Order');
