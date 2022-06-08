% %%
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
vValue = zeros(1,length(xValue));
vMaxValue = zeros(1,n);
%%
%计算挠度的多项式表达式
syms cTemp;

for i = 1: n
    vValue = zeros(1,length(xValue));
    eval(['syms c',num2str(i)]);
    eval(['cTemp = c',num2str(i)]);
    if (i == 1)
        vVariable = cTemp * xVariable * (l-xVariable);
%         vVariable = cTemp * xVariable * (l-xVariable);
    else
        vVariable = vVariable + cTemp * xVariable^i * (l^i - xVariable^i);
%         vVariable = vVariable + cTemp * xVariable^i * (l - xVariable)^i;
    end
    e1 = E*I* diff(vVariable,xVariable,2)^2/2;
    e2 = q * vVariable;
    eTotal = int(e1+e2,xVariable,[0,l]);
    Equs = [];
    for j = 1:i
        eval(['equ',num2str(j),' = diff(eTotal,c',num2str(j),') == 0;']);
        eval(['Equs = [Equs;equ',num2str(j),']']);
    end
    cValueStruct = solve(Equs);
    for j = 1:i
        if (i == 1)
            eval(['cValue(',num2str(j),') = double(cValueStruct)']);
        else
            eval(['cValue(',num2str(j),') = double(cValueStruct.c',num2str(i),')']);
        end
        vValue = vValue + cValue(j).*xValue.^j.*(l^i-xValue.^j); 
    end
    vMaxValue(i) = vValue(50);
    axes1 = subplot(2,2,1:2);
    plot(xValue,vValue,'LineWidth',1.5); grid on; hold on;
end
subplot(2,2,1:2);
xlabel('Length-m');ylabel('u_y/m');title('u_y - Length');
legend('n=1','n=2','n=3','n=4','n=5');
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.129880951597578 0.583478675333497 0.123928572211946 0.176190480050586]);

vMaxError = (vMaxValue - vMaxComp * ones(1,length(vMaxValue)) )./ vMaxComp;

subplot(2,2,3);scatter(1:n,vMaxValue,35,'MarkerEdgeColor',[0.0 0.0 0.0],'MarkerFaceColor',[0 .5 .5]);grid on;box on;
xlabel('Order-1');ylabel('u_y_ _m_a_x/m');title('u_y_ _m_a_x - Order');
subplot(2,2,4);plot(1:n,vMaxError * 100,'LineWidth',1.5);hold on; scatter(1:n,vMaxError * 100,'LineWidth',1.0); grid on; box on;
xlabel('Order-1');ylabel('Relative Error - %');title('Relative Error - Order');
