%% Gr?fica del campo magn?tico

close all
clear all
clc

%% Constantes
u0 = 4*pi*1e-7;
e0 = 8.8541878176*1e-12;
c = 1/sqrt(u0*e0);

v=8e7;

H = 4000;

D = 1e3;

% Tiempo de simulaci?n y particiones del tiempo
t=linspace(0,90e-6,1000);

% Matriz de elementos para almacenar los valores del campo magn?tico
B = zeros(1,length(t));
    
% Particiones para integraci?n con respecto a z y vector de pesos para
% integral simple
n = 100000;
w = zeros(1,(n+1)) + 2;
w(1,1) = 1;
w(1,(n+1)) = 1;

%% Integraci?n

% T?rminos independientes del tiempo
Z = 0:(H/n):H;
R = sqrt((D^2) + (Z.^2));
TH = (atan(Z./D)) + (pi/2);

% T?rminos independientes del tiempo en la integral de campo magn?tico
T1 = (sin(TH))./(R.^2);
T2 = (sin(TH))./(c*R);

for k = 1:length(t)
    % Campo magn?tico
    [I,DI] = piecewiseCurrent(Z,R,t(k),v);
    BB1 = T1.*(I); 
    BB2 = T2.*(DI);
    B(1,k) = ((u0/(2*pi))*((H/(2*n))*(BB1*w'))) + ((u0/(2*pi))*((H/(2*n))*(BB2*w')));
end

%% Gr?fica del campo magn?tico
figure
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
plot(t,B,'LineWidth',2,'Color','k'), grid on, hold on;
%plot(t,(B(2,:)*10),'LineWidth',2,'Color','b');
%plot(t,(B(3,:)*100),'LineWidth',2,'Color','r');
%xlabel('Time [s]','Interpreter','LaTeX','FontSize',30), ylabel('$||\vec{B}|| \times 10^7$ Field [A/m]','Interpreter','LaTex','FontSize',30)
title('Magnetic Field','Interpreter','LaTeX','FontSize',30)
%set(gca,'FontName','Times New Roman','FontSize',30)
%axis([0 tf 1.2*min([min(B(1,:)) min(B(2,:)*10) min(B(3,:)*100)]) 1.2*max([max(B(1,:)) max(B(2,:)*10) max(B(3,:)*100)])])
%lg = legend('$D = 1$ $km$ ($\times 1$)','$D = 10$ $km$ ($\times 10$)','$D = 100$ $km$ ($\times 100$)');
%set(lg,'Interpreter','LaTeX')