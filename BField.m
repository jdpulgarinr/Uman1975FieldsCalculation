%% Gr?fica del campo magn?tico

close all
clear all
clc

tic

%% Constantes
c = 299792458;
v = 8e7;
u0 = 4*pi*(1e-7);
e0 = 8.8541878176e-12;
H = 4000;

D = [1e3 10e3 100e3];

% Tiempo de simulaci?n y particiones del tiempo
tf = (D(3)/c)+(100e-6);
nt = 2000;
t = 0:(tf/nt):tf;

% Matriz de elementos para almacenar los valores del campo magn?tico
B = zeros(length(D),(nt+1));

% Matriz de elementos para almacenar los valores del campo el?ctrico
E = B;
    
% Particiones para integraci?n con respecto a z y vector de pesos para
% integral simple
n = 100000;
w = zeros(1,(n+1)) + 2;
w(1,1) = 1;
w(1,(n+1)) = 1;

%% Integraci?n
for s = 1:length(D)

    % T?rminos independientes del tiempo
    Z = 0:(H/n):H;
    R = sqrt((D(s)^2) + (Z.^2));
    TH = (atan(Z./D(s))) + (pi/2);
    
    % T?rminos independientes del tiempo en la integral de campo magn?tico
    T1 = (sin(TH))./(R.^2);
    T2 = (sin(TH))./(c*R);

    for k = 1:length(t)
        % Campo magn?tico
        [I,DI] = piecewiseCurrent(Z,R,t(k),v);
        BB1 = T1.*(I); 
        BB2 = T2.*(DI);
        B(s,k) = ((u0/(2*pi))*((H/(2*n))*(BB1*w'))) + ((u0/(2*pi))*((H/(2*n))*(BB2*w')));
    end
end

B = B*(1e7);

toc;

%% Gr?fica del campo magn?tico
figure
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
plot(t,B(1,:),'LineWidth',2,'Color','k'), grid on, hold on;
plot(t,(B(2,:)*10),'LineWidth',2,'Color','b');
plot(t,(B(3,:)*100),'LineWidth',2,'Color','r');
xlabel('Time [s]','Interpreter','LaTeX','FontSize',30), ylabel('$||\vec{B}|| \times 10^7$ Field [A/m]','Interpreter','LaTex','FontSize',30)
title('Magnetic Field','Interpreter','LaTeX','FontSize',30)
set(gca,'FontName','Times New Roman','FontSize',30)
axis([0 tf 1.2*min([min(B(1,:)) min(B(2,:)*10) min(B(3,:)*100)]) 1.2*max([max(B(1,:)) max(B(2,:)*10) max(B(3,:)*100)])])
lg = legend('$D = 1$ $km$ ($\times 1$)','$D = 10$ $km$ ($\times 10$)','$D = 100$ $km$ ($\times 100$)');
set(lg,'Interpreter','LaTeX')