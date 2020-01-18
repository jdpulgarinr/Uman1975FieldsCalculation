%% Gráfica del campo eléctrico

close all
clear all
clc

tic

%% Constantes
c = 299792458;
v = 8e7;
u0 = 4*pi*(1e-7);
e0 = 8.8541878176*(1e-12);
H = 4000;

D = [1e3 10e3 100e3];

% Tiempo de simulación y particiones del tiempo
tf = (D(3)/c)+(100e-6);
nt = 2000;
t = 0:(tf/nt):tf;

% Matriz de elementos para almacenar los valores del campo eléctrico
E = zeros(length(D),(nt+1));

% Particiones para integración con respecto a z y vector de pesos para
% integral simple
n = 10000;
w = zeros(1,(n+1)) + 2;
w(1,1) = 1;
w(1,(n+1)) = 1;

% Particiones para integración con respecto al tiempo y vectores de pesos
% para integral doble
m = 1000;
Q = w;
P = zeros((m+1),1) + 2;
P(1,1) = 1;
P((m+1),1) = 1;

F = zeros((m+1),(n+1));

%% Integración
for s = 1:length(D)

    % Términos independientes del tiempo
    Z = 0:(H/n):H;
    R = sqrt((D(s)^2) + (Z.^2));
    TH = (atan(Z./D(s))) + (pi/2);
    
    % Términos independientes del tiempo en la integral de campo eléctrico
    T3 = ((2 - (3*((sin(TH)).^2)))./(R.^3)); 
    T4 = ((2 - (3*((sin(TH)).^2)))./(c*(R.^2)));
    T5 = (((sin(TH)).^2)./((c^2)*R));

    for k = 1:length(t)
        % Campo eléctrico
        [I,DI] = piecewiseCurrent(Z,R,t(k),v);
        
        TAU = 0:(t(k)/m):t(k);
        
        for d = 1:length(TAU)
            [Itau,Dtau] = piecewiseCurrent(Z,R,TAU(d),v);
            F(d,:) = T3.*Itau;
        end
        
        BB4 = T4.*(I);
        BB5 = T5.*(DI);
        E(s,k) = (1/(2*pi*e0))*((((H*t(k))/(4*m*n))*(Q*F'*P)) + (((H/(2*n))*(BB4*w'))) - (((H/(2*n))*(BB5*w'))));
    end
end

E = E*(-1e-2);
toc;

%% Gráfica del campo eléctrico
figure
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
plot(t,E(1,:),'LineWidth',2,'Color','k'), grid on, hold on;
plot(t,(E(2,:)*10),'LineWidth',2,'Color','b');
plot(t,(E(3,:)*100),'LineWidth',2,'Color','r');
xlabel('Time [s]','Interpreter','LaTeX','FontSize',30), ylabel('$||\vec{E}|| \times 10^{-2}$ Field [V/m]','Interpreter','LaTex','FontSize',30)
title('Electric Field','Interpreter','LaTeX','FontSize',30)
set(gca,'FontName','Times New Roman','FontSize',30)
axis([0 tf 1.2*min([min(E(1,:)) min(E(2,:)*10) min(E(3,:)*100)]) 1.2*max([max(E(1,:)) max(E(2,:)*10) max(E(3,:)*100)])])
lg = legend('$D = 1$ $km$ ($\times 1$)','$D = 10$ $km$ ($\times 10$)','$D = 100$ $km$ ($\times 100$)');
set(lg,'Interpreter','LaTeX')