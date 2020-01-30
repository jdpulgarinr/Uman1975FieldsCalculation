function [timeVector, EFieldVector] = EField(calculationDistance)
%EField This function calculates the electric field at a distant point from
%a lightning channel by means of numerical integration of Equation (9) of
%the paper: The Electromagnetic Radiation from a Finite Antenna, by Martin
%A. Uman, D. Kenneth McLain, E. Philip Krider. 1975.
%
% Syntax:  [timeVector, EFieldVector] = EField(calculationDistance)
%
% Inputs:
%    calculationDistance:   This is the horizontal distance from the
%    lightning channel to the observation point of the field.
%
% Outputs:
%    timeVector:     The time vector used for the calculation of the field.
%
%    EFieldVector:   The electric field vector calculated at each element of
%    the time vector.
%
% Other m-files required: piecewiseCurrent.m
% Subfunctions: none
% MAT-files required: none

% Author: Juan Pulgarin
% email: jdpulgarinr@gmail.com
% June, 2014.

%------------- BEGIN CODE --------------

%% CONSTANTS
    % Empty Space
    u0 = 4*pi*1e-7;
    e0 = 8.8541878176*1e-12;
    c = 1/sqrt(u0*e0);

    % Speed of the Return-Stroke Front
    v=8e7;

    % Channel Height
    H = 4000;

    % Observation Distance from Channel
    D = calculationDistance;

%% VECTORS INITIALIZATION
    % Time Vector
    timeVector=linspace(0,90e-6,1000);

    % Electric Field Array Initialization
    E = zeros(1,length(timeVector));
    
%% INTEGRATION
    % Channel Partitions and Weights Vector for Single Integration
    n = 10000;
    w = zeros(1,(n+1)) + 2;
    w(1,1) = 1;
    w(1,(n+1)) = 1;

    % Time Partitions and Wights Vector for Double Integration
    m = 1000;
    Q = w;
    P = zeros((m+1),1) + 2;
    P(1,1) = 1;
    P((m+1),1) = 1;

    % Function Evaluations Array Initialization
    F = zeros((m+1),(n+1));

    % Time-Independent Terms
    Z = 0:(H/n):H;
    R = sqrt((D^2) + (Z.^2));
    TH = (atan(Z./D)) + (pi/2);
    
    % Time-Independent Terms in the Electric Field Integral
    T3 = ((2 - (3*((sin(TH)).^2)))./(R.^3)); 
    T4 = ((2 - (3*((sin(TH)).^2)))./(c*(R.^2)));
    T5 = (((sin(TH)).^2)./((c^2)*R));

    % Time-Marching Loop for Integration
    for k = 1:length(timeVector)
        [I,DI] = piecewiseCurrent(Z,R,timeVector(k),v);

        TAU = 0:(timeVector(k)/m):timeVector(k);

        for d = 1:length(TAU)
            [Itau,~] = piecewiseCurrent(Z,R,TAU(d),v);
            F(d,:) = T3.*Itau;
        end

        BB4 = T4.*(I);
        BB5 = T5.*(DI);
        E(1,k) = (1/(2*pi*e0))*((((H*timeVector(k))/(4*m*n))*(Q*F'*P)) + (((H/(2*n))*(BB4*w'))) - (((H/(2*n))*(BB5*w'))));
    end
    EFieldVector = E;
end