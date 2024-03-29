function [timeVector, BFieldVector] = BField(calculationDistance)
%BField This function calculates the magnetic field at a distant point from
%a lightning channel by means of numerical integration of Equation (7) of
%the paper: The Electromagnetic Radiation from a Finite Antenna, by Martin
%A. Uman, D. Kenneth McLain, E. Philip Krider. 1975.
%
% Syntax:  [timeVector, BFieldVector] = BField(calculationDistance)
%
% Inputs:
%    calculationDistance:   This is the horizontal distance from the
%    lightning channel to the observation point of the field.
%
% Outputs:
%    timeVector:     The time vector used for the calculation of the field.
%    BFieldVector:   The magnetic field vector calculated at each element of
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

    % Magnetic Field Array Initialization
    B = zeros(1,length(timeVector));

    %% INTEGRATION

    % Channel Partitions and Weight Vector for Integration
    n = 100000;
    w = zeros(1,(n+1)) + 2;
    w(1,1) = 1;
    w(1,(n+1)) = 1;

    % Time-Independent Terms
    Z = 0:(H/n):H;
    R = sqrt((D^2) + (Z.^2));
    TH = (atan(Z./D)) + (pi/2);

    % Time-Independent Terms in Integral
    T1 = (sin(TH))./(R.^2);
    T2 = (sin(TH))./(c*R);

    % Time-Marching Integration
    for k = 1:length(timeVector)
        [I,DI] = piecewiseCurrent(Z,R,timeVector(k),v);
        BB1 = T1.*(I); 
        BB2 = T2.*(DI);
        B(1,k) = ((u0/(2*pi))*((H/(2*n))*(BB1*w'))) + ((u0/(2*pi))*((H/(2*n))*(BB2*w')));
    end  
    
    BFieldVector = B;
end