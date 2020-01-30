function [currentVector,currentDerivative] = piecewiseCurrent(calculationHeight,radialDistance,calculationTime,speedRSFront)
%piecewiseCurrent This function calculates the current vector at a time
%instant and for different altitude points along a lightning channel, using
%a TL model and the channel base current of Figure 3a of the paper: The
%Electromagnetic Raadiation from a Finite Antenna, by Martin A. Uman, D.
%Kenneth McLain, E. Philip Krider. 1975. The function can also calculates
%the current vector at only one particular altitude and for a time vector.
%
% Syntax:  [current,currentDerivative] = piecewiseCurrent(ZZ,RR,TT,vv)
%
% Inputs:
%   calculationHeight:  A vector with the different heights along the
%   lightning channel at which the current is required for the time
%   instant calculationTime or a scalar with the height of the channel at
%   wich the current is required for the time vector calculationTime.
%   
%   radialDistance:     A vector with the radial distances from each point
%   of the calculationHeight vector to the observation point or a scalar
%   with the radial distance from the height calculationHeight to the
%   observation point.
%
%   calculationTime:    A scalar with the value of time at which the
%   current for the different points of calculationHeight is required or a
%   vector with the times at which current is required for the only height
%   of calculationHeight.
%
%   speedRSFront:       Speed of the return-stroke front.
%
% Outputs:
%    currentVector:     The current vector.
%
%    currentDerivative: The current vector derivative.
%
% Other m-files required: none.
% Subfunctions: none
% MAT-files required: none

% Author: Juan Pulgarin
% email: jdpulgarinr@gmail.com
% June, 2014.

%------------- BEGIN CODE --------------

    u0 = 4*pi*1e-7;
    e0 = 8.8541878176*1e-12;
    c = 1/sqrt(u0*e0);

    T1 = ((radialDistance/c) + (calculationHeight/speedRSFront));
    T2 = (T1 + (1e-6));
    T3 = (T1 + (25e-6));

    I1 = ((10e9*calculationTime) - ((10e9)*((radialDistance/c)+(calculationHeight/speedRSFront))));
    I2 = ((-5e9/12)*calculationTime) + ((5e9/12)*((radialDistance/c) + (calculationHeight/speedRSFront) + (1e-6))) + 10e3;

    M1 = 10e9;
    M2 = -5e9/12;

    currentVector = ((I1.*((calculationTime>=T1)&(calculationTime<=T2))) + (I2.*((calculationTime>T2)&(calculationTime<=T3))));
    currentDerivative = ((M1*((calculationTime>=T1)&(calculationTime<=T2))) + (M2*((calculationTime>T2)&(calculationTime<=T3))));
end