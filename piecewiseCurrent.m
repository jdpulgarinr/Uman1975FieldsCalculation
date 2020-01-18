function [I,DI] = piecewiseCurrent(ZZ,RR,TT,vv)

c = 299792458;

T1 = ((RR/c) + (ZZ/vv));
T2 = (T1 + (1e-6));
T3 = (T1 + (25e-6));

I1 = ((10e9*TT) - ((10e9)*((RR/c)+(ZZ/vv))));
I2 = ((-5e9/12)*TT) + ((5e9/12)*((RR/c) + (ZZ/vv) + (1e-6))) + 10e3;

M1 = 10e9;
M2 = -5e9/12;

I = ((I1.*((TT>=T1)&(TT<=T2))) + (I2.*((TT>=T2)&(TT<=T3))));
DI = ((M1*((TT>=T1)&(TT<=T2))) + (M2*((TT>=T2)&(TT<=T3))));