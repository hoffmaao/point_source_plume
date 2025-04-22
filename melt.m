%----------------------------------------------------------------------%
% MELT.M                                                 [Version 1.00]
% ******
%----------------------------------------------------------------------%
% Function to calculate melt rate and associated variables, given ISW
% properites and depth of ice shelf base.
%
% Input variables are:
%     D  : thickness of plume.
%     U  : velocity of plume.
%     T  : temperature of plume.
%     S  : salinity of plume.
%     De : depth of ice shelf base.
% Output variables are:
%     M  : melt rate (negative for freezing).
%     Tb : temperature at ice shelf base.
%     Sb : salinity at ice shelf base.
%     GT : heat transfer coefficient.
%     GS : salt transfer coefficient.
%     Tf : freezing point of ISW plume.
%
% Version 1.00
% Created 17/10/97 by A Jenkins.
%----------------------------------------------------------------------%

function [M,Tb,Sb,GT,GS,Tf] = melt(D,U,T,S,De,Tstar)

global K Ueps L Ci C FPa FPb FPc GamT GamS Ti Si

%----------------------------------------------------------------------%
% Put lower bound on U.
% This guards against taking the log of zero.

U = max([abs(U) Ueps]);

% Calculate turbulent exchange coefficients

GT = U*GamT;
GS = U*GamS;

%----------------------------------------------------------------------%
% Calculate the in situ freezing point.

Tf = FPa*S + FPb + FPc*De;

% If the ISW is above freezing melting will occur.
% Set a flag to indicate when this occurs.

MFlag = fix((sign(T-Tf)+1)/2);

% Calculate the "freezing point" of the ice.
% The ice is fresh if melting, but has salinity Si if freezing.

Tfi = FPa*(1-MFlag)*Si + FPb + FPc*De;

% Calculate melt rate.
% Heat conduction is only active if melting is occurring.

Q1 = L/C + MFlag*(Ci/C)*(Tfi-Ti);
Q2 = GS*(L/C + MFlag*(Ci/C)*(Tf-Ti)) + GT*(Tfi-T);
Q3 = GS*GT*(Tf-T);
M = -(Q2-sqrt((Q2^2)-4*Q1*Q3))/(2*Q1);

%M=.25*(Tstar-Tf).^(4/3)/1e6;
%if real(M)<real(Mstar)
%    M=Mstar;
%else
%end

% Calculate boundary temperature and salinity

Tb = (GT*T + MFlag*(Ci/C)*M*Ti - (L/C)*M) / (GT + MFlag*(Ci/C)*M);
Sb = (Tb - FPb - FPc*De)/FPa;
%----------------------------------------------------------------------%