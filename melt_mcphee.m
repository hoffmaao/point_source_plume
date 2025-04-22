%----------------------------------------------------------------------%
% MELT.M                                                 [Version 1.01]
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
% Version 1.01
% Created 17/10/97 by A Jenkins.
% Revised 08/08/05 by A Jenkins.
%     Added rms tidal current to friction velocity calculation.
%----------------------------------------------------------------------%

function [M,Tb,Sb,GT,GS,Tf] = melt_mcphee(D,U,T,S,De)

global Ueps L Ci C FPa FPb FPc GamTS Ti

%----------------------------------------------------------------------%
% Add in rms tidal current and put lower bound on U.
% This guards against taking the log of zero.

U = max([U Ueps]);

% Calculate turbulent exchange coefficients

GT = U*GamTS;
GS = U*GamTS;

%----------------------------------------------------------------------%
% If the ISW is above freezing melting will occur.
% Set a flag to indicate when this occurs.

Tf = FPa*S + FPb + FPc*De;
Tb = Tf;
Sb = S;

% If the ISW is above freezing melting will occur.
% Set a flag to indicate when this occurs.

MFlag = fix((sign(T-Tf)+1)/2);

% Calculate melt rate.
% Heat conduction is only active if melting is occurring.

M = GT*C*(T-Tf)/(L-MFlag*Ci*(Ti-Tf));

%----------------------------------------------------------------------%