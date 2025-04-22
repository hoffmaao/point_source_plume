%----------------------------------------------------------------------%
% ICEOCEAN.M                                             [Version 1.00]
% **********
%----------------------------------------------------------------------%
% Function to calculate along plume gradients of thickness, velocity,
% temperature and salinity for use in the ordinary differential
% equation solver ODEMOD.M.
%
% Input variables are:
%     Z    : independent variable (distance along ice shelf base).
%     X    : vector of dependent variables
%                 (X(1) is plume thickness,
%                  X(2) is plume velocity,
%                  X(3) is plume temperature,
%                  X(4) is plume salinity).
% Output variables are:
%     Xdot : vector of derivatives of dependent variables (dX/dZ).
%     Test : flag to indicate that integration should continue
%                 (switches off when plume is no longer buoyant).
%
% Version 1.00
% Created 17/10/97 by A Jenkins.
%----------------------------------------------------------------------%

function [Xdot,Test] = iceocean(Z,X)

global E0 K BS BT G C IDdata ISdata TSdata

%----------------------------------------------------------------------%
% Calculate depth and slope of ice shelf base.

De = interp1(IDdata(:,1),IDdata(:,2),Z);
SinT = interp1(ISdata(:,1),ISdata(:,2),Z);
CosT = sqrt(1-SinT^2);

%----------------------------------------------------------------------%
% Calculate ambient water properties at mid-depth of plume.

Tstar= trapz(TSdata(:,1),TSdata(:,4))/(TSdata(end,1)-TSdata(1,1));

TS = interp1(TSdata(:,1),TSdata(:,2:3),...
             max([De-abs(X(1))*CosT/2; TSdata(1,1)]));
Ta = TS(1,1);
Sa = TS(1,2);

% Calculate density contrast between plume and ambient water using a
% linear equation of state.

Rho = BS*(Sa-X(4))-BT*(Ta-X(3));

%----------------------------------------------------------------------%
% Calculate melt rate and associated variables using function MELT.M.

[M,Tb,Sb,GT,GS,Tf] = melt(X(1),X(2),X(3),X(4),De,Tstar);
% [M,Tb,Sb,GT,GS,Tf] = melt_mcphee(X(1),X(2),X(3),X(4),De);

%----------------------------------------------------------------------%
% Calculate ambient water properties at bottom of plume.

TS = interp1(TSdata(:,1),TSdata(:,2:3),...
             max([De-abs(X(1))*CosT; TSdata(1,1)]));
Ta = TS(1,1);
Sa = TS(1,2);

% Calculate entrainment rate based on Ellison and Turner (1959).

E = E0*X(2)*SinT;

%----------------------------------------------------------------------%
% Solve differential equations for D, U, T and S

Xdot(1) = 2*E/X(2) + 4*M/(pi*X(2)) - X(1)*G*SinT*Rho/(2*X(2)^2)+2*K/pi; %K + 2*(E+M)/X(2) - Rho*G*SinT*X(1)/(X(2)^2);
Xdot(2) = -2*E/X(1) - 4*M/(pi*X(1)) + G*SinT*Rho/X(2) - 4*K*X(2)/(pi*X(1)); %Rho*G*SinT/X(2) - K*X(2)/X(1) - (E+M)/X(1);
Xdot(3) = 2*E/X(2)*(Ta-X(3))/X(1) + 4*M*(Tb-X(3))/(pi*X(1)*X(2)) - 4*GT*sqrt(K)*(X(3)-Tb)/(pi*X(1)); %(Ta-X(3))*E/(X(1)*X(2)) + (Tb-X(3))*(GT+M)/(X(1)*X(2));
Xdot(4) = 2*E/X(2)*(Sa-X(4))/X(1) + 4*M*(Sb-X(4))/(pi*X(1)*X(2)) - 4*GS*sqrt(K)*(X(4)-Sb)/(pi*X(1)); %(Sa-X(4))*E/(X(1)*X(2)) + (Sb-X(4))*(GS+M)/(X(1)*X(2));

%----------------------------------------------------------------------%
% Check plume is still buoyant, and set flag to zero if it isn't.
% Do not do this if Rho is complex (step will be recomputed anyway).

% if isreal(Rho) & isreal(X(2))
%     Test = max([(sign(Rho)+1)/2 1-(sign(X(2))+1)/2]);
if isreal(Rho)
    Test = (sign(Rho)+1)/2;
else
    Test = 1;
end

%----------------------------------------------------------------------%