%----------------------------------------------------------------------%
% PLUME.M                                                [Version 2.00]
% *******
%----------------------------------------------------------------------%
% Script file to solve ordinary differential equations (ODE's) for an
% ISW plume.
%
% User is prompted for two parameters to define the intial conditions.
%      DU : initial flux of ice shelf water in m2/s.
%      cf : determines the initial temperature and salinity of the
%           plume, which must lie somewhere between the ambient
%           temperature and salinity and the freezing point (values
%           obtained by setting cf to 0 and 1 respectively).
%
% Results are contained in the following arrays:
%      X    : distance along plume path (m).
%      D    : plume thickness (m).
%      U    : plume velocity (m/s).
%      T    : plume temperature (deg C).
%      S    : plume salinity (psu).
%      De   : elevation (negative) of ice shelf base (m).
%      SinT : slope of ice shelf base.
%      Ta   : ambient temperature (deg C).
%      Sa   : ambient salinity (psu).
%      Rho  : dimensionless density contrast.
%      M    : melt rate (m/s converted to m/yr for output).
%      Tb   : temperature at ice shelf base (deg C).
%      Sb   : salinity at ice shelf base (psu).
%      GT   : heat exchange coefficient (m/s).
%      GS   : salt exchange coefficient (m/s).
%      Tf   : in situ freezing point of ISW plume (deg C).
%      E    : entrainment rate (m/s).
%
% Version 2.00
% [Compatible with MATLAB versions 6 and 7]
% Created 23/10/97 by A Jenkins.
% Initial inflow now freshwater, cooled/freshened ambient water now
% only used for restart plumes, where cooling fraction is set
% theoretically rather than chosen.
%
% Version 1.00
% [Compatible with MATLAB versions 4.2 and 5.0]
% Created 23/10/97 by A Jenkins.
%----------------------------------------------------------------------%

% Declare data variables as global.

global IDdata ISdata TSdata GamT GamS

%----------------------------------------------------------------------%

% Clear variables that are not overwitten.

X = [];
D = [];
U = [];
T = [];
S = [];

%----------------------------------------------------------------------%
% Define physical constants and make them global variables.
% To do this call the script file DEFINE.M.

define

%----------------------------------------------------------------------%
% Load ice shelf base and ambient water column data.
% To do this call the script file DATALOAD.M.

dataload

%----------------------------------------------------------------------%
% Prepare for integration.

% Prompt user for initial conditions.
% If values do not exist set defaults.

if ~exist('DUf')
   DUf = 0.0;
end

% Check if user wishes to use current values.

fprintf('\n');
disp(['Initial freshwater flux currently set to ',num2str(DUf),' m2/s.']);
newDUf = input('   To change this enter a new value:   ');
if ~isempty(newDUf)
   DUf = newDUf;
end

%----------------------------------------------------------------------%
% Check if user wishes to increase the ambient temperature.

Tstar = trapz(TSdata(:,1),TSdata(:,4))/(TSdata(end,1)-TSdata(1,1));

fprintf('\n');
disp(['Ambient water currently ',num2str(Tstar),' deg above '...
      'surface freezing point.']);
Tinc = input('   To change this enter a temperature increment:   ');
if isempty(Tinc)
   Tinc = 0;
end

TSdata(:,2) = TSdata(:,2) + Tinc;
TSdata(:,4) = TSdata(:,4) + Tinc;

%----------------------------------------------------------------------%
% Check if user wishes to generate repeated plumes until the end of the
% domain is reached.  Default is NO.

fprintf('\n');
disp('Continue integration to end of domain using multiple plumes?');
rep = input('   [Default=N]:   ','s');
if isempty(rep)
   rep = 'N';
end
fprintf('\n');

%----------------------------------------------------------------------%
% Define limits of domain.

Zf = IDdata(size(IDdata,1),1);
Z = 0;

%----------------------------------------------------------------------%
%----------------------------------------------------------------------%
% Enter solution loop one or multiple times.

IP = 0;

while (Z(size(Z,1)) < Zf) & (IP < 1 | rep == 'Y' | rep == 'y')

   IP = IP+1;

%----------------------------------------------------------------------%
% Set initial conditions for this plume.

% Set starting point and look up depth of ice shelf base, slope of ice
% shelf base and ambient properties at this point.

   Z0 = Z(size(Z,1));
   De = interp1(IDdata(:,1),IDdata(:,2),Z0);
   SinT = interp1(ISdata(:,1),ISdata(:,2),Z0);
   TS = interp1(TSdata(:,1),TSdata(:,2:3),De);
   Ta = TS(1,1);
   Sa = TS(1,2);
   
% Set the initial temperature and salinity by cooling and freshening
% the ambient water along the line defined by Gade (1979).  The factor
% "cf" determines how close to the in situ freezing point to go.

   Q1 = FPa*(1-Ci/C);
   Q2 = (FPb+FPc*De)*(1-Ci/C)-Ta-L/C+(Ci/C)*(FPa*Sa+Ti);
   Q3 = Sa*(L/C+(Ci/C)*(FPb+FPc*De-Ti));
   S0 = -(Q2+sqrt((Q2^2)-4*Q1*Q3))/(2*Q1);
   T0 = FPa*S0+FPb+FPc*De;
   
% Calculate the cooling fraction as the ratio of Stanton number to
% entrainment constant.

   Stant = GamT/(1-(FPa*S0*C*GamT)/(L*GamS));
   cf = Stant/(E0*SinT+Stant);
   
   if IP > 1

% Calculate the initial flux, temperature and salinity.
   
       DU0 = DU;
       S0 = Sa + cf*(S0-Sa);
       T0 = Ta + cf*(T0-Ta);
       
   else
       
% Set the initial salinity to zero and the intitial temperture to
% the pressure freezing point.

       DU0 = DU + DUf;
       S0 = DU * (Sa + cf*(S0-Sa)) / DU0;
       T0 = (DU * (Ta + cf*(T0-Ta)) + DUf * (FPb+FPc*De)) / DU0;

   end
   
% Set initial velocity and thickness.  This is done by assuming a
% balance between buoyancy and friction at the start.

   Rho = BS*(Sa-S0)-BT*(Ta-T0);
   U0 = max([(Rho*G*SinT*DU/(E0*SinT+K))^0.333 Ueps]);
   D0 = DU0/U0;

%----------------------------------------------------------------------%
% Form the vector of initial conditions used by the ODE solver.

   Y0 = [D0 U0 T0 S0]';

% Solve ODE's using the Runge-Kutta solver in ODEMOD.M.
% The derivatives are evaluated by the function ICEOCEAN.M.

   [Z,Y] = odemod('iceocean',Z0,Zf,Y0);

% Inform the user of progress.

   fprintf('\n');
   disp(['Plume ',num2str(IP),'; endpoint ',num2str(Z(size(Z,1))),...
         ' m.']);

% If the plume failed to get anywhere try again with increased initial
% flux.

   if Z(size(Z,1)) == 0
       IP = IP-1;
       DU = DU*10;
   else

%----------------------------------------------------------------------%
% Add the solution on to those from earlier plumes.

       X = [X;Z];
       D = [D;Y(:,1)];
       U = [U;Y(:,2)];
       T = [T;Y(:,3)];
       S = [S;Y(:,4)];
       
   end

% End of while loop.

end

%----------------------------------------------------------------------%
%----------------------------------------------------------------------%
% Using the solution for plume thickness, velocity, temperature and
% salinity, recalculate other important variables.

% Create empty arrays for results that may be plotted.

Res = zeros(size(X));
De = Res;
SinT = Res;
Ta = Res;
Sa = Res;
Rho = Res;
M = Res;
Tb = Res;
Sb = Res;
GT = Res;
GS = Res;
Tf = Res;
E = Res;
P = Res;

% Fill up the arrays ready for plotting.

for i = 1:size(X,1)
   De(i) = interp1(IDdata(:,1),IDdata(:,2),X(i));
   SinT(i) = interp1(ISdata(:,1),ISdata(:,2),X(i));
   CosT = sqrt(1-SinT(i)^2);
   TS = interp1(TSdata(:,1),TSdata(:,2:3),...
                max([De(i)-D(i)*CosT/2; TSdata(1,1)]));
   Tad = TS(1,1);
   Sad = TS(1,2);
   Rho(i) = BS*(Sad-S(i))-BT*(Tad-T(i));
   [M(i),Tb(i),Sb(i),GT(i),GS(i),Tf(i)] = melt(D(i),U(i),T(i),S(i),...
                                               De(i));
%    [M(i),Tb(i),Sb(i),GT(i),GS(i),Tf(i)] = melt_mcphee(D(i),U(i),T(i),S(i),...
%                                                       De(i));
   TS = interp1(TSdata(:,1),TSdata(:,2:3),...
                max([De(i)-D(i)*CosT; TSdata(1,1)]));
   Ta(i) = TS(1,1);
   Sa(i) = TS(1,2);
   E(i) = E0*U(i)*SinT(i);
end

% Convert distance to km and entrainment and melt rate to m/yr.
% All other variables remain in SI.

X = X/1000;
E = E*3.16E7;
M = M*3.16E7;

%----------------------------------------------------------------------%
% Inform user of successful completion.

   fprintf('\nIntegration completed successfully.\n\n');

%----------------------------------------------------------------------%
% Tidy up workspace.

clear CosT IP Q1 Q2 Q3 Res Sad TS Tad Y Y0 Z i newDU newcf rep

%----------------------------------------------------------------------%