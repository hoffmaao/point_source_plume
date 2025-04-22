%----------------------------------------------------------------------%
% DEFINE.M                                               [Version 1.00]
% ********
%----------------------------------------------------------------------%
% Script file to define physical constants used in ISW plume model.
%
% Version 1.00
% Created 23/10/97 by A Jenkins.
%----------------------------------------------------------------------%

% Define physical constants and make them global variables.

global E0 K Ueps DU BS BT Rho0 G L Ci C FPa FPb FPc GamT GamS GamTS

% Entrainment and drag coefficients.
% (These depend on the amount of turbulent mixing).

E0 = 0.036;
K = 2.5E-3;

% Minimum velocity and flux, used to avoid a divide by zero.

Ueps = 1e-6;
DU = 1e-6;

% Themal expansion and haline contraction coefficients.
% (These depend on the water properties).

BS = 7.86E-4;
BT = 3.87E-5;
Rho0 = 1030;

% Acceleration due to gravity.

G = 9.81;

% Thermal properties of ice and seawater.
% (These are relatively constant).

L = 3.35E5;
Ci = 2009.0;
C = 3974.0;
 
% Coefficients in the freezing point of seawater.
% (These depend on the water properties).

FPa = -0.0573;
FPb = 0.0832;
FPc = 7.61E-4;

% Thermal and haline Stanton numbers.
% (These depend on the amount of turbulent mixing).

GamT = 0.0011;
GamS = 3.1E-5;
GamTS = 5.9E-4;

%----------------------------------------------------------------------%
% Choose ice properties.
% These are also global variables.

global Ti Si
 
% Temperature and salinity of ice.
% The temperature determines the amount of heat conducted into the ice
% shelf.  It only affects melting, if freezing occurs the conduction
% term is assumed to be zero.
% The salinity determines the amount of salt trapped in the ice on
% freezing.  For melting the salinity is always zero.

Ti = -15.0;
Si = 0.0;

%----------------------------------------------------------------------%