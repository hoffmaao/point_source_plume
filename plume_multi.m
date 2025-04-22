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


%----------------------------------------------------------------------%

% Clear variables that are not overwitten.


%----------------------------------------------------------------------%
% Load ice shelf base and ambient water column data.
% To do this call the script file DATALOAD.M.

temp_directory='/Users/andrew/Downloads/Plumemodel_point_source/tsdata/beardmore/'
temp_files = dir(fullfile(temp_directory,'*.dat'));
discharge_vec=[logspace(-5,1,20)]

M10km = zeros(length(temp_files),length(discharge_vec))
%MTstar10km = zeros(length(temp_files),length(discharge_vec))
%Mcombined = zeros(length(temp_files),length(discharge_vec))
for file=1:length(temp_files)

%----------------------------------------------------------------------%
% Define physical constants and make them global variables.
% To do this call the script file DEFINE.M.

    for dis_ind=1:length(discharge_vec)
        X = [];
        D = [];
        U = [];
        T = [];
        S = [];
        define
        tsfile= fullfile(temp_files(file).folder, temp_files(file).name);
        DUf = discharge_vec(dis_ind)
        disp(tsfile)

        %----------------------------------------------------------------------%
        % If no input data file exists, select file interactively.
        global IDdata ISdata TSdata GamT GamS



        infile = '/Users/andrew/Downloads/Plumemodel_original/indata/beardmore/bm_ideal_gl900.dat'
        INdata = load(fullfile(infile))
        IDdata = zeros(101,3);
        IDdata(:,1) = [0:1:100]'*(INdata(end,1)-INdata(1,1))/100;
        if INdata(end,1)-INdata(1,1) > 0
            IDdata(:,2) = interp1(INdata(:,1),INdata(:,2),IDdata(:,1),'pchip');
        else
            IDdata(:,2) = interp1(INdata(:,2)-INdata(1,2),INdata(:,2),...
                                  [0:1:100]'*(INdata(end,2)-INdata(1,2))/100,'linear');
        end
        IDdata(:,3) = (-110/920)*IDdata(:,2)+12;

        %----------------------------------------------------------------------%
        % Convert horizontal distance to distance along plume path.

        temp1 = IDdata(1,1);
        IDdata(1,1) = 0.0;
        for i = 2:size(IDdata,1)
           temp2 = IDdata(i,1);
           IDdata(i,1) = IDdata(i-1,1) + ...
                         sqrt((IDdata(i,2)-IDdata(i-1,2))^2 + ...
                              (temp2-temp1)^2);
           temp1 = temp2;
        end

        %----------------------------------------------------------------------%
        % Calculate basal slope and store it.

        ISdata = zeros(2*(size(IDdata,1)-1),2);
        for i = 1:size(IDdata,1)-1
           Slope = (IDdata(i+1,2)-IDdata(i,2))/(IDdata(i+1,1)-IDdata(i,1));
           ISdata(2*i-1:2*i,:) = [IDdata(i,1)+1.0 Slope; IDdata(i+1,1) Slope];
        end
        ISdata(1,1) = IDdata(1,1);
        ISdata(:,2) = min(ISdata(:,2),1);

        %TSdata = load(tsfile);
        TSdata = load(fullfile(tsfile))

        if TSdata(1,1) > IDdata (1,2)
           TSdata = [IDdata(1,2) TSdata(1,2:end); TSdata];
        end

        %----------------------------------------------------------------------%
        % Add new columns with water temperature and pressure freezing point.
        %TSdata = [TSdata(:,1) ...                       % Depth (z)
        %          TSdata(:,2)+FPa*TSdata(:,3)+FPb ...   % 
        %          TSdata(:,3) ...                       % Ambient salinity (Sa)
        %          TSdata(:,2) ...                       % Ambient temperature (Ta)
        %          FPa*TSdata(:,3)+FPb+FPc*TSdata(:,1)]; % Pressure freezing point at depth z (Tf)
                                                        % = Temperature at ice-ocean interface (Tb)
        TSdata = [TSdata(:,1) ...                       % Depth (z)
                  TSdata(:,2) ...                       % Ambient temperature (Ta)
                  TSdata(:,3:end) ...                   % Ambient salinity (Sa)
                  TSdata(:,2)-FPa*TSdata(:,3)-FPb ...   % T difference w.r.t to surface freezing point
                  FPa*TSdata(:,3)+FPb+FPc*TSdata(:,1)]; % Pressure freezing point at depth z (Tf)
                                                        % = Temperature at ice-ocean interface (Tb)


        clear INdata Slope

        %----------------------------------------------------------------------%
        % Prepare for integration.

        % Prompt user for initial conditions.
        % If values do not exist set defaults.

        %----------------------------------------------------------------------%
        % Check if user wishes to increase the ambient temperature.

        Tstar = trapz(TSdata(:,1),TSdata(:,4))/(TSdata(end,1)-TSdata(1,1))
        Temp_ambient=-str2num(append(tsfile(75),'.',tsfile(77)))        %if num2str(tsfile())
        
        if Temp_ambient~=Tstar
            Tinc = Tstar-Temp_ambient
            disp("temp isn't right")
            
        else
            Tinc = 0;
        end
        Tinc=0

        TSdata(:,2) = TSdata(:,2) + Tinc;
        TSdata(:,4) = TSdata(:,4) + Tinc;

        %----------------------------------------------------------------------%
        % Check if user wishes to generate repeated plumes until the end of the
        % domain is reached.  Default is NO.

        fprintf('\n');
        disp('Continue integration to end of domain using multiple plumes?');
        rep = 'Y'

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
           [D0,U0]=inlet(1000.0,900.0,DUf)

           %Rho = BS*(Sa-S0)-BT*(Ta-T0);
           %U0 = max([(Rho*G*SinT*DU/(E0*SinT+K))^0.333 Ueps]);
           %D0 = DU0/U0;

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
                                                       De(i),Tstar);
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
        
        [ d, ix ] = min(abs(X-10.));
        
        Mstar=1/3.16e-8*.25*(Tstar-Tf).^(4/3)*1e-6

        %----------------------------------------------------------------------%
        M10km(file,dis_ind)= trapz(X(1:ix),M(1:ix))/X(ix);
        MTstar10km(file,dis_ind)=trapz(X(1:ix),Mstar(1:ix))/X(ix);
        Mcombined(file,dis_ind)=trapz(X(1:ix),max(M(1:ix),Mstar(1:ix)))/X(ix);
        clearvars -except discharge_vec temp_files M10km MTstar10km Mcombined file disind IDdata IDdata ISdata TSdata GamT GamS
    end
end