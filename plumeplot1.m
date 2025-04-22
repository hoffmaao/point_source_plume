% PlumePlot.m
% -----------
%
% Script to plot basic parameters for plume.
% A. Jenkins, 29/08/07

% Ensure that the variables holding the figure handles exist.

if ~exist('fig1') | ~exist('fig2') | ~exist('fig3')
    fig1 = [];
    fig2 = [];
    fig3 = [];
end


% If none of the figures exist initialise the header, colour sequence and
% average accumulator.

if (~ishandle(fig1) & ~ishandle(fig2) & ~ishandle(fig3)) | ...
   (isempty(fig1) & isempty(fig2) & isempty(fig3))
    
    switch infile(1:5)
        case 'ideal'
            titlestring = 'Idealised Ice Shelf';
        case 'amery'
            titlestring = 'Amery Ice Shelf';
        case 'brunt'
            titlestring = 'Brunt Ice Shelf';
        case 'erebu'
            titlestring = 'Erebus Glacier Tongue';
        case 'ekstr'
            titlestring = 'Ekstrom Ice Shelf';
        case 'frisc'
            titlestring = 'Central Ronne Ice Shelf';
        case 'frisw'
            titlestring = 'Western Ronne Ice Shelf';
        case 'georg'
            titlestring = 'Northern George VI Ice Shelf';
        case 'pinei'
            titlestring = 'Pine Island Glacier';
        case 'riise'
            titlestring = 'Riiser-Larsen Ice Shelf';
        case 'rossj'
            titlestring = 'Ross Ice Shelf';
    end
    
    spectrum = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 0 1; 0 1 1; 1 1 0];
    colour = 0;
    
    meandata = [];

end

% Find x axis limits for these plots.

IDlim = [IDdata(1,1) IDdata(end,1)]/1000;

% Increment plot colour.

colour = mod(colour,7)+1;

% Calculate dervied values and means.

Flux = D.*U;
dRho = Rho0*Rho;
Melt = cumtrapz(X*1000,M/3.16E7);
Amb = cumtrapz(X*1000,E/3.16E7);
Meltf = Melt(2:end)./(Melt(2:end)+Amb(2:end));
Meltf = [0; Meltf];
MeltFlux = Flux.*Meltf*3.16E7/1e6;
Force1 = Rho0*Rho.*G.*SinT.*D;
Force2 = -Rho0*K.*U.*U;
Heat1 = Rho0*C*(Ta-T).*E/3.16E7;
Heat2 = Rho0*C*(Tb-T).*(GT+M/3.16E7);
Buoyf = Rho0*G*Flux.*Rho;

meandata = [meandata; ...
            De(1) trapz(X,SinT)/X(end) trapz(X,Ta-FPa*Sa-FPb)/X(end) ...
            trapz(X,M)/X(end) trapz(X,D)/X(end) trapz(X,U)/X(end) ...
            trapz(X,T-Tb)/X(end) trapz(X,Ta-T)/X(end) ...
            trapz(X,S-Sb)/X(end) trapz(X,Sa-S)/X(end)];

% Plot the ice shelf thickness profile.

if ishandle(fig1)
    
    figure(fig1);

else

    fig1 = figure('Units','centimeters','Position',[1 2 20 20],...
                  'Name',[titlestring ' 1'],'NumberTitle','off');
              
    axea = axes('Position',[0.1 0.56 0.4 0.4],'Box','on','XGrid','on',...
                'YGrid','on','NextPlot','add','XLim',IDlim,...
                'YLim',[min(IDdata(:,2))*1.1 max(IDdata(:,3))*1.1]);
    xlabel('Distance (km)');
    ylabel('Depth (m)');
    
    axeb = axes('Position',[0.1 0.375 0.4 0.125],'Box','on','XGrid','on',...
                'YGrid','on','NextPlot','add','XLim',IDlim,...
                'YLim',[min([0 min(ISdata(:,2))*1.1]) max(ISdata(:,2))*1.1]);
    xlabel('Distance (km)');
    ylabel('Basal slope');
    
    axec = axes('Position',[0.6 0.56 0.15 0.4],'Box','on','XGrid','on',...
                'YGrid','on','NextPlot','add',...
                'XLim',[min([0 min(TSdata(:,5))*1.1]) ...
                        max([0 max(TSdata(:,2))*1.1])],...
                'YLim',[min(IDdata(:,2))*1.1 max(IDdata(:,3))*1.1]);
    xlabel('Temperature (\circC)');
    ylabel('Depth (m)');
    
    axed = axes('Position',[0.84 0.56 0.15 0.4],'Box','on','XGrid','on',...
                'YGrid','on','NextPlot','add',...
                'XLim',[min(TSdata(:,3))*0.975 max(TSdata(:,3))*1.025],...
                'YLim',[min(IDdata(:,2))*1.1 max(IDdata(:,3))*1.1]);
    xlabel('Salinity');
    ylabel('Depth (m)');

    axee = axes('Position',[0.1 0.06 0.4 0.25],'Box','on','NextPlot','add',...
                'XGrid','on','YGrid','on');
    xlabel('Ambient water temperature w.r.t. surface freezing point (\circC)');
    ylabel('Mean basal melt rate (m/yr)');
    
    axef = axes('Position',[0.6 0.32 0.15 0.175],'Box','on','NextPlot','add',...
                'XGrid','on','YGrid','on');
    xlabel('T_{a} w.r.t. T_{f} (\circC)');
    ylabel('Mean(D) (m)');
    
    axeg = axes('Position',[0.84 0.32 0.15 0.175],'Box','on','NextPlot','add',...
                'XGrid','on','YGrid','on');
    xlabel('T_{a} w.r.t. T_{f} (\circC)');
    ylabel('Mean(U) (m/s)');
    
    axeh = axes('Position',[0.6 0.06 0.15 0.175],'Box','on','NextPlot','add',...
                'XGrid','on','YGrid','on');
    xlabel('T_{a} w.r.t. T_{f} (\circC)');
    ylabel('Mean(T-T_{b} & T_{a}-T) (\circC)');
    
    axei = axes('Position',[0.84 0.06 0.15 0.175],'Box','on','NextPlot','add',...
                'XGrid','on','YGrid','on');
    xlabel('T_{a} w.r.t. T_{f} (\circC)');
    ylabel('Mean(S-S_{b} & S_{a}-S)');
    
end

axes(axea);
xlim = get(axea,'XLim');
xlim = [min([xlim(1) IDlim(1)]) max([xlim(2) IDlim(2)])];
set(axea,'XLim',xlim);
ylim = get(axea,'YLim');
ylim = [min([ylim(1) min(IDdata(:,2))*1.1/1000]) ...
        max([ylim(2) max(IDdata(:,3))*1.1/1000])];
set(axea,'YLim',ylim);
plot(IDdata(:,1)/1000,IDdata(:,2),'-','Color',spectrum(colour,:));
plot(IDdata(:,1)/1000,IDdata(:,3),'-','Color',spectrum(colour,:));
      
axes(axeb);
set(axeb,'XLim',xlim);
ylim = get(axeb,'YLim');
ylim = [min([ylim(1) min(ISdata(:,2))*1.1]) max([ylim(2) max(ISdata(:,2))*1.1])];
set(axeb,'YLim',ylim);
plot(ISdata(:,1)/1000,ISdata(:,2),'-','Color',spectrum(colour,:));

axes(axec);
tlim = get(axec,'XLim');
tlim = [min([tlim(1) min(TSdata(:,5))*1.1]) max([tlim(2) max(TSdata(:,2))*1.1])];
set(axec,'XLim',tlim);
ylim = get(axec,'YLim');
ylim = [min([ylim(1) min(IDdata(:,2))*1.1/1000]) ...
        max([ylim(2) max(IDdata(:,3))*1.1/1000])];
set(axec,'YLim',ylim);
plot(TSdata(:,2),TSdata(:,1),'-','Color',spectrum(colour,:));
plot(TSdata(:,5),TSdata(:,1),'--','Color',spectrum(colour,:));

axes(axed);
slim = get(axed,'XLim');
slim = [min([slim(1) min(TSdata(:,3))*0.975]) max([slim(2) max(TSdata(:,3))*1.025])];
set(axed,'XLim',slim);
ylim = get(axed,'YLim');
ylim = [min([ylim(1) min(IDdata(:,2))*1.1/1000]) ...
        max([ylim(2) max(IDdata(:,3))*1.1/1000])];
set(axed,'YLim',ylim);
plot(TSdata(:,3),TSdata(:,1),'-','Color',spectrum(colour,:));

axes(axee);
plot(meandata(end,3),meandata(end,4),'*','Color',spectrum(colour,:));

axes(axef);
plot(meandata(end,3),meandata(end,5),'*','Color',spectrum(colour,:));

axes(axeg);
plot(meandata(end,3),meandata(end,6),'*','Color',spectrum(colour,:));

axes(axeh);
plot(meandata(end,3),meandata(end,7),'*','Color',spectrum(colour,:));
plot(meandata(end,3),meandata(end,8),'*','Color',spectrum(colour,:));

axes(axei);
plot(meandata(end,3),meandata(end,9),'*','Color',spectrum(colour,:));
plot(meandata(end,3),meandata(end,10),'*','Color',spectrum(colour,:));


if ishandle(fig2)
    
    figure(fig2);

else
    
    fig2 = figure('Units','centimeters','Position',[2 3 27 18],...
                  'Name',[titlestring ' 2'],'NumberTitle','off');

    axe1 = axes('Position',[0.075 0.55 0.25 0.4],'Box','on','XGrid','on',...
                'YGrid','on','NextPlot','add','XLim',IDlim,...
                'YLim',[min([0 min(D)*1.1]) max(D)*1.1]);
    xlabel('Distance (km)');
    ylabel('Plume thickness (m)');

    axe2 = axes('Position',[0.4 0.55 0.25 0.4],'Box','on','XGrid','on',...
                'YGrid','on','NextPlot','add','XLim',IDlim,...
                'YLim',[min([0 min(U)*1.1]) max(U)*1.1]);
    xlabel('Distance (km)');
    ylabel('Plume velocity (m/s)');

    axe3 = axes('Position',[0.725 0.55 0.25 0.4],'Box','on','XGrid','on',...
                'YGrid','on','NextPlot','add','XLim',IDlim,...
                'YLim',[min([0 min(E)*1.1]) max(E)*1.1]);
    xlabel('Distance (km)');
    ylabel('Entrainment rate (m/yr)');

    axe4 = axes('Position',[0.075 0.075 0.25 0.4],'Box','on','XGrid','on',...
                'YGrid','on','NextPlot','add','XLim',IDlim,...
                'YLim',[33 max([35 max(Sb)*1.01])]);
    xlabel('Distance (km)');
    ylabel('Plume/boundary salinity');

    axe5 = axes('Position',[0.4 0.075 0.25 0.4],'Box','on','XGrid','on',...
                'YGrid','on','NextPlot','add','XLim',IDlim,...
                'YLim',[min([0 min(Tb)*1.1]) T(end)*(1+sign(T(end))/10)]);
    xlabel('Distance (km)');
    ylabel('Plume/boundary temperature (\circC)');

    axe6 = axes('Position',[0.725 0.075 0.25 0.4],'Box','on','XGrid','on',...
                'YGrid','on','NextPlot','add','XLim',IDlim,...
                'YLim',[min([0 min(M)*1.1]) max(M)*1.1]);
    xlabel('Distance (km)');
    ylabel('Ice shelf basal melt rate (m/yr)');
    if exist([infile '_melt.mdat'])==2;
        obsmelt = load([infile '_melt.mdat']);
        plot(obsmelt(:,1),obsmelt(:,2),'k+');
    end

    
end
    
axes(axe1);
set(axe1,'XLim',xlim);
ylim = get(axe1,'YLim');
ylim = [min([ylim(1) min(D)*1.1]) max([ylim(2) max(D)*1.1])];
set(axe1,'YLim',ylim);
plot(X,D,'-','Color',spectrum(colour,:));

axes(axe2);
set(axe2,'XLim',xlim);
ylim = get(axe2,'YLim');
ylim = [min([ylim(1) min(U)*1.1]) max([ylim(2) max(U)*1.1])];
set(axe2,'YLim',ylim);
plot(X,U,'-','Color',spectrum(colour,:));

axes(axe3);
set(axe3,'XLim',xlim);
ylim = get(axe3,'YLim');
ylim = [min([ylim(1) min(E)*1.1]) max([ylim(2) max(E)*1.1])];
set(axe3,'YLim',ylim);
plot(X,E,'-','Color',spectrum(colour,:));

axes(axe4);
set(axe4,'XLim',xlim);
ylim = get(axe4,'YLim');
ylim = [33 max([ylim(2) max(Sb)*1.01])];
set(axe4,'YLim',ylim);
plot(X,S,'-','Color',spectrum(colour,:));
plot(X,Sb,'--','Color',spectrum(colour,:));

axes(axe5);
set(axe5,'XLim',xlim);
ylim = get(axe5,'YLim');
ylim = [min([ylim(1) min(Tf)*1.1]) max([ylim(2) T(end)*(1+sign(T(end))/10)])];
set(axe5,'YLim',ylim);
plot(X,T,'-','Color',spectrum(colour,:));
plot(X,Tb,'--','Color',spectrum(colour,:));

axes(axe6);
set(axe6,'XLim',xlim);
ylim = get(axe6,'YLim');
ylim = [min([ylim(1) min(M)*1.1]) max([ylim(2) max(M)*1.1])];
set(axe6,'YLim',ylim);
plot(X,M,'-','Color',spectrum(colour,:));


if ishandle(fig3)
    
    figure(fig3);
    
else
        
    fig3 = figure('Units','centimeters','Position',[3 2 27 18],...
                  'Name',[titlestring ' 3'],'NumberTitle','off');

    axe7 = axes('Position',[0.075 0.55 0.25 0.4],'Box','on','XGrid','on',...
                'YGrid','on','NextPlot','add','XLim',IDlim,...
                'YLim',[min([0 min(Flux)*1.1]) max(Flux)*1.1]);
    xlabel('Distance (km)');
    ylabel('Plume flux per unit width (m^{2}/s)');

    axe8 = axes('Position',[0.4 0.55 0.25 0.4],'Box','on','XGrid','on',...
                'YGrid','on','NextPlot','add','XLim',IDlim,...
                'YLim',[0 0.5]);
    xlabel('Distance (km)');
    ylabel('Density difference (kg/m^{3})');

    axe9 = axes('Position',[0.725 0.55 0.25 0.4],'Box','on','XGrid','on',...
                'YGrid','on','NextPlot','add','XLim',IDlim,...
                'YLim',[min([0 min(Meltf)*1.1]) max(Meltf)*1.1]);
    xlabel('Distance (km)');
    ylabel('Meltwater fraction');

%     axe9 = axes('Position',[0.725 0.55 0.25 0.4],'Box','on','XGrid','on',...
%                 'YGrid','on','NextPlot','add','XLim',IDlim,...
%                 'YLim',[min([0 min(Buoyf)*1.1]) max(Buoyf)*1.1]);
%     xlabel('Distance (km)');
%     ylabel('Meltwater fraction');

    axe10 = axes('Position',[0.075 0.075 0.25 0.4],'Box','on','XGrid','on',...
                 'YGrid','on','NextPlot','add','XLim',IDlim,...
                 'YLim',[min(Force2)*1.1 max(Force1)*1.1]);
    xlabel('Distance (km)');
    ylabel('Force per unit area (N/m^{2})');

    axe11 = axes('Position',[0.4 0.075 0.25 0.4],'Box','on','XGrid','on',...
                 'YGrid','on','NextPlot','add','XLim',IDlim,...
                 'YLim',[min(Heat2)*1.1 max(Heat1)*1.1]);
    xlabel('Distance (km)');
    ylabel('Heat flux (W/m^{2})');

    axe12 = axes('Position',[0.725 0.075 0.25 0.4],'Box','on','XGrid','on',...
                 'YGrid','on','NextPlot','add','XLim',IDlim,...
                 'YLim',[min([0 min(MeltFlux)*1.1]) max(MeltFlux)*1.1]);
    xlabel('Distance (km)');
    ylabel('Meltwater flux per unit width (km^{2}/yr)');

end

axes(axe7);
set(axe7,'XLim',xlim);
ylim = get(axe7,'YLim');
ylim = [min([ylim(1) min(Flux)*1.1]) max([ylim(2) max(Flux)*1.1])];
set(axe7,'YLim',ylim);
plot(X,Flux,'-','Color',spectrum(colour,:));

axes(axe8);
set(axe8,'XLim',xlim);
plot(X,dRho,'-','Color',spectrum(colour,:));

axes(axe9);
set(axe9,'XLim',xlim);
ylim = get(axe9,'YLim');
ylim = [min([ylim(1) min(Meltf)*1.1]) max([ylim(2) max(Meltf)*1.1])];
set(axe9,'YLim',ylim);
plot(X,Meltf,'-','Color',spectrum(colour,:));

% axes(axe9);
% set(axe9,'XLim',xlim);
% ylim = get(axe9,'YLim');
% ylim = [min([ylim(1) min(Buoyf)*1.1]) max([ylim(2) max(Buoyf)*1.1])];
% set(axe9,'YLim',ylim);
% plot(X,Buoyf,'-','Color',spectrum(colour,:));

axes(axe10);
set(axe10,'XLim',xlim);
ylim = get(axe10,'YLim');
ylim = [min([ylim(1) min(Force2)*1.1]) max([ylim(2) max(Force1)*1.1])];
set(axe10,'YLim',ylim);
plot(X,Force1+Force2,'-','Color',spectrum(colour,:));
plot(X,Force1,':','Color',spectrum(colour,:));
plot(X,Force2,'--','Color',spectrum(colour,:));

axes(axe11);
set(axe11,'XLim',xlim);
ylim = get(axe11,'YLim');
ylim = [min([ylim(1) min(Heat2)*1.1]) max([ylim(2) max(Heat1)*1.1])];
set(axe11,'YLim',ylim);
plot(X,Heat1+Heat2,'-','Color',spectrum(colour,:));
plot(X,Heat1,':','Color',spectrum(colour,:));
plot(X,Heat2,'--','Color',spectrum(colour,:));

axes(axe12);
set(axe12,'XLim',xlim);
ylim = get(axe12,'YLim');
ylim = [min([ylim(1) min(MeltFlux)*1.1]) max([ylim(2) max(MeltFlux)*1.1])];
set(axe12,'YLim',ylim);
plot(X,MeltFlux,'-','Color',spectrum(colour,:));