%----------------------------------------------------------------------%
% DATALOAD.M                                             [Version 2.00]
% **********
%----------------------------------------------------------------------%
% Script to load data needed by ISW plume model from specified files.
%
% User must provide two files of input data:
%      1) Shape of ice shelf base.
%           1st column: distance from grounding line in m.
%           [For a vertical ice face set all rows in this column to 0]
%           2nd column: elevation of ice shelf base with respect to
%                       sea level (ie negative) in m.
%           1st row:  grouding line.
%           Last row: ice front.
%      2) Properties of ambient water column.
%           1st column: elevation with respect to sea level
%                       (ie negative) in m.
%           2nd column: water temperature in degrees C.
%           3rd column: water salinity in psu.
%           1st row:  bottom of water column.
%           [For a uniform water column only this row is required]
%           Last row: top of water column.
% Script creates five variables:
%      infile : name of file 1.
%      tsfile : name of file 2.
%      IDdata : elevation of ice shelf base.
%      ISdata : slope of ice shelf base.
%      TSdata : ambient temperature and salinity.
%
% Version 2.00
% Modified 28/08/07 by A Jenkins.
% Added interactive selection of data files to replace manual entry of
% file names.
% Temperature in file is now given w.r.t. surface freezing point.
%
% Version 1.00
% Created 23/10/97 by A Jenkins.
%----------------------------------------------------------------------%

% Declare data variables as global.

global IDdata ISdata TSdata

%----------------------------------------------------------------------%
% If no input data file exists, select file interactively.

if exist('infile')
    fprintf('\n');
    disp(['Use ice shelf basal elevation data from "' infile ...
          '.dat" (Y/N)?']);
    useoldindata = input('   [Default=Y]:   ','s');
    if isempty(useoldindata)
        useoldindata = 'Y';
    end
else
    useoldindata = 'N';
end

%----------------------------------------------------------------------%
% Read basal elevation data if a new file was selected.
    
if useoldindata == 'N' | useoldindata == 'n'
    
    [infile,inpath] = uigetfile('*.dat',...
                                ['Choose basal elevation data' ...
                                 'file from indata directory']);
    infile = infile(1:end-4);
    eval(['load ', inpath, infile,'.dat;']);
    eval(['INdata = ', infile,';']);
                            
%----------------------------------------------------------------------%
% Eliminate sharp changes in basal slope by use of piecewise cubic
% hermite interpolation onto 101 equally spaced points.

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
    
end

%----------------------------------------------------------------------%
% If no T/S data file exists, select file interactively.

if exist('tsfile')
    fprintf('\n');
    disp(['Use ice temperature/salinty data from "' tsfile ...
          '.dat" (Y/N)?']);
    useoldtsdata = input('   [Default=Y]:   ','s');
    if isempty(useoldtsdata)
        useoldtsdata = 'Y';
    end
else
    useoldtsdata = 'N';
end

%----------------------------------------------------------------------%
% Read ambient water properties if a new file was selected.
    
if useoldtsdata == 'N' | useoldtsdata == 'n'
    
    [tsfile,tspath] = uigetfile('*.dat',...
                                ['Choose temperature/salinity '...
                                 'data file from tsdata directory']);
    tsfile = tsfile(1:end-4);
    eval(['load ', tspath, tsfile,'.dat;']);
    eval(['TSdata = ', tsfile,';']);
                            
%----------------------------------------------------------------------%
% Ensure that the ambient water column is deep enough.
% If it is not, extend it to the depth of the grounding line.

    if TSdata(1,1) > IDdata (1,2)
       TSdata = [IDdata(1,2) TSdata(1,2:end); TSdata];
    end

%----------------------------------------------------------------------%
% Add new columns with water temperature and pressure freezing point.

    if isempty(strfind(tspath,'real'))

        TSdata = [TSdata(:,1) TSdata(:,2)+FPa*TSdata(:,3)+FPb ...
                  TSdata(:,3) TSdata(:,2) ...
                  FPa*TSdata(:,3)+FPb+FPc*TSdata(:,1)];
              
    else
        
        TSdata = [TSdata(:,1) TSdata(:,2) ...
                  TSdata(:,3:end) TSdata(:,2)-FPa*TSdata(:,3)-FPb ...
                  FPa*TSdata(:,3)+FPb+FPc*TSdata(:,1)];
              
    end
          
end

%----------------------------------------------------------------------%

% Tidy up workspace.

clear INdata Slope

%----------------------------------------------------------------------%