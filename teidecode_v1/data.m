clear;close all;
set(groot,'DefaultAxesFontSize',11);
set(groot,'DefaultTextInterpreter','latex');

%% Thermal Subsytem

%% Constants
Tc = 32+273.15; %[K]
Re = 6371e3; % Earth radius [m]
g = 9.81; %[m/s2]
a = 0.4; %Earth albedo
sigma = 5.67E-8; %[W/m2/K4]
Gs0 = 1371; %[W/m2]
Gp = 237; %[W/m2] Value at surface. Could adjust this for distance on cold critical cases.

%% Orbit and orientation
temp.orbit_params = readtable("Thermal_Data.xlsx", 'Sheet',"Orbital Parameters", VariableNamingRule="preserve");
temp.orbit_params = table2array(temp.orbit_params);

% Attitude parameters [rad]
Att.yaw = deg2rad(temp.orbit_params(1, 1)); 
Att.pitch = deg2rad(temp.orbit_params(1, 2));
Att.roll = deg2rad(temp.orbit_params(1, 3));

% Orbital parameters
Orb.omega = deg2rad(temp.orbit_params(1, 4)); % Longitude of Ascending node[rad] 
Orb.i = deg2rad(temp.orbit_params(1, 5)); % Inclination [rad]
Orb.radi = Re + 1000*temp.orbit_params(1, 6); % Radius[m]
T0 = 2*pi/Re*sqrt(Orb.radi^3/g); % Orbital period [s] 
Orb.n = 2*pi/T0; %Angular speed/mean motion [rad/s]
Orb.nu0 = 0; % Starting anomaly [rad] 
nu_f = @(t) Orb.nu0 + Orb.n*t; % True anomaly function [rad]

% Rotation matrices
R.yaw = Rmatrix(3, -Att.yaw);
R.pitch = Rmatrix(2, -Att.pitch);
R.roll = Rmatrix(1, -Att.roll);
R.Att = R.yaw*R.pitch*R.roll; % Attitude rotation matrix
R.omega = Rmatrix(2, Orb.omega);
R.i = Rmatrix(3, Orb.i);
R.Orb = R.omega*R.i; % Orbital rotation matrix
R.nu_f = @(nu) Rmatrix(2, nu); % Mean anomaly rotation matrix function
R.pos_f = @(Rnu) R.Orb*Rnu; % Position rotation matrix (effect of orbit+anomaly)

% Alignment vectors
up = R.Att*[0; 0; 1]; % Earth pointing vector
us_f = @(Rpos) Rpos*R.Att*[0; 0; -1]; % Sun pointing vector

% Eclipse flag 
cosbeta_f = @(us) dot(up, us); % beta is the angle between the sun and earth pointing vectors
cosdelta = sqrt(1-Re^2/Orb.radi^2); % delta is the max value of beta while on eclipse
eclipse_flag = @(us) dot(up,us)>cosdelta;


%% Orbit
%orb_case = input('Input orbit case to run: ');
%switch orb_case
%    case 1
%        theta_SC = deg2rad(0); %[rad]
%        phi_SC = 0; %[rad]
%    case 2
%        theta_SC = deg2rad(45); %[rad]
%        phi_SC = 0; %[rad]
%    otherwise
%        orb_case = 1;
%        theta_SC = deg2rad(0); %[rad]
%        phi_SC = 0; %[rad]
%end

%All of these are orbital parameters
%ha = 400e3; %[m]
%r = ha + Re; %[m]
%T0 = 2*pi/Re*sqrt(r^3/g); %[s]
%omega0 = 2*pi/T0; %[rad/s]
%gamma0 = 0; %[rad]
%gamma_f = @(t) gamma0 + 2*pi/T0.*t; %[rad]
%beta_f = @(gamma) abs(wrapToPi(gamma)); %[rad]
%up_f = @(theta,phi) [-sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];
%us_f = @(gamma,up) Rmatrix(2,theta_SC).'*[-sin(gamma); 0; cos(gamma)];
%delta = acos(Re/r);
%eclipse_flag = @(gamma) (pi-abs(wrapToPi(gamma)))<delta;



 %%  Nodes On/Off condition
 %missing obc and power subsystems
% Nota: A las 16:00 de cada día el satélite entra en active mode durante 20 minutos para self test y copiar las memorias 
 %day 1
 
 %(85800:86400) reorientation mode 23:50
 %(86400:86410) payload mode, La siguiente ventana es a las 2:00 y seguimos el mismo proceso?
 %SC(25200:25500) radio on on day 1 at 7:00
 
 %day 2
 %90000:90300 radio mode day 2 1:00 p.m.
 %115200:115500 radio mode day 2 8:00 p.m.
 
 %day 3 idle
 %day 4 idle
 
 %day 5
 %435600:435900 radio mode day 5 1:00 p.m.
 %460800:461100 radio mode day 5 8:00 p.m.
 
 %day 6 same as day 1
 
 

   
%% Spacecraft data
A = 0.02; %Areas [m^2]
W = sqrt(A); %Width [m]

%structure nodes 
% Import the data
temp.NodesData = readtable("Thermal_Data.xlsx", 'Sheet',"Nodes thermal properties", VariableNamingRule="preserve");

% Convert to output type
temp.NodesNames = temp.NodesData{ :, 1};
temp.NodesData = temp.NodesData(1:end,2:end); %Crop text header column
temp.NodesData = table2array(temp.NodesData);

% Assign into struct element and save number

N = size(temp.NodesData, 1);

SC = struct('name', [], 'm', [], 'Cp', [], 'k', [], 'A', [], 'e', [], 'a', [], 'n', [], 'radiates', []);

for i = 1:N
    SC(i).name = temp.NodesNames{i, 1}; % Name
    SC(i).m = temp.NodesData(i, 1); % Mass
    SC(i).Cp = temp.NodesData(i, 2); % Specific Heat
    SC(i).k = temp.NodesData(i, 3); % Conductivity
    SC(i).A = temp.NodesData(i, 4); % Area
    SC(i).e = temp.NodesData(i, 5); % Emissivity
    SC(i).a = temp.NodesData(i, 6); % Absorptivity
    SC(i).n = [temp.NodesData(i, 7); temp.NodesData(i, 8); temp.NodesData(i, 9)]; % Normal vector to surface
    SC(i).radiates = temp.NodesData(i, 10); % Boolean to see if it radiates
end

%% Operation Conditions

% Read Operation modes settings
    temp.ModesOp = readtable("Thermal_Data.xlsx", 'Sheet',"Operation");
    temp.ModesTimes = temp.ModesOp(1, 3:end);
    temp.ModesTimes = table2array(temp.ModesTimes);
    temp.ModesHeats = temp.ModesOp(6:end, 3:end);
    temp.ModesHeats = table2array(temp.ModesHeats);
    temp.ModesNames = temp.ModesOp.Properties.VariableNames; 
    temp.ModesNames = temp.ModesNames(1, 3:end);
Mode = struct('name', [], 'time', [], 'heats', []);

if size(temp.ModesHeats, 1) ~= N
error("Node properties and number of nodes per operation mode have different sizes. Check data")
end 

for i = 1:size(temp.ModesNames, 2)
    Mode(i).name = temp.ModesNames{1, i};
    Mode(i).time = temp.ModesTimes(1, i);
    Mode(i).HeatGen = temp.ModesHeats(1:end, i);
end



%% Solver config
%these parameters will define the duration fo the simulation, if tf=10*T0
%the simulation will run for 10 orbits
% Read Startup conditions
temp.SolverConfig = readtable("Thermal_Data.xlsx", 'Sheet',"Startup Parameters", VariableNamingRule="preserve");
temp.SolverConfig = table2array(temp.SolverConfig);
config = struct('SolRad', [], 'EnvRad', []);
config.SolRad =  temp.SolverConfig(1, 1);
config.EnvRad =  temp.SolverConfig(1, 2);
EnvT =  temp.SolverConfig(1, 3);
InitialT =  temp.SolverConfig(1, 4);
dt = temp.SolverConfig(1, 5); %[s]
t0 = 0;
%tf = 10*T0; %T0*77.85; %5 days, with 15,57 orbits per day aproximately
%tf = 432000; %time to complete five days
tf = Mode(end).time-1;
ti = (t0:dt:tf).';

clear temp % Clear all temporary variables

 %leds power generation
 SCledqgen = zeros(1,length(ti)); 
 %SCledqgen (1,44300:44350)=50; %W
 %SCledqgen (1,86400:86410)=50; %W dato cotejado con juan (aproximacion)
 SCledqgen (1,1:length(ti))=5.121; %W for stationary tests

 %radio mode, heat applied to trasceiver and antenas?
 SCradiomode = zeros(1,length(ti)); 
 SCradiomode (1,25200:25500) = 5; %W at day 1 supuestos
 SCradiomode (1,115200:115500) = 5; %W at day 2
 SCradiomode (1,435600:435900) = 5; %W at day 5
 SCradiomode (1,460800:461100) = 5; %W at day 5

 %reorientation mode 
  SCreomode = zeros(1,length(ti));
  SCreomode (1,85800:86400) = 5; %W at day 1 (i should ask someone about this shit)

%H45 = sigma/(SC(4).R + 1/A + SC(5).R);

% K(1,2) = SC(1).k*A*2*SC(1).L;
% K(2,3) = SC(3).k*A*2/SC(3).L;
% K(3,4) = 1/(SC(3).L/2/SC(3).k/A + SC(4).L/2/SC(4).k/A);
% K(1,6:9) = 4e-2;
% K(2,6:9) = 0.2;
% K(4,6:9) = 4e-4;
% K(6,[7,9]) = 0.4;
% K(7,8) = 0.4;
% K(8,9) = 0.4;
% K = K.'+K;

% Conductances

 %K = ThermalDataS2(ThermalDataS1);
 
 
%% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: C:\Users\erjav\Desktop\TeideThermal\TeideThermal\teidecode_v1\Thermal_Data.xlsx
%    Worksheet: Nodes properties
%
% Auto-generated by MATLAB on 29-Oct-2023 13:45:26

%% Setup the Import Options and import the data (struct) 
opts = spreadsheetImportOptions("NumVariables", 11);

% Specify sheet and range
opts.Sheet = "Nodes properties";
opts.DataRange = "B40:L69";

% Specify column names and types
opts.VariableNames = ["VarName2", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11", "VarName12"];
opts.VariableTypes = ["categorical", "string", "double", "categorical", "categorical", "double", "double", "string", "double", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, ["VarName3", "VarName9"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["VarName2", "VarName3", "VarName5", "VarName6", "VarName9"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, ["VarName4", "VarName7", "VarName8", "VarName10", "VarName11", "VarName12"], "FillValue", 0);

%SC = ThermalData; % THIS DOES NOT WO,RK, COMO IMPORTO LA TABLA A UN Struct
%% Clear temporary variables
clear opts
 
%% Import data from spreadsheet (MATRIZ K)
% Script for importing data from the following spreadsheet:
%
%    Workbook: C:\Users\erjav\Desktop\TeideThermal\TeideThermal\teidecode_v1\Thermal_Data.xlsx
%    Worksheet: Conductances_between_nodes
%
% Auto-generated by MATLAB on 07-Apr-2023 12:17:07


%% Setup the Import Options and import the data
%opts = spreadsheetImportOptions("NumVariables", 29);

% Specify sheet and range
%opts.Sheet = "Conductances_between_nodes";
%opts.DataRange = "B5:AD33";

% Specify column names and types
%opts.VariableNames = ["node1", "node2", "node3", "node4", "node5", "node6", "node7", "node8", "node9", "node10", "node11", "node12", "node13", "node14", "node15", "node16", "node17", "node18", "node19", "node20", "node21", "node22", "node23", "node24", "node25", "node26", "node27", "node28", "node29"];
%opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Import the data
ThermalDataS1 = readtable("Thermal_Data.xlsx", 'Sheet',"Conductances_between_nodes", VariableNamingRule="preserve");

% Convert to output type
ThermalDataS1 = ThermalDataS1(1:end,2:end); %Crop text headers
ThermalDataS1 = table2array(ThermalDataS1);

% Asign and check size
clear opts
K = ThermalDataS1;
if size(K, 1) ~= N
error("Node properties and conductance matrix have different sizes. Check data")
end

%{ 
 %K = zeros(N); %[W/K]
 K(1,3:6)=0.37;
 K(2,3:6)=0.4;
 K(1,26)=0.8;
 K(2,27)=0.8;
 K(3,[7,12,13,20,21,23])=0.2;
 K(4,[7,12,13,20,21,23])=0.2;
 K(5,[12,13,20,21,23,22])=0.2;
 K(6,[12,13,20,21,23,22])=0.2;
 K(7,8:11)= 0.9; %danger
 K(12,8:11)= 0.4; %danger
 K(13,14:19)=0.8;
 K(20,[21,23])=0.2;
 K(1,26)=0.8;
 K(2,27)=0.8;
 K(24,[1,2,4,5])=0.8;
 K(25,[1,2,3,6])=0.8;
 K(28,[1, 2 ,5, 6])=0.3;
 K(28,2)=0.8;
 K(29,[3, 4, 5, 6, 8 ,9 ,10, 11])=0.8;
 K = K.'+K;
%}
%% Auxiliar functions
function [R] = Rmatrix(raxis,rangle)
%Rotation matrix:
%   raxis: rotation axis
%       x -> 1
%       y -> 2
%       z -> 3
%   rangle: rotation angle (rad)

    R = zeros(3,3);
    S = sin(rangle);
    C = cos(rangle);
    A = [C,-S;S,C];
    if raxis == 2
        A = A.';
    end
    v = zeros(1,2);
    
    k = 0;
    for i = 1:3
        if i == raxis
            R(i,i) = 1;
        else
            k = k + 1;
            v(k) = i;
        end
    end
    
    n = 0;
    for i = v
        n = n + 1;
        l = 0;
        for j = v
            l = l + 1;
            R(i,j) = A(n,l);
        end
    end
end

