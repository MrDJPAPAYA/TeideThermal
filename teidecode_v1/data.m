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
Gp = 237; %[W/m2] Value at surface. Adjusted below

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
disp(strcat('Orbital period: ', string(round(T0)), ' seconds'));
Orb.n = 2*pi/T0; %Angular speed/mean motion [rad/s]
Orb.nu0 = 0; % Starting anomaly [rad] 
nu_f = @(t) Orb.nu0 + Orb.n*t; % True anomaly function [rad]
Gp = Gp*(Re/Orb.radi)^2;

% Rotation matrices
R.yaw = Rmatrix(3, -Att.yaw);
R.pitch = Rmatrix(2, -Att.pitch);
R.roll = Rmatrix(1, -Att.roll);
R.Att = R.yaw*R.pitch*R.roll; % Attitude rotation matrix
R.omega = Rmatrix(2, Orb.omega);
R.i = Rmatrix(3, Orb.i);
R.Orb = R.i*R.omega; % Orbital rotation matrix
R.nu_f = @(nu) Rmatrix(2, nu); % Mean anomaly rotation matrix function
R.pos_f = @(Rnu) Rnu*R.Orb; % Position rotation matrix (effect of orbit+anomaly)

% Alignment vectors
up = R.Att*[0; 0; 1]; % Earth pointing vector
us_f = @(Rpos) R.Att*Rpos*[0; 0; -1]; % Sun pointing vector

% Eclipse flag 
cosbeta_f = @(us) dot(up, us); % beta is the angle between the sun and earth pointing vectors
cosdelta = sqrt(1-Re^2/Orb.radi^2); % delta is the max value of beta while on eclipse
eclipse_flag = @(us) dot(up,us)>cosdelta;
 
%% Spacecraft data (Useless Right now)
A = 0.02; %Areas [m^2]
W = sqrt(A); %Width [m]

%% Import the nodes data
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
    SC(i).display = temp.NodesData(i, 11); % Boolean to see if it should be displayed on the plot
end

%% Import operation Conditions

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

%% Import solver config
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
%tf = 432000; %time to complete five days
tf = Mode(end).time-1;
ti = (t0:dt:tf).';

clear temp % Clear all temporary variables
 
%% Import conductance matrix.
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

