clc;clear;close all;
set(groot,'DefaultAxesFontSize',11);
set(groot,'DefaultTextInterpreter','latex');

%% Thermal Subsytem

%% Constants
Tc = 32+273.15; %[K]
Re = 6371e3; %[m]
g = 9.81; %[m/s2]
a = 0.4; %Earth albedo
sigma = 5.67E-8; %[W/m2/K4]
Gs0 = 1371; %[W/m2]
Gp = 237; %[W/m2]

%% Orbit
orb_case = input('Input orbit case to run: ');
switch orb_case
    case 1
        theta_SC = deg2rad(0); %[rad]
        phi_SC = 0; %[rad]
    case 2
        theta_SC = deg2rad(45); %[rad]
        phi_SC = 0; %[rad]
    otherwise
        orb_case = 1;
        theta_SC = deg2rad(0); %[rad]
        phi_SC = 0; %[rad]
end
ha = 400e3; %[m]
r = ha + Re; %[m]
T0 = 2*pi/Re*sqrt(r^3/g); %[s]
omega0 = 2*pi/T0; %[rad/s]
gamma0 = 0; %[rad]
gamma_f = @(t) gamma0 + 2*pi/T0.*t; %[rad]
beta_f = @(gamma) abs(wrapToPi(gamma)); %[rad]
up_f = @(theta,phi) -[sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];
us_f = @(gamma,up) Rmatrix(2,theta_SC).'*[-sin(gamma); 0; cos(gamma)];
delta = acos(Re/r);
eclipse_flag = @(gamma) (pi-abs(wrapToPi(gamma)))<delta;

%% Solver config
%these parameters will define the duration fo the simulation, if tf=10*T0
%the simulation will run for 10 orbits

dt = 1; %[s]
t0 = 0;
%tf = 10*T0; %T0*77.85; %5 days, with 15,57 orbits per day aproximately
tf = 432000; %time to complete five days
ti = (t0:dt:tf).';

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
 
 
 %leds power generation
 SCledqgen = zeros(1,length(ti)); 
 %SCledqgen (1,44300:44350)=50; %W
 SCledqgen (1,86400:86410)=50; %W dato cotejado con juan (aproximacion)
 
 %radio mode, heat applied to trasceiver and antenas?
 SCradiomode = zeros(1,length(ti)); 
 SCradiomode (1,25200:25500) = 5; %W at day 1 supuestos
 SCradiomode (1,115200:115500) = 5; %W at day 2
 SCradiomode (1,435600:435900) = 5; %W at day 5
 SCradiomode (1,460800:461100) = 5; %W at day 5

 %reorientation mode 
  SCreomode = zeros(1,length(ti));
  SCreomode (1,85800:86400) = 5; %W at day 1 (i should ask someone about this shit)
 
  
 %% Try to read thy excel
  
%IMMA = readmatrix('Thermal_Data.xlsx')
%opts = detectImportOptions('Thermal_Data.xlsx');
%opts = xmlImportOptions('Conductances_between_nodes')
%preview('Thermal_Data.xlsx',opts)
  
%% Spacecraft data
A = 0.02; %Areas [m^2]
W = sqrt(A); %Width [m]

% Nodes
N = 29;
SC = struct([]);

%structure nodes 
SC(1).name = 'Frame 1, +X';
SC(1).m = 0.04028;
SC(1).Cp = 800; %consultar 
SC(1).A = 0.0052; %calculo a mano aprox
SC(1).n = [1;0;0];
SC(1).e = 0.8;
SC(1).a = 0.8;
SC(1).coupling = [0, 3, 4, 5, 6]; 
SC(1).qgen=0;


SC(2).name = 'Frame 2, -X';
SC(2).m = 0.04028;
SC(2).Cp = 800;
SC(2).A = 0.0052; %calculo a mano aprox
SC(2).n = [-1;0;0];
SC(2).e = 0.8;
SC(2).a = 0.8;
SC(2).coupling = [0, 3, 4, 5, 6];
SC(2).qgen=0;


SC(3).name = 'Beam -Y +Z';
SC(3).m = 0.007291;
SC(3).Cp = 800;
SC(3).A = 0.0021; %calculo a mano aprox
SC(3).n = [0;-1;1];
SC(3).e = 0.8;
SC(3).a = 0.8;
SC(3).coupling = [0]; 
SC(3).qgen=0;


SC(4).name = 'Beam +Y +Z';
SC(4).m = 0.007291;
SC(4).Cp = 800;
SC(4).A = 0.0021; %calculo a mano aprox
SC(4).n = [0;1;1];
SC(4).e = 0.8;
SC(4).a = 0.8;
SC(4).coupling = [0];
SC(4).qgen=0;


SC(5).name = 'Beam +Y -Z';
SC(5).m = 0.006425;
SC(5).Cp = 800;
SC(5).A = 0.0021; %calculo a mano aprox
SC(5).n = [0;1;-1];
SC(5).e = 0.8;
SC(5).a = 0.8;
SC(5).coupling = [0];
SC(5).qgen=0;


SC(6).name = 'Beam -Y -Z';
SC(6).m = 0.006425;
SC(6).Cp = 800;
SC(6).A = 0.0021; %calculo a mano aprox
SC(6).n = [0;-1;-1];
SC(6).e = 0.8;
SC(6).a = 0.8;
SC(6).coupling = [0]; 
SC(6).qgen=0;

%payload components 

SC(7).name = 'panel +Z, payload plate ';
SC(7).m = 0.02192;
SC(7).Cp = 800; %check
SC(7).A = 0.0096;
SC(7).n = [0;0;1];
SC(7).e = 0.8;
SC(7).a = 0.8;
SC(7).coupling = [0, 1, 2, 3, 4]; 
SC(7).qgen=0;


SC(8).name = 'LED 1 ';
SC(8).m = 0.0018;
SC(8).Cp = 800; %check
SC(8).A = 5.7600e-04;
SC(8).n = [0;0;1];
SC(8).e = 0.8; %check
SC(8).a = 0.8;
SC(8).coupling = [0, 7];
SC(8).qgen= 0; %W


SC(9).name = 'LED 2';
SC(9).m = 0.0018;
SC(9).Cp = 800; %check
SC(9).A = 5.7600e-04;
SC(9).n = [0;0;1];
SC(9).e = 0.8; %check
SC(9).a = 0.8;
SC(9).coupling = [0, 7]; 
SC(9).qgen=0; %W


SC(10).name = 'LED 3';
SC(10).m = 0.0018;
SC(10).Cp = 800; %check
SC(10).A = 5.7600e-04;
SC(10).n = [0;0;1];
SC(10).e = 0.8; %check
SC(10).a = 0.8;
SC(10).coupling = [0, 7]; 
SC(10).qgen=0; %W


SC(11).name = 'LED 4 ';
SC(11).m = 0.0018; %0.2192;
SC(11).Cp = 800; %check
SC(11).A = 5.7600e-04; %0.0096;
SC(11).n = [0;0;1];
SC(11).e = 0.8;
SC(11).a = 0.8;
SC(11).coupling = [0, 7];
SC(11).qgen=0;


SC(12).name = 'Node 12, BMS'; % Sub eps low power, Battery managment, this goes next to the batteries , not on payload plate
SC(12).m = 0.02375;
SC(12).Cp = 600; %isolatros should have more right?
SC(12).A = 0.01; %?
SC(12).coupling = [12:19 ,4, 5]; 
SC(12).qgen=0; 


SC(13).name = 'Node 13, Batteries Isolation';
SC(13).m = 0.125; %yoquese
SC(13).Cp = 800; %isolation should have more right?
SC(13).A = 0.01; %?
SC(13).coupling = [3, 4, 5, 6]; 
SC(13).qgen=0; 

SC(14).name = 'Node 14, LiFePo4 1';
SC(14).m = 0.039; 
SC(14).Cp = 600;
SC(14).A = 0.01; %?
SC(14).coupling = [13]; 
SC(14).qgen=0.5; %??


SC(15).name = 'Node 15, LiFePo4 3';
SC(15).m = 0.039;
SC(15).Cp = 600;
SC(15).A = 0.1; %?
SC(15).coupling = [13];
SC(15).qgen=0.5; %?? 


SC(16).name = 'Node 16, LiFePo4 2';
SC(16).m = 0.039;
SC(16).Cp = 600;
SC(16).A = 0.1; %?
SC(16).coupling = [13];
SC(16).qgen=0.5; %?? 

SC(17).name = 'Node 17, LiFePo4 4';
SC(17).m = 0.039;
SC(17).Cp = 600;
SC(17).A = 0.1; %?
SC(17).coupling = [13];
SC(17).qgen=0.5; %?? 


SC(18).name = 'Node 18, Cell Litio 1';
SC(18).m = 0.05; %checked
SC(18).Cp = 600;
SC(18).A = 0.1; %?
SC(18).coupling = [13];
SC(18).qgen=0.1; %solo funca al principio me cago en la puta


SC(19).name = 'Node 19, Cell Litio 2';
SC(19).m = 0.05;
SC(19).Cp = 600;
SC(19).A = 0.1; %?
SC(19).coupling = [13];
SC(19).qgen=0; %supuestamente añadido


SC(20).name = 'Node 20, LoMo';
SC(20).m = 0.03849;
SC(20).Cp = 600; %directed by quentin inventino
SC(20).A = 0.0086; %HMM ?
SC(20).coupling = [3, 4, 5, 6];
SC(20).qgen=0.5; %god knows

SC(21).name = 'Node 21, Transceiver(Tx-Rx)';
SC(21).m = 0.075;
SC(21).Cp = 600; %directed by quentin inventino
SC(21).A = 0.0086; %HMM ?
SC(21).coupling = [3, 4, 5, 6];
SC(21).qgen=0.5;


SC(22).name = 'Node 22, Antena module';
SC(22).m = 0.085;
SC(22).Cp = 600;
SC(22).coupling = [1, 2, 5, 6];
SC(22).qgen=2; 
SC(22).k = 15;
SC(22).L = 0.7; %ancho
SC(22).A = 0.08; %not sure


SC(23).name = 'Node 23, OBC';
SC(23).m = 0.024; %putos gramos cuidado unidades
SC(23).Cp = 600;
SC(23).A = 0.33; 
SC(23).coupling = [3, 4, 5, 6, 21];
SC(23).qgen=1.5; 
%SC(1).qgenx= SCled4qgen; %W

SC(24).name = 'Node 24, Solar Panel (Y)';
SC(24).m = 0.053;
SC(24).Cp = 600;
SC(24).coupling = [1, 2, 4, 5];
SC(24).qgen=2.4; 
SC(24).e = 0.8;
SC(24).a = 0.8;
SC(24).k = 15;
SC(24).L = 0.0022; %ancho unidades? (2.2 mm)
SC(24).A = 0.008094; %m2?
SC(24).n = [0;1;0]; 
%SC(1).qgenx= SCled4qgen; %W

SC(25).name = 'Node 25, Solar Panel (-Y)';
SC(25).m = 0.044;
SC(25).Cp = 600;
SC(25).coupling = [1, 2, 3, 6];
SC(25).qgen=2.4; 
SC(25).e = 0.8;
SC(25).a = 0.8;
SC(25).k = 15;
SC(25).L = 0.0022; %ancho unidades? (2.2 mm)
SC(25).A = 0.008094; %m2?
SC(25).n = [0;-1;0]; 
%SC(1).qgenx= SCled4qgen; %W

SC(26).name = 'Node 26, Solar Panel (X)';
SC(26).m = 0.045;
SC(26).Cp = 600;
SC(26).coupling = [1];
SC(26).qgen=2.4; 
SC(26).e = 0.8;
SC(26).a = 0.8;
SC(26).k = 15;
SC(26).L = 0.0022; %ancho unidades? (2.2 mm)
SC(26).A = 0.008094; %m2?
SC(26).n = [1;0;0]; 
%SC(26).qgenx= SCled4qgen; %W

SC(27).name = 'Node 27, Solar Panel (-X)RBF';
SC(27).m = 0.053;
SC(27).Cp = 600;
SC(27).coupling = [2];
SC(27).qgen=2.4; 
SC(27).e = 0.8;
SC(27).a = 0.8;
SC(27).k = 15;
SC(27).L = 0.0022; %ancho unidades? (2.2 mm)
SC(27).A = 0.008094; %m2?
SC(27).n = [-1;0;0]; 
%SC(27).qgenx= SCled4qgen; %W

SC(28).name = 'Node 28, Solar Panel (-Z(antenna))';
SC(28).m = 0.0575;
SC(28).Cp = 600;
SC(28).coupling = [1, 2 ,5, 6, 22];
SC(28).qgen=2.4; 
SC(28).e = 0.8;
SC(28).a = 0.8;
SC(28).k = 15;
SC(28).L = 0.0022; %ancho unidades? (2.2 mm)
SC(28).A = 0.008447; %m2?
SC(28).n = [0;0;-1]; 
%SC(28).qgenx= SCled4qgen; %W

SC(29).name = 'Node 29, EPS de alta'; 
SC(29).m = 0.04;
SC(29).Cp = 600; 
SC(29).A = 0.01; %?
SC(29).coupling = [3, 4, 5, 6, 8 ,9 ,10, 11]; 
SC(29).qgen=0; 
%SC(29).qgenx= SCled4qgen; %W


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
opts = spreadsheetImportOptions("NumVariables", 29);

% Specify sheet and range
opts.Sheet = "Conductances_between_nodes";
opts.DataRange = "B5:AD33";

% Specify column names and types
opts.VariableNames = ["node1", "node2", "node3", "node4", "node5", "node6", "node7", "node8", "node9", "node10", "node11", "node12", "node13", "node14", "node15", "node16", "node17", "node18", "node19", "node20", "node21", "node22", "node23", "node24", "node25", "node26", "node27", "node28", "node29"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Import the data
ThermalDataS1 = readtable("Thermal_Data.xlsx", 'Sheet',"Conductances_between_nodes");

%% Convert to output type
ThermalDataS1=ThermalDataS1(1:end,2:end); %Crop text headers
ThermalDataS1 = table2array(ThermalDataS1);

%% Clear temporary variables
clear opts
K = ThermalDataS1;

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
    A = [C,S;-S,C];
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