close all;

%% Thermal Subsytem
%Property of TEIDESAT
%Base code created by Leyre Hernández Palacios
%Modified for TEIDESAT-I by Javier González Vilar
%contact: teidesat13@ull.edu.es

%%
%Self notes
%probar el codigo con un uninodo de un bloque solido a ver si se aproxima
%al resultado dado
%corregir parametros deentrada, leerlos desde un excel, conseguido con la K
%falta con el struct

%% Path Config
path = "../saves/";


%% Load or solve
action = input("Choose action (load/solve): ", "s");

%% Solve problem and plot
if action == "solve"
%Solve problem
solver;
%Save data
if input("Save data? (y/n): ", "s")== "y"
    data = {N, T, T0, t};
    filename = input("Save as: ", "s");
    save(append(path,filename), "data")
end
%Plot data
plotter;
end
%% Load data and plot
if action == "load"
    %Load and unpack
    filename = input("Enter file name: ", "s");
    if isfile(append(path,filename,".mat"))
    load(append(path,filename), 'data');
    N = data{1};
    T = data{2};
    T0 = data{3};
    t = data{4};
    %plot
    plotter;
    %If file couldnt be found
    else
    disp(['File ' filename ' does not exist.']);
    end
end
