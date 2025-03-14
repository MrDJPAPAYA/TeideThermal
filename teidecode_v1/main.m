close all; clear

%% Thermal Subsytem
%Property of TEIDESAT
%Base code created by Leyre Hernández Palacios
%Modified for TEIDESAT-I by Javier González Vilar
%Further modified by Eduardo Andrés Navarro Santos
%contact: teidesat03@ull.edu.es

%% Load or solve
action = input("Choose action (load/solve): ", "s");

%% Solve problem and plot
if action == "solve"
% Load data
data;
%Solve problem
solver;
%Plot data
plotter;

%% Load data and plot
elseif action == "load"
    %Load and unpack
    path = "../saves/"; %Folder location
    filename = input("Enter file name: ", "s");
    filename = append(filename, ".mat");
    filepath = fullfile(path,filename);
    if isfile(filepath)
    load(filepath);
    %If file couldnt be found
    else
    disp(['File ' filename ' does not exist.']);
    end
else
    error("Invalid Command")
end    

%% Save data
SaveVar = [0, 0];
if input("Save data as mathlab file? (y/n): ", "s")== "y"
    SaveVar(1) = 1;
end
if input("Save data as excel file? (y/n): ", "s")== "y"
    SaveVar(2) = 1;
end
if any(SaveVar)
path = "../saves/"; %Folder location
    filename = input("Save as: ", "s");
    filepath = fullfile(path,filename);
end
if SaveVar(1)
   save(filepath)
end
if SaveVar(2)
   disp("Saving on excel... Please, don't open the file until finalization to avoid errors.")
   SaveResultsToExcel(filepath, 'Temperatures', Results.t, Results.T, {SC.name});
   SaveResultsToExcel(filepath, 'Net heat transfer', Results.t, Results.NetCond, {SC.name});
   SaveResultsToExcel(filepath, 'Heat generation', Results.t, Results.HeatGen, {SC.name});
   if config.SolRad
       SaveResultsToExcel(filepath, 'Heat by sun and albedo', Results.t, Results.SunRad, {SC.name});
   end
   if config.EnvRad
       SaveResultsToExcel(filepath, 'Heat radiated', Results.t, Results.EnvRad, {SC.name});
   end
   disp("Done saving!")
end

%% Save to excel auxiliary function

function SaveResultsToExcel(filename, sheetname, time, result, names)
    % Extract time instances and element names
    timeVector = time;   % Column vector of time instances
    elementNames = names; % Row vector of element names
    
    % Ensure proper dimensions
    timeVector = timeVector(:);  % Force column format
    elementNames = elementNames(:)'; % Force row format
    
    % Iterate over each result type
    resultType = string(sheetname);  % Get current result type (e.g., 'heat', 'temperature')

    % Extract the corresponding result matrix
    dataMatrix = transpose(result);
        
    % Combine time vector with result matrix
    fullMatrix = [timeVector, dataMatrix];

    % Prepare the header row (including time column label)
    headerRow = ['Time (s)', elementNames];

    % Convert to cell array for writing
    outputData = [headerRow; num2cell(fullMatrix)];
   
    % Write to Excel in a new sheet
    writetable(cell2table(outputData), strcat(filename,".xlsx"), "Sheet", resultType, 'WriteVariableNames', false);
end