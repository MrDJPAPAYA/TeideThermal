%% Thermal Subsytem

%% Solver
Ti = ones(N,1)*InitialT; %[K] Define starting temperature
%Ti = innit; %[K] Define starting temperature
%Ti(2) = Tc;
D = ([SC.m].*[SC.Cp]).'; %Vector with heat capacities (N sized)
%D(2,2) = 0;
recordt = 1; % Time at wich it will start recording data
t = find(ti>recordt);  
T = zeros(N,length(t));
Results.t = find(ti>recordt)*dt; 
Results.T = zeros(N,length(t));
Results.HeatGen = zeros(N,length(t));
Results.NetCond = zeros(N,length(t));

if config.SolRad
Results.SolRad = zeros(N,length(t));
end
if config.EnvRad
Results.EnvRad = zeros(N,length(t));
end

%% Main calculation Loop
%Startup

C = K - eye(N).*sum(K,2); %Proper conductance matrix
disp("Main loop startup")
%Loop
j = 0; % Keeps track of the position in the result matrixes.
m = 1; %Keeps track of the current mode
changetime = Mode(1).time; %Time at which it will change mode.
for k = 1:length(ti)

    %Parameter calculation
    nu = nu_f(ti(k)); % Current anomaly
    Rnu = R.nu_f(nu); % Anomaly rotation matrix
    Rpos = R.pos_f(Rnu); % Position rotation matrix
    us = us_f(Rpos); % Sun pointing vector
    Gs = Gs0*(~eclipse_flag(us)); % Solar radiation if not on eclipse
    F = albedo_F_f(acos(cosbeta_f(us))); % Visibility factor

    %Boundary conditions
    SolRad = zeros(N,1); % Radiation from the Sun and Earth
    EnvRad = zeros(N,1); % Radiation dissipated to the enviroment
    for i = 1:N           
        if SC(i).radiates %External loads
            if config.SolRad
            cos_s = us.'*SC(i).n; cos_s = cos_s*(cos_s>0);
            cos_p = up.'*SC(i).n; cos_p = cos_p*(cos_p>0);
            SolRad(i) = SC(i).A*(SC(i).a*Gs*(cos_s + cos_p*a*F) + SC(i).e*cos_p*Gp); % Heat due to albedo, earth IR and sun.
            %SolRad(i) = SC(i).A*(SC(i).e*cos_p*Gp); % Heat due to albedo
            end

            if config.EnvRad
            EnvRad(i) = -SC(i).A*SC(i).e*sigma*(Ti(i)^4-(EnvT)^4); % Heat dissipation via radiation
            end
        end
    end

    while ti(k) > changetime % Changes mode when time is due and defines new time
        m = m + 1;
        changetime = Mode(m).time;
    end

    B = SolRad + EnvRad + Mode(m).HeatGen(:, :); % Define Boundary condition.

    % Define net transfer from conduction
    NetCond = C*Ti; 

    % Update temperature
    Ti = Ti+dt*(B + NetCond)./D;

    %Records data for later usage
    if  ti(k)>recordt
        j = j + 1;
        T(:,j) = Ti; 
        Results.T(:, j) = Ti; % Records temperatures
        Results.HeatGen(:, j) = Mode(m).HeatGen(:, :); % Records heat generation
        Results.NetCond(:, j) = NetCond; % Records net heat conduction
        if config.SolRad
            Results.SolRad(:, j) = SolRad; % Records Radiation to Earth and Sun
        end
        if config.EnvRad
            Results.EnvRad(:, j) = EnvRad; % Records Enviroment radiation
        end
    end
end
disp("Main loop end")