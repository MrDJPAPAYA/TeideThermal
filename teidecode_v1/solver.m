%% Thermal Subsytem

%% Load data
data;

%% Solver
Ti = ones(N,1)*InitialT; %[K] Define starting temperature
%Ti(2) = Tc;
D = ([SC.m].*[SC.Cp])/dt*ones(N,1); %Vector with heat capacities/dt (N sized)
%D(2,2) = 0;

t = find(ti>0); %find(ti>9*T0); 
T = zeros(N,length(t));
q = zeros(N,N,length(t));
Q = zeros(N,length(t));

%% Main calculation Loop
%Startup

C = -K + eye(N).*sum(K,2); %Proper conductance matrix
disp("Main loop startup")
%Loop
j = 0;
m = 1; %Keeps track of the current mode
changetime = Mode(1).time; %Time at wich it will change mode.
for k = 1:length(ti)
    %Parameter calculation
    gamma = gamma_f(ti(k));
    up = up_f(theta_SC,phi_SC);
    us = us_f(gamma,up);
    Gs = Gs0*(~eclipse_flag(gamma));
    beta = beta_f(gamma);
    F = albedo_F_f(beta);
    B = zeros(N,1);
    P = zeros(N,1);
    %Boundary conditions
    for i = 1:N
                       
        if SC(i).radiates %External loads
            cos_s = us.'*SC(i).n; cos_s = cos_s*(cos_s>0);
            cos_p = up.'*SC(i).n; cos_p = cos_p*(cos_p>0);

            if config.SolRad
            B(i) = SC(i).A*(SC(i).a*Gs*(cos_s + cos_p*a*F) + SC(i).e*cos_p*Gp); %Heat due to albedo and sun
            end

            if config.EnvRad
            P(i) = SC(i).A*SC(i).e*sigma*(Ti(i)^4-(EnvT)^4); %Heat dissipation via radiation
            end
        end
    end
    if k > changetime % Changes mode when time is due and defines new time
        m = m + 1;
        changetime = Mode(m).time;
    end

    B = B + Mode(m).heats(:, :); % Add Heat generation
    B = B-P;
    %Lineal heat tranfers ahead
    %B(2) = Tc; %Boundary condition/Temperature constraint in node 2
    %K(4,5) = H45*(Ti(5)^2+Ti(4)^2)*(Ti(5)+Ti(4)); K(5,4) = K(4,5); %Radiative conductances
    %C(2,:) = 0; C(2,2) = 1; %Boundary condition/Temperature constraint in node 2
    
    
    Ti = Ti+(B-C*Ti)/D;%Update temperature

    %Records stuff for later usage
    if  ti(k)>0 %ti(k)>9*T0 % %Set to record data, when sim time is greater than x number of orbits
        j = j + 1;
        T(:,j) = Ti; %Records temperatures in big matrix
        
        for i = 1:N %For every node
            q(i,:,j) = -K(i,:).*(Ti-Ti(i)).';
            q(i,i,j) = B(i)*(i~=2) - P(i);
        end
        Q(:,j) = sum(q(:,:,j)); %Definetly creates something related to heat... Maybe I wrote this line myself and now can't rememberwhat it does?
    end
end
disp("Main loop end")