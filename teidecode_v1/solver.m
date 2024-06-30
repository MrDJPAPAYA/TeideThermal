%% Thermal Subsytem

%% Load data
data;

%% Solver
Ti = ones(N,1)*300; %[K]
%Ti(2) = Tc;
D = ([SC.m].*[SC.Cp])/dt*ones(N,1); %Vector with heat capacities/dt (N sized)
%D(2,2) = 0;

t = find(ti>T0); %find(ti>9*T0); 
T = zeros(N,length(t));
q = zeros(N,N,length(t));
Q = zeros(N,length(t));

%% Main calculation Loop
%Startup
B = zeros(N,1);
P = zeros(N,1);
C = -K + eye(N).*sum(K,2); %Proper conductance matrix
disp("Main loop startup")
%Loop
j = 0;
for k = 1:length(ti)
    %Parameter calculation
    gamma = gamma_f(ti(k));
    up = up_f(theta_SC,phi_SC);
    us = us_f(gamma,up);
    Gs = Gs0*(~eclipse_flag(gamma));
    beta = beta_f(gamma);
    F = albedo_F_f(beta);

    %Boundary conditions
    for i = 1:N
                       
        if ismember(0,SC(i).coupling) %External loads
            cos_s = us.'*SC(i).n; cos_s = cos_s*(cos_s>0);
            cos_p = up.'*SC(i).n; cos_p = cos_p*(cos_p>0);

            B(i) = SC(i).A*(SC(i).a*Gs*(cos_s + cos_p*a*F) + SC(i).e*cos_p*Gp) + SC(i).qgen; %Heat due to albedo and sun (Probably)
            P(i) = SC(i).A*SC(i).e*sigma*Ti(i)^4; %Heat dissipation via radiation
            
            if (8<=i)&&(i<=11) %if the loop is in the leds it shall add the power defined in data
                B(i)=B(i)+SCledqgen(k);
            end
            if (21<=i)&&(i<=23) %if the loop is in the Tx-Tx and the module
                B(i)=B(i)+SCradiomode(k);
            end  
            if i==20 %if the loop is in Lomo
                B(i)=B(i)+SCradiomode(k);
            end
        end
    end
    B = B-P;
    %Lineal heat tranfers ahead
    %B(2) = Tc; %Boundary condition/Temperature constraint in node 2
    %K(4,5) = H45*(Ti(5)^2+Ti(4)^2)*(Ti(5)+Ti(4)); K(5,4) = K(4,5); %Radiative conductances
    %C(2,:) = 0; C(2,2) = 1; %Boundary condition/Temperature constraint in node 2
    
    
    Ti = Ti+(B-C*Ti)/D;%Update temperature

    %Tbh I dont know what this is
    if  ti(k)>T0 %ti(k)>9*T0 % 
        j = j + 1;
        T(:,j) = Ti;
        
        for i = 1:N
            q(i,:,j) = -K(i,:).*(Ti-Ti(i)).';
            q(i,i,j) = B(i)*(i~=2) - P(i);
        end
        Q(:,j) = sum(q(:,:,j));
    end
end
disp("Main loop end")