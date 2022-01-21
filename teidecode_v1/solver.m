%% Thermal Subsytem

%% Load data
data;

%% Solver
Ti = ones(N,1)*300; %[K]
%Ti(2) = Tc;
D = ([SC.m].*[SC.Cp])/dt.*eye(N);
%D(2,2) = 0;

t = find(ti>9*T0);
T = zeros(N,length(t));
q = zeros(N,N,length(t));
Q = zeros(N,length(t));

%% Loop
j = 0;
for k = 1:length(ti)
    gamma = gamma_f(ti(k));
    up = up_f(theta_SC,phi_SC);
    us = us_f(gamma,up);
    Gs = Gs0*(~eclipse_flag(gamma));
    beta = beta_f(gamma);
    F = albedo_F_f(beta);
    B = zeros(N,1);
    K0 = zeros(N,1);
    
     %LEDs On/Off condition
     
%       if 0<ti(k)<10
%          SC(8).qgen=100;
%          SC(9).qgen=100;
%          SC(10).qgen=100;
%          SC(11).qgen=100;
%       end
%       
%       if 10<ti(k)<5548 %final primera orbita
%          SC(8).qgen=0;
%          SC(9).qgen=0;
%          SC(10).qgen=0;
%          SC(11).qgen=0; 
%       end
      
%       if 5548<ti(k)<5558 
%          SC(8).qgen=100;
%          SC(9).qgen=100;
%          SC(10).qgen=100;
%          SC(11).qgen=100;
%       end

      if 54300<ti(k)<54301
         SC(8).qgen=100;
         SC(9).qgen=100;    
         SC(10).qgen=100;
         SC(11).qgen=100;
      end
      
      if 54301<ti(k)
         SC(8).qgen=0;
         SC(9).qgen=0;
         SC(10).qgen=0;
         SC(11).qgen=0;
       end
      
    for i = 1:N
                       
        if ismember(0,SC(i).coupling) %External loads
            cos_s = us.'*SC(i).n; cos_s = cos_s*(cos_s>0);
            cos_p = up.'*SC(i).n; cos_p = cos_p*(cos_p>0);
             

            B(i) = SC(i).A*(SC(i).a*Gs*(cos_s + cos_p*a*F) + SC(i).e*cos_p*Gp) + SC(i).qgen;
            K0(i) = SC(i).A*SC(i).e*sigma*Ti(i)^3;
            
        end
    end
    %B(2) = Tc; %Boundary condition/Temperature constraint in node 2
    %K(4,5) = H45*(Ti(5)^2+Ti(4)^2)*(Ti(5)+Ti(4)); K(5,4) = K(4,5); %Radiative conductances
    C = -K + eye(N).*(sum(K,2) + K0);
    %C(2,:) = 0; C(2,2) = 1; %Boundary condition/Temperature constraint in node 2
    
    
    Ti = (D+C)\(D*Ti+B);
    
    if ti(k)>9*T0
        j = j + 1;
        T(:,j) = Ti;
        
        for i = 1:N
            q(i,:,j) = -K(i,:).*(Ti-Ti(i)).';
            q(i,i,j) = B(i)*(i~=2) - K0(i)*Ti(i);
        end
        Q(:,j) = sum(q(:,:,j));
    end
end