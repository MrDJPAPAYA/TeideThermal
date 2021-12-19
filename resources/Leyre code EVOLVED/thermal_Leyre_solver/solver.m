%% Thermal Subsytem
%Property of TESSERA 
%Base code created by Leyre Hernández Palacios
%Modified for TESSERA mission ASTRID by Javier González Vilar
%contact: 100445014@alumnos.uc3m.es
%% Load data
data;

%% Solver
Ti = ones(N,1)*300; %[K]
%Ti(2) = Tc;
D = ([SC.m].*[SC.Cp])/dt.*eye(N);
%D(2,2) = 0;

t = ti;% find(ti>9*T0);
T = zeros(N,length(t));
q = zeros(N,N,length(t));
Q = zeros(N,length(t));

%% Loop
j = 0;
for k = 1:length(ti)
    
    gamma = gamma_f(ti(k));
    up = up_f(theta_SC,phi_SC);
    us = us_f(gamma,up);
    Gs = Gs0*(1/SC2Sun_AU(k))^2;%(~eclipse_flag(gamma));
    beta = beta_f(gamma);
    F = albedo_F_f(beta);
    B = zeros(N,1);
    K0 = zeros(N,1);
    for i = 1:N
        if ismember(0,SC(i).coupling) %External loads
            
            
             if k<((1/4)*(length(ti))) %reorientacion of the panel at the initial phase of the trip
                
           SC(7).n=[sin(30);cos(70);0];
           SC(8).n=[sin(30);cos(70);0];   
             else
           SC(7).n=[0;1;0];
           SC(8).n=[0;1;0];       
             end
           
            cos_s = us.'*SC(i).n; cos_s = cos_s*(cos_s>0);
            cos_p = up.'*SC(i).n; cos_p = cos_p*(cos_p>0);
            B(i) = SC(i).A*(SC(i).a*Gs*(cos_s + cos_p*a*F) + SC(i).e*cos_p*Gp);
            K0(i) = SC(i).A*SC(i).e*sigma*Ti(i)^3;
            
            
            
            K0(7) = SC(7).A*0.5*2*sigma*Ti(7)^3; %añadido parte de atras de los desplegables
            K0(8) = SC(7).A*0.5*2*sigma*Ti(7)^3;
            
        end
        B(i)=B(i)+SC(i).qgen;
    end
    %B(2) = Tc; %Boundary condition/Temperature constraint in node 2
    %K(4,5) = H45*(Ti(5)^2+Ti(4)^2)*(Ti(5)+Ti(4)); K(5,4) = K(4,5); %Radiative conductances
    C = -K + eye(N).*(sum(K,2) + K0);
    %C(2,:) = 0; C(2,2) = 1; %Boundary condition/Temperature constraint in node 2
    
    Ti = (D+C)\(D*Ti+B);
    
%     if ti(k)>9*T0
      j = j + 1;
      T(:,j) = Ti;
        for i = 1:N
             q(i,:,j) = -K(i,:).*(Ti-Ti(i)).';
            q(i,i,j) = B(i)*(i~=2) - K0(i)*Ti(i);
        end
         Q(:,j) = sum(q(:,:,j));
%     end
end
%% Max and mins Temps
T_MAX=zeros (1,22);
for i=1:22
T_MAX(i)= (max (T(i,:)))-273.15;  
end

 T_MIN=zeros (1,22);
 for i=1:22
 T_MIN(i)= (min (T(i,:))) -273.15;  
 end