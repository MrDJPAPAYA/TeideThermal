clc;clear;close all;
set(groot,'DefaultAxesFontSize',11);
set(groot,'DefaultTextInterpreter','latex');

%% Thermal Subsytem
%Property of TESSERA 
%Base code created by Leyre Hernández Palacios
%Modified for TESSERA mission ASTRID by Javier González Vilar
%contact: 100445014@alumnos.uc3m.es


%% Constants
Tc = 32+273.15; %[K]
Re = 800; %[m] Dydimos radious
g = 0.001; %[m/s2]
a = 0.15; %dydimos albedo
sigma = 5.67E-8; %[W/m2/K4]

load('Dist2Sun.mat');
d_sun=2; %SC2Sun_AU; %AU
Gs0 = 1366;% *(1/d_sun)^2; %[W/m2] roughly reduced
Gp = 0; %[W/m2] supposed

orb_case=1;
%% Orbit
%orb_case = input('Input orbit case to run: ');
switch orb_case
    otherwise 
        theta_SC = deg2rad(0); %[rad]
        phi_SC = 0; %[rad]
%     case 2
%         theta_SC = deg2rad(45); %[rad]
%         phi_SC = 0; %[rad]
%     otherwise
%         orb_case = 1;
%         theta_SC = deg2rad(0); %[rad]
%         phi_SC = 0; %[rad]
end
ha = 1000; %[m]
r = ha + Re; %[m]
 T0 = 3600*24*365.25; %2*pi/Re*sqrt(r^3/g); %[s]
 omega0 = 2*pi/T0; %[rad/s]º
 gamma0 = 0; %[rad]
 gamma_f = @(t) gamma0 + 2*pi/T0.*t; %[rad]
 beta_f = @(gamma) abs(wrapToPi(gamma)); %[rad]
 up_f = @(theta,phi) -[sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)]; %sun pointing vector, sun sycronous
 us_f = @(gamma,up)[0;1;0]; %@(gamma,up) Rmatrix(2,theta_SC).'*[-sin(gamma); 0; cos(gamma)];
 delta = acos(Re/r);
 eclipse_flag = @(gamma) 0;%(pi-abs(wrapToPi(gamma)))<delta;

%% Solver config
dt = 1; %[s] i think if should be days now
t0 = 0;
tf = 10*T0;
ti =t_days*24*60*60; %(t0:dt:tf).';

%% Spacecraft data
A1 = 0.1026; %Area [m^2]
A2 = 0.051; %Area [m^2]
W = sqrt(A1); %Width [m]

% Nodes
N = 22;
SC = struct([]);
for i=1:N
SC(i).qgen=0;
end

SC(1).name = 'Node 1, side face +X';
SC(1).m = 0.5;
SC(1).Cp = 887;
SC(1).e = 0.8;
SC(1).a = 0.88;
SC(1).k = 15;
SC(1).L = 0.454;
SC(1).A = A1;
SC(1).n = [1;0;0]; %vector normal a la superficie
SC(1).coupling = [0, 4, 3, 6, 5];%checked

SC(2).name = 'Node 2, back face -X';
SC(2).m = 0.5;
SC(2).Cp = 887;
SC(2).e = 0.8;
SC(2).a = 0.88;
SC(2).k = 15;
SC(2).L = 0.454;
SC(2).A = A1;
SC(2).n = [-1;0;0]; 
SC(2).coupling = [0, 3, 4, 5, 6]; %checked

SC(3).name = 'Node 3, side face +Y';
SC(3).m = 0.5;
SC(3).Cp = 887;
SC(3).e = 0.8;
SC(3).a = 0.88;
SC(3).k = 15;
SC(3).L = 0.454;
SC(3).A = A1;
SC(3).n = [0;1;0]; 
SC(3).coupling = [0, 1, 2, 5, 6];%checked

SC(4).name = 'Node 4, side face -Y';
SC(4).m = 0.5;
SC(4).Cp = 887;
SC(4).e = 0.8;
SC(4).a = 0.88;
SC(4).k = 15;
SC(4).L = 0.454;
SC(4).A = A1;
SC(4).n = [0;-1;0]; %posicion relativa?
SC(4).coupling = [0, 1, 2, 5, 6];

SC(5).name = 'Node 5, top face +Z';
SC(5).m = 0.4;
SC(5).Cp = 887;
SC(5).e = 0.8;
%SC(5).R = (1-SC(5).e)/(A*SC(5).e); what the fuck is that
SC(5).a = 0.88;
SC(5).L = 0.226;
SC(5).A = A2;
SC(5).n = [0;0;1];
SC(5).coupling = [0, 1, 2, 3, 4, 7, 8]; %checked

SC(6).name = 'Node 6, top face -Z';
SC(6).m = 0.4;
SC(6).Cp = 887;
SC(6).e = 0.8;
SC(6).a = 0.88;
SC(6).L = 0.226;
SC(6).A = A2;
SC(6).n = [0;0;-1];
SC(6).coupling = [0, 1, 2, 3, 4];%checked

SC(7).name = 'Node 7, Deployable panel 1';
SC(7).m = 0.8;
SC(7).Cp = 887;
SC(7).A = 0.77*4;
SC(7).n = [0;1;0];
SC(7).e = 0.8;
SC(7).a = 0.88;
SC(7).coupling = [0,1];%checked

SC(8).name = 'Node 8, Deployable panel 2';
SC(8).m = 0.8;
SC(8).Cp = 887;
SC(8).A = 0.77*4;
SC(8).n = [0;1;0];
SC(8).e = 0.8;
SC(8).a = 0.88;
SC(8).coupling = [0,2];%checked

SC(9).name = 'Node 9, OBC';
SC(9).m = 0.076;
SC(9).Cp = 600;
SC(9).A = (96*90)*(10^-3)^2;
SC(9).coupling = [1, 3, 10, 15]; 
SC(9).qgen=0.4; 

SC(10).name = 'Node 10, Transceiver';
SC(10).m = 0.180;
SC(10).Cp = 600; %directed by quentin inventino
SC(10).A = 0.10; %que area es esta, que opero yo con ella bro xddd
SC(10).coupling = [1, 3, 9, 13, 15]; %checked
SC(10).qgen=13;

SC(11).name = 'Node 11, Payload';
SC(11).m = 1.1;
SC(11).Cp = 600;
SC(11).A = 0.1;
SC(11).coupling = [2, 4, 5];  %checked
SC(11).qgen=0; %5.7;

SC(12).name = 'Node 12, batteries with PDU';
SC(12).m = 1.3;
SC(12).Cp = 600;
SC(12).A = 0.1;
SC(12).coupling = [2, 4, 16]; %checked
SC(12).qgen=0.5; %no info in the datasheet,very low generation take care mi man
SC(13).name = 'Node 13, HPA';
SC(13).m = 0.283;
SC(13).Cp = 600;
SC(13).A = 0.1; %?
SC(13).coupling = [1, 3, 10];%check
SC(13).qgen=1.2; 

SC(14).name = 'Node 14, Omnidirectional antenna';
SC(14).m = 0.130;
SC(14).Cp = 600;
SC(14).e = 0.8;
SC(14).a = 0.8;
SC(14).k = 15;
SC(14).L = 0.454;
SC(14).A = 0.08;
SC(14).n = [0;1;0]; %vector normal a la superficie, where is antena bro
SC(14).coupling = [0, 3]; %checked
SC(14).qgen=3;

SC(15).name = 'Node 15, 4RW array';
SC(15).m = 0.665;
SC(15).Cp = 600;
SC(15).A = 0.1; %?
SC(15).coupling = [1, 3, 9];%checked
SC(15).qgen=14; 

SC(16).name = 'Node 16, Star tracker';
SC(16).m = 0.35;
SC(16).Cp = 600;
SC(16).A = 0.1; %?
SC(16).coupling = [2, 4, 12];%checked
SC(16).qgen=1.2; 

SC(17).name = 'Node 17, Back-up Star tracker';
SC(17).m = 0.085;
SC(17).Cp = 600;
SC(17).A = 0.1; %?
SC(17).coupling = [2, 4, 6];%checked
SC(17).qgen=1.2; 

%to add from AOCS, low mass, no compsumtion specified
% 2	Fine Sun sensor
% 4	Coarse Sun sensor
% 1	Gyro
% 1	Backup Gyro
% 4	Acc

SC(18).name = 'Node 18, Trushter';
SC(18).m = 1.9; %checked
SC(18).Cp = 600;
SC(18).A = 0.1; %?
SC(18).coupling = [1, 3, 6, 21];%checked
SC(18).qgen=0.1; %solo funca al principio me cagoi en la puta


SC(19).name = 'Node 19, Tank 1';
SC(19).m = 4.5;
SC(19).Cp = 600;
SC(19).A = 0.1; %?
SC(19).coupling = [1, 4];%checked
SC(19).qgen=0; %supuestamente añadido 

SC(20).name = 'Node 20, Tank 2';
SC(20).m = 4.5;
SC(20).Cp = 600;
SC(20).A = 0.1; %?
SC(20).coupling = [2, 3];%checked
SC(20).qgen=0; %supuestamente añadido 


SC(21).name = 'Node 21, MEMs'; %son las imus los mems? quien soy yo?
SC(21).m = 0.05;
SC(21).Cp = 600;
SC(21).A = 0.1; %?
SC(21).coupling = [1, 3];%checked
SC(21).qgen=0.5; 

SC(22).name = 'Node 22, Tank He';
SC(22).m = 0.21;
SC(22).Cp = 600;
SC(22).A = 0.1; %?
SC(22).coupling = [1, 3, 5];%checked
SC(22).qgen=0; %supuestamente añadido



%H45 = sigma/(SC(4).R + 1/A + SC(5).R);

% Conductances   FIXEAR, me duele la cabeza, ahora no :(, se me acaba el
% tiempo me quiero morir
K = zeros(N); %[W/K]
%K(1,2) = SC(1).k*A*2*SC(1).L;
%K(2,3) = SC(3).k*A*2/SC(3).L;
%K(3,4) = 1/(SC(3).L/2/SC(3).k/A + SC(4).L/2/SC(4).k/A);
%K(1,6:9) = 4e-2;
%K(2,6:9) = 0.2;
%K(4,6:9) = 4e-4;
%K(6,[7,9]) = 0.4;
%K(7,8) = 0.4;
%K(8,9) = 0.4;

K(1,3:6)=0.4; %posiblemente aumentable
K(1,[13,10,9,15,21])=0.4;
K(1,18)=1.2;
K(1,19)=0.4;
K(1,22)=0.1;

K(2,3:6)=0.4;
K(2,[11,12,16,17])=0.4;


K(3,[1,2,5,6,10,13,9,15,21])=0.4;
K(3,14)=1.6;
K(3,20)=0.6;
K(3,18)=1.2;
k(3,22)=1;


K(4,[1,2,5,6])=0.4;
K(4,[11,12,16,17])=0.4;
K(4,19)=0.2;


K(5,1:4)=0.4;
K(5,[7,8])=0.4;
K(5,[11,22])=0.2;


K(6,1:4)=0.4;
K(6,18)=0.4;

K(7,3)=0.4;
K(8,3)=0.4;

K = K.'+K; 

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