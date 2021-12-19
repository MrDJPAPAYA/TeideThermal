function F = albedo_F_f(beta)
    %Approximate albedo view factor at ha=400 km 
    %Created by: Antonio Acosta Iborra. 12 March 2021. Carlos III University of
    %Madrid
    %INPUT ARGUMENT: 
    %beta [rad], angle between solar and planet unit vectors
    %OUTPUT ARGUMENT: 
    %F [-], albedo view factor

    Albedom=[0	1.15
        30	1
        60	0.65
        70	0.5
        80	0.3
        90	0.14
        100	4e-2
        110	7e-4
        120 0
        180 0];

    Albedobetav=Albedom(:,1)*pi/180;
    AlbedoFv=Albedom(:,2);

    F=interp1(Albedobetav,AlbedoFv,beta,'pchip');
end