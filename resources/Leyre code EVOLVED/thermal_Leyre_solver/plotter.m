close all;

%% Thermal Subsytem
%Property of TESSERA 
%Base code created by Leyre Hernández Palacios
%Modified for TESSERA mission ASTRID by Javier González Vilar
%contact: 100445014@alumnos.uc3m.es

%% Solve problem
solver;

%% Plot results
figdir = ['../Figures/Orb_case' sprintf('%i',orb_case) '_'];

%% Temperature evolution
figure(1);
colours = hsv(N);
markers = {'-o','-s','-d','-v','-^','->','-<','-x','-*','-o','-*','-o','->','-o','-s','-d','-v','-x','-*','-o','-*','-o'};
for i=1:N
    p = plot((t-t(1))/T0,T(i,:)-273.15,markers{i},'Color',colours(i,:),'LineWidth',0.6);
    p.MarkerIndices = 1:500:length(t);
    hold on;
end
legend({SC.name},'Interpreter','latex','Location','eastoutside');
xlabel('$t/T_0$');
ylabel('T [$^{\circ}C$]');
grid on; grid minor;
set(gcf,'Position',[488.2000  296.2000  844.8000  465.6000]);
title(['Temperature evolution in one period' newline' '$\theta_{sc} = ' ...
    sprintf('%.0f',rad2deg(theta_SC)) '^{\circ}$ and $\phi_{sc} = ' ...
    sprintf('%.0f',rad2deg(phi_SC)) '^{\circ}$']);

%% Temperature hot case
figure(2);
colours = hsv(N);
markers = {'-o','-s','-d','-v','-^','->','-<','-x','-*','-o','-*','-o','->','-o','-s','-d','-v','-x','-*','-o','-*','-o'};

for i=1:N
    p = plot((t(1:5000)-t(1))/T0,T(i,1:5000)-273.15,markers{i},'Color',colours(i,1:3),'LineWidth',0.6);
    p.MarkerIndices = 1:500:5000;
    hold on;
end

%legend({SC.name},'Interpreter','latex','Location','eastoutside');
xlabel('$t/T_0$');
ylabel('T [$^{\circ}C$]');
grid on; grid minor;
%set(gcf,'Position',[488.2000  296.2000  844.8000  465.6000]);
%title(['Temperature evolution Hot case' newline' '$\theta_{sc} = ' ...
   % sprintf('%.0f',rad2deg(theta_SC)) '^{\circ}$ and $\phi_{sc} = ' ...
    %sprintf('%.0f',rad2deg(phi_SC)) '^{\circ}$']);

%% Temperature cold case
figure(3);
colours = hsv(N);
markers = {'-o','-s','-d','-v','-^','->','-<','-x','-*','-o','-*','-o','->','-o','-s','-d','-v','-x','-*','-o','-*','-o'};

for i=1:N
    p = plot((t(5000:13052)-t(1))/T0,T(i,5000:13052)-273.15,markers{i},'Color',colours(i,1:3),'LineWidth',0.6);
    p.MarkerIndices = 5000:504:13052;
    hold on;
end

%legend({SC.name},'Interpreter','latex','Location','eastoutside');
xlabel('$t/T_0$');
ylabel('T [$^{\circ}C$]');
grid on; grid minor;
%set(gcf,'Position',[488.2000  296.2000  844.8000  465.6000]);
%title(['Temperature evolution Cold case' newline' '$\theta_{sc} = ' ...
    %sprintf('%.0f',rad2deg(theta_SC)) '^{\circ}$ and $\phi_{sc} = ' ...
    %   sprintf('%.0f',rad2deg(phi_SC)) '^{\circ}$']);

% %% Temperature at outer surface
% figure(2);
% j = []; 
% xt = [0 pi/2 pi 3*pi/2 2*pi];
% xtl = {'$0$', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'};
% for i = 1:N
%     if ismember(0,SC(i).coupling)
%         p = plot(wrapTo2Pi(gamma_f(t(1:end-1))),T(i,1:end-1),markers{i},'Color',colours(i,:),'LineWidth',0.6);
%         j = [j i];
%         p.MarkerIndices = round(linspace(1,length(t)-1,length(xt)));
%         hold on;
%     end
%     legend({SC(j).name},'Interpreter','latex','Location','eastoutside');
% end
% xlabel('$\gamma$ [rad]');
% ylabel('T [K]');
% grid on; grid minor;
% set(gcf,'Position',[488.2000  296.2000  844.8000  465.6000]);
% ax = gca;
% xlim([xt(1) xt(end)]);
% ax.XTick = xt;
% ax.XTickLabel = xtl;
% ax.TickLabelInterpreter = 'latex';
% title(['Outer surface temperature in one period' newline '$\theta_{sc} = ' ...
%     sprintf('%.0f',rad2deg(theta_SC)) '^{\circ}$ and $\phi_{sc} = ' ...
%     sprintf('%.0f',rad2deg(phi_SC)) '^{\circ}$']);


