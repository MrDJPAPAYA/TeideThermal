close all;

%% Thermal Subsytem
%Property of TEIDESAT
%Base code created by Leyre Hernández Palacios
%Modified for TEIDESAT-I by Javier González Vilar
%contact: teidesat13@ull.edu.es
%% Solve problem
solver;

%% Plot results
figdir = ['../Figures/Orb_case' sprintf('%i',orb_case) '_'];

%% Temperature evolution
figure(1);
colours = lines(N);
markers = {'-o','-s','-d','-v','-^','->','-<','-x','-*','-o','-*','-o','->','-o','-s','-d','-v','-x','-*','-o','-*','-o','-s','-d','-v','-^','->','-<','-x'};
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
title(['Temperature evolution' newline' '$\theta_{sc} = ' ...
    sprintf('%.0f',rad2deg(theta_SC)) '^{\circ}$ and $\phi_{sc} = ' ...
    sprintf('%.0f',rad2deg(phi_SC)) '^{\circ}$']);
%saveas(gcf,[figdir 'Tev' '.png']);

%% Maximum temperature
Tmax = max(T,[],'all');
[n_Tmax, t_Tmax] = find(T==Tmax);
t_Tmax = t(t_Tmax);

hold on;
plot((t_Tmax-t(1))/T0,Tmax-273.15,'ro','HandleVisibility','off');
text((t_Tmax-t(1))/T0+0.17,Tmax-290,...
    ['Maximum temperature: $' sprintf('%0.2f', Tmax-273.15) '^{\circ}C$' ...
    newline 'found at node ' sprintf('%i', n_Tmax') ', ' SC(n_Tmax).name ',' ...
    newline 'and orbit angle $\gamma = ' ...
    sprintf('%0.2f', rad2deg(mod(gamma_f(t_Tmax),2*pi))) '^{\circ}$.']);
%saveas(gcf,[figdir 'Tmax' '.png']);

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
%saveas(gcf,[figdir 'Tout' '.png']);

% %% 3D plot
% k = 3;
% for i = round(linspace(1,length(t),length(xt)))
%     k = k + 1;
%     up = up_f(theta_SC,phi_SC);
%     us = us_f(xt(k-3),up);
%     figure(k);
%     Sat3DnodeT_v1_f(W,L,SC(5).L,T(:,i),us,up);
%     sgtitle(['Outer temperature distribution for $\gamma$=' xtl{k-3} ...
%         newline '$\theta_{sc} = ' sprintf('%.0f',rad2deg(theta_SC)) ...
%         '^{\circ}$ and $\phi_{sc} = ' sprintf('%.0f',rad2deg(phi_SC)) '^{\circ}$'],...
%         'Interpreter','latex');
%     %saveas(gcf,[figdir '3Dgamma' sprintf('%i',k-3) '.png']);
% end
% 


%% Repeat?
% switch input('Run another case? ','s')
%     case {'1', 'y', 'Y', 'yes', 'Yes', 'Sure!', 'true', 'T', 'True'}
%         plotter;
%     otherwise
%         if input('Close plots? ')
%             close all;
%         end
%         if input('Clear results and cmd? ')
%             clear;clc;
%         end
% end
