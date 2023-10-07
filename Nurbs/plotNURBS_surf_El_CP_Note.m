function plotNURBS_surf_El_CP_Note(p,q,U,V,CP)
% plots the surface, elements and control points

[X,Y,Z] = create_surf(p,q,U,V,CP);

% geometry
figure('Color','w','Units','normalized','Outerposition',[0 0 0.40 0.95],'visible','on'); axis equal
% figure('Color','w','Units','normalized','Outerposition',[0 0 0.33 0.91],'visible','on'); axis equal
% figure('Color','w','Units','normalized','Outerposition',[0 0 0.51 0.95],'visible','on'); axis equal
surf(X,Y,Z,'FaceColor','none','EdgeColor','none','LineStyle','--', 'LineWidth',1.2);
% xlabel('x'); ylabel('y'); zlabel('z');
hold on;

%% Square Plate
view(2)
xlim([0, 1]); ylim([0, 1]);
set(gca, 'YTick', [0.95], 'XTick', [0.95], 'TickLength', [0.001 0.001], 'TickLabelInterpreter', 'latex')
pgca = gca;
pgca.YTickLabel = ['$a$'];
pgca.XTickLabel = ['$a$'];
pgca.XRuler.TickLabelGapOffset = 0;
pgca.YRuler.TickLabelGapOffset = 10;
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
set(gca, 'box', 'off', 'GridLineStyle', 'none')
axp = get(gca,'Position');
axis off
annotation('arrow', [axp(1)+0.028 axp(1)+axp(3)+0.04], [axp(2) axp(2)]+0.11, 'linewidth', 1.5);
annotation('arrow', [axp(1) axp(1)]+0.032, [axp(2)+0.11 axp(2)+axp(4)-0.05], 'linewidth', 1.5);
annotation('textbox', [0.89 0.12 0.1 0.1],'EdgeColor','none','String','\textbf{$x$}','FitBoxToText','on','fontsize',26,'Interpreter','latex');
annotation('textbox', [0.07 0.80 0.1 0.1],'EdgeColor','none','String','\textbf{$y$}','FitBoxToText','on','fontsize',26,'Interpreter','latex');
annotation('textbox', [0.07 0.12 0.1 0.1],'EdgeColor','none','String','\textbf{$O$}','FitBoxToText','on','fontsize',26,'Interpreter','latex');
set(gca,'fontsize', 26)

%% Rectangular Plate
% view(2)
% xlim([0, 1]); ylim([0, 1]);
% set(gca, 'YTick', [1.45], 'XTick', [0.95], 'TickLength', [0.001 0.001], 'TickLabelInterpreter', 'latex')
% pgca = gca;
% pgca.YTickLabel = ['$b$'];
% pgca.XTickLabel = ['$a$'];
% pgca.XRuler.TickLabelGapOffset = 0;
% pgca.YRuler.TickLabelGapOffset = 10;
% set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
% set(gca, 'box', 'off', 'GridLineStyle', 'none')
% axp = get(gca,'Position');
% annotation('arrow', [axp(1) axp(1)+axp(3)+0.07], [axp(2) axp(2)]+0.026, 'linewidth', 1);
% annotation('arrow', [axp(1) axp(1)]+0.03, [axp(2) axp(2)+axp(4)+0.07], 'linewidth', 1);
% annotation('textbox', [0.92 0.02 0.1 0.1],'EdgeColor','none','String','\textbf{$x$}','FitBoxToText','on','fontsize',26,'Interpreter','latex');
% annotation('textbox', [0.07 0.92 0.1 0.1],'EdgeColor','none','String','\textbf{$y$}','FitBoxToText','on','fontsize',26,'Interpreter','latex');
% annotation('textbox', [0.07 0.02 0.1 0.1],'EdgeColor','none','String','\textbf{$O$}','FitBoxToText','on','fontsize',26,'Interpreter','latex');
% set(gca,'fontsize', 26)

%% Circular Plate
% view(2)
% xlim([-1.2, 1.2]); ylim([-1.2, 1.2]);
% set(gca, 'YTick', [-0.88, 0.88], 'XTick', [-0.88, 0.88], 'TickLength', [0.001 0.001], 'TickLabelInterpreter', 'latex')
% pgca = gca;
% pgca.YTickLabel = ["$-R$", "$R$"];
% pgca.XTickLabel = ["$-R$", "$R$"];
% pgca.XRuler.TickLabelGapOffset = 4;
% pgca.YRuler.TickLabelGapOffset = 8;
% set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
% set(gca, 'box', 'off', 'GridLineStyle', 'none')
% axp = get(gca,'Position');
% annotation('arrow', [axp(1) axp(1)+axp(3)+0.03], [axp(2)+axp(4)/2 axp(2)+axp(4)/2], 'linewidth', 1);
% annotation('arrow', [axp(1)+axp(3)/2 axp(1)+axp(3)/2], [axp(2) axp(2)+axp(4)+0.03], 'linewidth', 1);
% annotation('textbox', [0.88 0.41 0.1 0.1],'EdgeColor','none','String','\textbf{$x$}','FitBoxToText','on','fontsize',24,'Interpreter','latex');
% annotation('textbox', [0.46 0.85 0.1 0.1],'EdgeColor','none','String','\textbf{$y$}','FitBoxToText','on','fontsize',24,'Interpreter','latex');
% annotation('textbox', [0.46 0.41 0.1 0.1],'EdgeColor','none','String','\textbf{$O$}','FitBoxToText','on','fontsize',24,'Interpreter','latex');
% set(gca,'fontsize', 24)

% element edges
create_el_edges(p,q,U,V,CP)

axis equal;

hold off;

end