function write_FIS_diagnostic_plots(plotfolder,plot_num,par,hplot,Rplot,T_surface_plot,To_plot, ...
    S,dhdt_cond,dhdt_odiff,S_init,div_hv,kappa_del2_h,v_np1,nan_mask)
% Write the standard FIS main-output and source diagnostic plots.

jb = 2;
je = par.nj-1;
plot_idx = (je:-1:jb).';
lat = 90-par.theta(plot_idx);

hplot = hplot(:);
Rplot = Rplot(:);
T_surface_plot = T_surface_plot(:);
To_plot = To_plot(:);
S = S(:);
dhdt_cond = dhdt_cond(:);
dhdt_odiff = dhdt_odiff(:);
S_init = S_init(:);
div_hv = div_hv(:);
kappa_del2_h = kappa_del2_h(:);
v_np1 = v_np1(:);
nan_mask = nan_mask(:);

hplot(hplot==0) = NaN;
hplot2 = Rplot*par.Hcr;

fig = figure('Visible','off');
subplot(4,1,1)
[ax,hl1,hl2] = plotyy(lat,hplot(plot_idx).*nan_mask(plot_idx), ...
    lat,-v_np1(plot_idx).*nan_mask(plot_idx)*par.year);
set([hl1,hl2],'linewidth',2);
valid_v = abs(v_np1(~isnan(v_np1)));
if isempty(valid_v)
    vmax = 0;
else
    vmax = par.year*max(valid_v);
end
h2b = ylabel(ax(1),'h (m)');
h2a = ylabel(ax(2),'v (m/yr)');
h3 = title(sprintf('hcr = %.2d,;v max=%3.0fm/yr',par.Hcr,vmax));
set(ax(2),'xticklabel','');
xlim(ax(1),[90-par.theta(je),90-par.theta(jb)]);
xlim(ax(2),[90-par.theta(je),90-par.theta(jb)]);
set([gca,h2a,h2b,h3],'fontsize',10);
set_axis_limits(ax(1),lat,hplot(plot_idx));
set_axis_limits(ax(2),lat,v_np1(plot_idx)*par.year);
set(gca,'box','off');
set(ax(2),'XAxisLocation','top','linewidth',1,'XTickLabel',[])

subplot(4,1,2)
hl1 = plot(lat,T_surface_plot(plot_idx)-par.T_f,lat,To_plot(plot_idx)-par.T_f,':');
set(hl1,'linewidth',1.5);
grid on
xlim([90-par.theta(je),90-par.theta(jb)]);
legend({'Ts','To'});
h2 = ylabel('(C)');
h3 = title('Surface Temperature');
set([gca,h2,h3],'fontsize',10);

subplot(4,1,3)
hl1 = plot(lat,S(plot_idx).*par.s_col(plot_idx).*nan_mask(plot_idx)*par.year,'r');
hold on
hl2 = plot(lat,div_hv(plot_idx).*par.s_col(plot_idx).*nan_mask(plot_idx)*par.year,'--g');
hl3 = plot(lat,kappa_del2_h(plot_idx).*par.s_col(plot_idx).*nan_mask(plot_idx)*par.year,'--c');
rhs = div_hv-kappa_del2_h;
hl4 = plot(lat,rhs(plot_idx).*par.s_col(plot_idx).*nan_mask(plot_idx)*par.year,'--b');
xlim([90-par.theta(je),90-par.theta(jb)]);
h1 = xlabel('latitude');
h2 = ylabel('m/yr');
h3 = title(sprintf('terms in h eqn'));
h4 = legend('S','\nabla vh','k\nabla^2h','rhs','Location','southeast');
set([hl1,hl2,hl3,hl4],'linewidth',2)
set([gca,h1,h2,h3],'fontsize',10);
set([h4],'fontsize',8);

subplot(4,1,4)
[ax,hl1,hl2] = plotyy(lat,Rplot(plot_idx),lat,hplot2(plot_idx));
hl1.LineStyle = ':';
hl2.LineStyle = '--';
set([hl1,hl2],'linewidth',2);
h1 = xlabel(ax(1),'latitude');
h2b = ylabel(ax(1),'R');
h2a = ylabel(ax(2),'h (m)');
h3 = title('Fractional Ice Cover (R) and height if uniformly spread over cell');
set([gca,h1,h2a,h2b,h3],'fontsize',10);

set(fig, 'PaperUnits', 'inches');
set(fig, 'PaperSize', [5.5 6]);
set(fig, 'PaperPosition', [0 0 5.5 6]);
saveas(fig,sprintf('%s/main-output-%.2d.pdf',plotfolder,plot_num));
close(fig);

fig = figure('Visible','off');
plot(lat,S(plot_idx)*par.year*100,'-+k',lat,dhdt_cond(plot_idx)*par.year*100,'--r', ...
    lat,dhdt_odiff(plot_idx)*par.year*100,'-+b',lat,S_init(plot_idx)*par.year*100,'--g');
yl = ylabel('Source Term (cm/yr)');
xl = xlabel('Latitude');
tl = title('Source Terms');
legend({'Stot','C','Odiff','S-EBM'});
set([gca,xl,yl,tl],'fontsize',10);

saveas(fig,sprintf('%s/source-%.2d.png',plotfolder,plot_num));
close(fig);

end
