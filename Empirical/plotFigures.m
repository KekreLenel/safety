%% plotFigures.m
% Plots Figures 1 and 8

runTime = tic;

Interpreter = 'latex';
fontsize = 14;
% Define colors
blue        = [0       0.4470   0.7410];
lightblue   = [107,174,214]./255;
darkblue    = [8,69,148]./255;
lightred    = [0.9500  0.4250   0.2];
darkred     = [0.6350  0.0780   0.1840];
red         = [0.8500  0.3250   0.0980];
yellow      = [0.70    0.640    0.1250];
green       = [0.1     0.75     0.2];
grey        = [0.4     0.4      0.4];

%% Figure 1

datafile = importdata('output/fig_1.csv');
data = datafile.data;
names = datafile.colheaders; 

spreads = figure;
box on;
l1 = plot(data(1:end,strcmp('yearmonth',names)),data(1:end,strcmp('dcpnf3m_y',names)),'LineWidth',2,'Color',darkblue);
hold on;
l2 = plot(data(1:end,strcmp('yearmonth',names)),data(1:end,strcmp('cip_govt_3m_2020',names)),'LineWidth',2,'Color',blue);
set(gca,'ycolor','k');
ylabel('bp','Interpreter',Interpreter);
ylim([-50,400]);
xlim([data(1,strcmp('yearmonth',names)),data(end,strcmp('yearmonth',names))]);
hold on
plot([1997.3334 1998.5834 1999.6666 2000.25 2007.5834 2008.6666 2016.5 2020.1666;1997.3334 1998.5834 1999.6666 2000.25 2007.5834 2008.6666 2016.5 2020.1666], [ylim; ylim; ylim; ylim; ylim; ylim; ylim; ylim]', '--', 'Color', grey)
hold off
text([1997.2], [150], {'Asian', 'crisis'}, 'Interpreter',Interpreter, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'BackgroundColor', 'white')
text([1998.02], [220], {'LTCM', 'fails'}, 'Interpreter',Interpreter, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'BackgroundColor', 'white')
text([1999], [280], {'Y2K', 'risk'}, 'Interpreter',Interpreter, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'BackgroundColor', 'white')
text([1999.9], [175], {'Nasdaq', 'crash'}, 'Interpreter',Interpreter, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'BackgroundColor', 'white')
text([2006.8], [250], {'BNP', 'Paribas', 'freezes', 'funds'}, 'Interpreter',Interpreter, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'BackgroundColor', 'white')
text([2008.2], [390], {'Lehman', 'Brothers', 'fails'}, 'Interpreter',Interpreter, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'BackgroundColor', 'white')
text([2016], [157], {'Brexit', 'vote'}, 'Interpreter',Interpreter, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'BackgroundColor', 'white')
text([2019.8], [250], {'Covid-19', 'global', 'lockdowns'}, 'Interpreter',Interpreter, 'FontWeight', 'bold', 'VerticalAlignment', 'top', 'BackgroundColor', 'white')
leg = legend([l1,l2],{'3m AA-Tbill','3m swapped G10-Tbill'},'Location','southoutside','Orientation','horizontal','Interpreter',Interpreter);
set(leg,'Box','off');
ax2 = gca;
ax2.FontSize = fontsize;
ax2.TickLabelInterpreter = Interpreter;
ax2.YColor = [0 0 0];
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0.0 0.0 8.5 4.0],'PaperSize',[8.5 4.0]);
print(spreads,'output/fig_1','-dpdf', '-r600');

%% Figure 8

datafile = importdata('output/fig_8.csv');
data = datafile.data;
names = datafile.colheaders; 

eF = figure;
plot(data(:,strcmp('yearmonth',names)),zeros(1,size(data,1)),'-k');
hold on;
plot(data(:,strcmp('yearmonth',names)),data(:,strcmp('eF_base',names)),'LineWidth',2,'Color',darkblue);
xlim([data(1,strcmp('yearmonth',names)),data(size(data,1),strcmp('yearmonth',names))]);
ax2 = gca;
ax2.FontSize = fontsize;
ax2.TickLabelInterpreter = Interpreter;
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0.0 0.0 8.5 5.0],'PaperSize',[8.5 5.0]);
print(eF,'output/fig_8','-dpdf', '-r600');

disp(['Total runtime: ',num2str(toc(runTime)),'s']);