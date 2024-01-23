% -------------------------------------------------------------------------
% create_safety_fig
% -------------------------------------------------------------------------
% authors:         Rohan Kekre and Moritz Lenel
% for updates see: https://github.com/KekreLenel/safety
% -------------------------------------------------------------------------

% disaster shock data
    
    % disaster shock data
    data_tmp     = csvread('../src/data/fig_10.csv',1,0);
    T_data       = size(data_tmp, 1);

    omg_data        = data_tmp(:,1);
    omg_data_lc     = data_tmp(:,7);
    omg_data_uc     = data_tmp(:,13);
    
    log_q_data     = data_tmp(:,2);
    log_q_data_lc  = data_tmp(:,8);
    log_q_data_uc  = data_tmp(:,14);
    
    diff_i_data     = data_tmp(:,3);
    diff_i_data_lc  = data_tmp(:,9);
    diff_i_data_uc  = data_tmp(:,15);
    
    excA_data       = data_tmp(:,4);
    excA_data_lc    = data_tmp(:,10);
    excA_data_uc    = data_tmp(:,16);
    
    yh_data         = data_tmp(:,5);
    yh_data_lc      = data_tmp(:,11);
    yh_data_uc      = data_tmp(:,17);
    
    y_diff_data     = data_tmp(:,6);
    y_diff_data_lc  = data_tmp(:,12);
    y_diff_data_uc  = data_tmp(:,18);    
    
    Interpreter = 'latex';
    fontsize = 14;
    % Define colors
    blue        = [0       0.4470   0.7410];
    lightblue   = [107,174,214]./255;
    darkblue    = [8,69,148]./255;

    tmp_irf = collected_irfs(:,:,2,1);
    T_irf  = size(tmp_irf,1);
    fig_handle = figure;
    set(gcf,'Visible', 'off');
    subplot(2,3,1)
    plot((0:(T_data-1)),zeros(1,T_data),'-k');
    hold on 
    yplot = [omg_data_lc',omg_data_uc(end:-1:1)'];
    fill([0:(T_data-1),(T_data-1):-1:0],yplot,1,'facecolor', lightblue, 'edgecolor', 'none','facealpha',0.5);
    plot((0:(T_data-1)),omg_data,'LineWidth',2,'Color',blue);
    hold on
    irf_vec = collected_irfs(:,irf_idxs.omg,2,1)/omg_scaling;
    irf_scaling = omg_data(1)/irf_vec(1);
    plot(0:(T_irf-1),irf_scaling*irf_vec,'LineWidth',2,'Color',darkblue);
    xlim([0,11]);
    ylabel('bp','Fontsize',fontsize,'Interpreter',Interpreter);
    hold off
    ax2 = gca;
    ax2.FontSize = fontsize;
    ax2.TickLabelInterpreter = Interpreter;
    title({'3-mo swapped G10-Tbill','$\omega$ (scaled)'},'Fontsize',fontsize,'Interpreter',Interpreter);
    
    subplot(2,3,2)
    plot((0:(T_data-1)),zeros(1,T_data),'-k');
    hold on 
    yplot = [log_q_data_lc',log_q_data_uc(end:-1:1)'];
    fill([0:(T_data-1),(T_data-1):-1:0],yplot,1,'facecolor', lightblue, 'edgecolor', 'none','facealpha',0.5);
    plot((0:(T_data-1)),log_q_data,'LineWidth',2,'Color',blue);
    hold on
    irf_vec = collected_irfs(:,irf_idxs.log_qx,2,1);
    plot(0:(T_irf-1),irf_scaling*irf_vec,'LineWidth',2,'Color',darkblue);
    xlim([0,11]);
    ylabel('bp','Fontsize',fontsize,'Interpreter',Interpreter);
    hold off
    ax2 = gca;
    ax2.FontSize = fontsize;
    ax2.TickLabelInterpreter = Interpreter;
    title({'Real exchange rate',irf_titles{irf_idxs.log_qx}},'Fontsize',fontsize,'Interpreter',Interpreter);
    
    subplot(2,3,3)
    plot((0:(T_data-1)),zeros(1,T_data),'-k');
    hold on 
    yplot = [diff_i_data_lc',diff_i_data_uc(end:-1:1)'];
    fill([0:(T_data-1),(T_data-1):-1:0],yplot,1,'facecolor', lightblue, 'edgecolor', 'none','facealpha',0.5);
    h{1} =plot((0:(T_data-1)),diff_i_data,'LineWidth',2,'Color',blue);
    hold on
    irf_vec = 4*(collected_irfs(:,irf_idxs.if,2,1) - collected_irfs(:,irf_idxs.ih,2,1));
    h{2} = plot(0:(T_irf-1),irf_scaling*irf_vec,'LineWidth',2,'Color',darkblue);
    xlim([0,11]);
    ylabel('bp','Fontsize',fontsize,'Interpreter',Interpreter);
    hold off
    legend([h{:}],{'Data','Model'},'Fontsize',fontsize-4,'Interpreter',Interpreter,'Location','ne')
    ax2 = gca;
    ax2.FontSize = fontsize;
    ax2.TickLabelInterpreter = Interpreter;    
    title({'3-mo G10-Tbill','$4(i^\ast-i$)'},'Fontsize',fontsize,'Interpreter',Interpreter); 
    
    subplot(2,3,4)
    plot((0:(T_data-1)),zeros(1,T_data),'-k');
    hold on 
    yplot = [excA_data_lc',excA_data_uc(end:-1:1)'];
    fill([0:(T_data-1),(T_data-1):-1:0],yplot,1,'facecolor', lightblue, 'edgecolor', 'none','facealpha',0.5);
    plot((0:(T_data-1)),excA_data,'LineWidth',2,'Color',blue);
    hold on
    irf_vec = collected_irfs(:,irf_idxs.excA,2,1);
    plot(0:(T_irf-1),irf_scaling*irf_vec,'LineWidth',2,'Color',darkblue);
    xlim([0,11]);
    ylabel('bp','Fontsize',fontsize,'Interpreter',Interpreter);
    hold off
    ax2 = gca;
    ax2.FontSize = fontsize;
    ax2.TickLabelInterpreter = Interpreter;    
    title({'3-mo real excess MSCI ACWI',irf_titles{irf_idxs.excA}},'Fontsize',fontsize,'Interpreter',Interpreter);  
    
    subplot(2,3,5)
    plot((0:(T_data-1)),zeros(1,T_data),'-k');
    hold on 
    yplot = [yh_data_lc',yh_data_uc(end:-1:1)'];
    fill([0:(T_data-1),(T_data-1):-1:0],yplot,1,'facecolor', lightblue, 'edgecolor', 'none','facealpha',0.5);
    plot((0:(T_data-1)),yh_data,'LineWidth',2,'Color',blue);
    hold on
    irf_vec = collected_irfs(:,irf_idxs.log_yh,2,1);
    plot(0:(T_irf-1),irf_scaling*irf_vec,'LineWidth',2,'Color',darkblue);
    xlim([0,11]);
    ylabel('bp','Fontsize',fontsize,'Interpreter',Interpreter);
    hold off
    ax2 = gca;
    ax2.FontSize = fontsize;
    ax2.TickLabelInterpreter = Interpreter;    
    title({'U.S. IP',' $\log y $'},'Fontsize',fontsize,'Interpreter',Interpreter);  
    
    
    subplot(2,3,6)
    plot((0:(T_data-1)),zeros(1,T_data),'-k');
    hold on 
    yplot = [y_diff_data_lc',y_diff_data_uc(end:-1:1)'];
    fill([0:(T_data-1),(T_data-1):-1:0],yplot,1,'facecolor', lightblue, 'edgecolor', 'none','facealpha',0.5);
    plot((0:(T_data-1)),y_diff_data,'LineWidth',2,'Color',blue);
    hold on
    irf_vec = collected_irfs(:,irf_idxs.log_yf,2,1) - collected_irfs(:,irf_idxs.log_yh,2,1);
    plot(0:(T_irf-1),irf_scaling*irf_vec,'LineWidth',2,'Color',darkblue);
    xlim([0,11]);
    ylabel('bp','Fontsize',fontsize,'Interpreter',Interpreter);
    hold off
    ax2 = gca;
    ax2.FontSize = fontsize;
    ax2.TickLabelInterpreter = Interpreter;    
    title({'G10-U.S. IP',' $\log y^\ast -\log y $'},'Fontsize',fontsize,'Interpreter',Interpreter); 
    
    set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0.0 0.0 8.5 5.3],'PaperSize',[8.5 5.3]);
    print(fig_handle,[fig_path, '/fig_10'],'-dpdf', '-r600')

    
