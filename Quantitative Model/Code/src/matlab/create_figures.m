% -------------------------------------------------------------------------
% create_figures.m: create figures
% -------------------------------------------------------------------------
% authors:         Rohan Kekre and Moritz Lenel
% for updates see: https://github.com/KekreLenel/safety
% -------------------------------------------------------------------------

% length of figures
start_t = 1;
len     = 29;

% DISASTER PROBABILITY SHOCK FIGURES
% ---------------------------------------------------------------------------------

irf_idx = 1;

series_to_plot = [irf_idxs.p; irf_idxs.excA; irf_idxs.excrf; irf_idxs.log_qx; irf_idxs.thth; irf_idxs.log_yh];
y_array = {}; x_array = {}; title_array = {}; ylabel_array = {};
for iii = 1:size(series_to_plot,1)
    y_array{iii}      = squeeze(collected_irfs(start_t:len,series_to_plot(iii),irf_idx,:) );
    x_array{iii}      = repmat((1:len-start_t+1)', [1, n_comp]);
    title_array{iii}  = irf_titles{series_to_plot(iii)};
    ylabel_array{iii} = "bp";
end

try
    for ccc = 1:n_comp
        legend_array = {}; legend_loc1 = 1; legend_loc2 = "ne"; 
        fig_handle = make_figure([2,3],{ccc},{ccc},x_array,y_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2);
        set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0.0 0.0 8.5 5.0],'PaperSize',[8.5 5.0]);
        if ccc == ix_bm
        print(fig_handle,[fig_path, '/fig_2'],'-dpdf', '-r600')
        end
        print(fig_handle,[fig_path, '/addl_figs/p_', num2str(ccc), '_fig'],'-dpdf', '-r600')
    end
catch me 
end

            
% appendix figure 1:
series_to_plot = [irf_idxs.p;  irf_idxs.Erh; irf_idxs.Erf;  irf_idxs.excrk; irf_idxs.log_qk;       irf_idxs.excrf; ...
                  irf_idxs.log_qx; irf_idxs.log_E;   irf_idxs.thth; irf_idxs.log_kh;    irf_idxs.nfa_ybar; irf_idxs.val_ybar];
y_array = {}; x_array = {}; title_array = {}; ylabel_array = {};
for iii = 1:size(series_to_plot,1)
    y_array{iii}      = squeeze(collected_irfs(start_t:len,series_to_plot(iii),irf_idx,:) );
    x_array{iii}      = repmat((1:len-start_t+1)', [1, n_comp]);
    title_array{iii}  = irf_titles{series_to_plot(iii)};
    ylabel_array{iii} = "bp";
end

try
    for ccc = 1:n_comp
        legend_array = {}; legend_loc1 = 1; legend_loc2 = "ne"; 
        fig_handle = make_figure([3,4],{ccc},{ccc},x_array,y_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2);
        set(gcf,'PaperType','usletter','PaperOrientation','landscape','PaperPosition',[0.0 0.0 12.5 8.5],'PaperSize',[12.5 8.5]);
        if ccc == ix_bm
        print(fig_handle,[fig_path, '/fig_11'],'-dpdf', '-r600')
        end
        print(fig_handle,[fig_path, '/addl_figs/p_', num2str(ccc), '_1_fig'],'-dpdf', '-r600')
    end
catch me 
end

% appendix figure 2:
series_to_plot = [irf_idxs.nx_ybar;  irf_idxs.log_ch; irf_idxs.log_cf;  irf_idxs.log_x; irf_idxs.log_infl_h;       irf_idxs.log_wh; ...
                  irf_idxs.log_lh; irf_idxs.log_yh;   irf_idxs.log_infl_f; irf_idxs.log_wf;    irf_idxs.log_lf; irf_idxs.log_yf];
y_array = {}; x_array = {}; title_array = {}; ylabel_array = {};
for iii = 1:size(series_to_plot,1)
    y_array{iii}      = squeeze(collected_irfs(start_t:len,series_to_plot(iii),irf_idx,:) );
    x_array{iii}      = repmat((1:len-start_t+1)', [1, n_comp]);
    title_array{iii}  = irf_titles{series_to_plot(iii)};
    ylabel_array{iii} = "bp";
end

try
    for ccc = 1:n_comp
        legend_array = {}; legend_loc1 = 1; legend_loc2 = "ne"; 
        fig_handle = make_figure([3,4],{ccc},{ccc},x_array,y_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2);
        set(gcf,'PaperType','usletter','PaperOrientation','landscape','PaperPosition',[0.0 0.0 12.5 8.5],'PaperSize',[12.5 8.5]);
        if ccc == ix_bm
        print(fig_handle,[fig_path, '/fig_12'],'-dpdf', '-r600')
        end
        print(fig_handle,[fig_path, '/addl_figs/p_', num2str(ccc), '_2_fig'],'-dpdf', '-r600')
    end
catch me 
end


% SAFETY SHOCK FIGURES
% ---------------------------------------------------------------------------------

irf_idx = 2;

series_to_plot = [irf_idxs.omg; irf_idxs.excA; irf_idxs.excrf; irf_idxs.log_qx; irf_idxs.thth; irf_idxs.log_yh];
y_array = {}; x_array = {}; title_array = {}; ylabel_array = {};
for iii = 1:size(series_to_plot,1)
    y_array{iii}      = squeeze(collected_irfs(start_t:len,series_to_plot(iii),irf_idx,:) );
    x_array{iii}      = repmat((1:len-start_t+1)', [1, n_comp]);
    title_array{iii}  = irf_titles{series_to_plot(iii)};
    ylabel_array{iii} = "bp";
end
 
try
    for ccc = 1:n_comp
        legend_array = {}; legend_loc1 = 1; legend_loc2 = "ne"; 
        fig_handle = make_figure([2,3],{ccc},{ccc},x_array,y_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2);
        set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0.0 0.0 8.5 5.0],'PaperSize',[8.5 5.0]);
        print(fig_handle,[fig_path, '/addl_figs/omg_', num2str(ccc), '_fig'],'-dpdf', '-r600')
    end
catch me 
end

% compare benchmark, symmetry, symmetry no stickiness
% try
ix_vec = {ix_symm_flex, ix_symm, ix_bm};
legend_array = {['$\gamma=\gamma^\ast$,' newline '$\chi^{W}=0$'], '$\gamma=\gamma^\ast$', 'Model'}; legend_loc1 = 6; legend_loc2 = "se"; 

fig_handle = make_figure([2,3],ix_vec,ix_vec,x_array,y_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2);
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0.0 0.0 8.5 5.0],'PaperSize',[8.5 5.0]);
print(fig_handle,[fig_path, '/addl_figs/omg_comp3_fig'],'-dpdf', '-r600')
print(fig_handle,[fig_path, '/fig_3'],'-dpdf', '-r600')

fig_handle = make_figure([2,3],ix_vec(1:2),ix_vec,x_array,y_array, title_array, ylabel_array, legend_array(1:2), legend_loc1, legend_loc2);
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0.0 0.0 8.5 5.0],'PaperSize',[8.5 5.0]);
print(fig_handle,[fig_path, '/addl_figs/omg_comp2_fig'],'-dpdf', '-r600')

fig_handle = make_figure([2,3],ix_vec(1),ix_vec,x_array,y_array, title_array, ylabel_array, legend_array(1), legend_loc1, legend_loc2);
set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0.0 0.0 8.5 5.0],'PaperSize',[8.5 5.0]);
print(fig_handle,[fig_path, '/addl_figs/omg_comp1_fig'],'-dpdf', '-r600')
% catch me 
% end

% compare benchmark to alternative taylor rules

try

    series_to_plot = [irf_idxs.ih; irf_idxs.excA; irf_idxs.excrf; irf_idxs.log_qx; irf_idxs.thth; irf_idxs.log_yh];
    y_array = {}; x_array = {}; title_array = {}; ylabel_array = {};
    for iii = 1:size(series_to_plot,1)
        y_array{iii}      = squeeze(collected_irfs(start_t:len,series_to_plot(iii),irf_idx,:) );
        x_array{iii}      = repmat((1:len-start_t+1)', [1, n_comp]);
        title_array{iii}  = irf_titles{series_to_plot(iii)};
        ylabel_array{iii} = "bp";
    end
     
    ix_vec = {ix_tayl_rho, ix_tayl_y0, ix_bm, ix_symm, ix_symm_flex};
    legend_array = {['$\phi^y=0.5/4$,' newline '$\rho^i = 0.5$'], '$\phi^y=0.5/4$', 'Model'}; legend_loc1 = 6; legend_loc2 = "se"; 
    
    fig_handle = make_figure([2,3],ix_vec(1:3),ix_vec,x_array,y_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2);
    set(gcf,'PaperType','usletter','PaperOrientation','portrait','PaperPosition',[0.0 0.0 8.5 5.0],'PaperSize',[8.5 5.0]);
    print(fig_handle,[fig_path, '/fig_5'],'-dpdf', '-r600')

catch me 
end

% appendix figure 1:
series_to_plot = [irf_idxs.omg;  irf_idxs.Erh; irf_idxs.Erf;  irf_idxs.excrk; irf_idxs.log_qk;       irf_idxs.excrf; ...
                  irf_idxs.log_qx; irf_idxs.log_E;   irf_idxs.thth; irf_idxs.log_kh;    irf_idxs.nfa_ybar; irf_idxs.val_ybar];
y_array = {}; x_array = {}; title_array = {}; ylabel_array = {};
for iii = 1:size(series_to_plot,1)
    y_array{iii}      = squeeze(collected_irfs(start_t:len,series_to_plot(iii),irf_idx,:) );
    x_array{iii}      = repmat((1:len-start_t+1)', [1, n_comp]);
    title_array{iii}  = irf_titles{series_to_plot(iii)};
    ylabel_array{iii} = "bp";
end

try
    for ccc = 1:n_comp
        legend_array = {}; legend_loc1 = 1; legend_loc2 = "ne"; 
        fig_handle = make_figure([3,4],{ccc},{ccc},x_array,y_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2);
        set(gcf,'PaperType','usletter','PaperOrientation','landscape','PaperPosition',[0.0 0.0 12.5 8.5],'PaperSize',[12.5 8.5]);
        if ccc == ix_bm
        print(fig_handle,[fig_path, '/fig_19'],'-dpdf', '-r600')
        end
        print(fig_handle,[fig_path, '/addl_figs/omg_', num2str(ccc), '_1_fig'],'-dpdf', '-r600')
    end
catch me 
end

    
% appendix figure 2:
series_to_plot = [irf_idxs.nx_ybar;  irf_idxs.log_ch; irf_idxs.log_cf;  irf_idxs.log_x; irf_idxs.log_infl_h;       irf_idxs.log_wh; ...
                  irf_idxs.log_lh; irf_idxs.log_yh;   irf_idxs.log_infl_f; irf_idxs.log_wf;    irf_idxs.log_lf; irf_idxs.log_yf];
y_array = {}; x_array = {}; title_array = {}; ylabel_array = {};
for iii = 1:size(series_to_plot,1)
    y_array{iii}      = squeeze(collected_irfs(start_t:len,series_to_plot(iii),irf_idx,:) );
    x_array{iii}      = repmat((1:len-start_t+1)', [1, n_comp]);
    title_array{iii}  = irf_titles{series_to_plot(iii)};
    ylabel_array{iii} = "bp";
end

try
    for ccc = 1:n_comp
        legend_array = {}; legend_loc1 = 1; legend_loc2 = "ne"; 
        fig_handle = make_figure([3,4],{ccc},{ccc},x_array,y_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2);
        set(gcf,'PaperType','usletter','PaperOrientation','landscape','PaperPosition',[0.0 0.0 12.5 8.5],'PaperSize',[12.5 8.5]);
        if ccc == ix_bm
        print(fig_handle,[fig_path, '/fig_20'],'-dpdf', '-r600')
        end
        print(fig_handle,[fig_path, '/addl_figs/omg_', num2str(ccc), '_2_fig'],'-dpdf', '-r600')
    end
catch me 
end



% Z SHOCK FIGURES
% ---------------------------------------------------------------------------------

irf_idx = 3;

% appendix figure 1:
series_to_plot = [irf_idxs.log_z;  irf_idxs.Erh; irf_idxs.Erf;  irf_idxs.excrk; irf_idxs.log_qk;       irf_idxs.excrf; ...
                  irf_idxs.log_qx; irf_idxs.log_E;   irf_idxs.thth; irf_idxs.log_kh;    irf_idxs.nfa_ybar; irf_idxs.val_ybar];
y_array = {}; x_array = {}; title_array = {}; ylabel_array = {};
for iii = 1:size(series_to_plot,1)
    y_array{iii}      = squeeze(collected_irfs(start_t:len,series_to_plot(iii),irf_idx,:) );
    x_array{iii}      = repmat((1:len-start_t+1)', [1, n_comp]);
    title_array{iii}  = irf_titles{series_to_plot(iii)};
    ylabel_array{iii} = "bp";
end

try
    for ccc = 1:n_comp
        legend_array = {}; legend_loc1 = 1; legend_loc2 = "ne"; 
        fig_handle = make_figure([3,4],{ccc},{ccc},x_array,y_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2);
        set(gcf,'PaperType','usletter','PaperOrientation','landscape','PaperPosition',[0.0 0.0 12.5 8.5],'PaperSize',[12.5 8.5]);
        if ccc == ix_bm
        print(fig_handle,[fig_path, '/fig_15'],'-dpdf', '-r600')
        end
        print(fig_handle,[fig_path, '/addl_figs/z_', num2str(ccc), '_1_fig'],'-dpdf', '-r600')
    end
catch me 
end


    
% appendix figure 2:
series_to_plot = [irf_idxs.nx_ybar;  irf_idxs.log_ch; irf_idxs.log_cf;  irf_idxs.log_x; irf_idxs.log_infl_h;       irf_idxs.log_wh; ...
                  irf_idxs.log_lh; irf_idxs.log_yh;   irf_idxs.log_infl_f; irf_idxs.log_wf;    irf_idxs.log_lf; irf_idxs.log_yf];
y_array = {}; x_array = {}; title_array = {}; ylabel_array = {};
for iii = 1:size(series_to_plot,1)
    y_array{iii}      = squeeze(collected_irfs(start_t:len,series_to_plot(iii),irf_idx,:) );
    x_array{iii}      = repmat((1:len-start_t+1)', [1, n_comp]);
    title_array{iii}  = irf_titles{series_to_plot(iii)};
    ylabel_array{iii} = "bp";
end

try
    for ccc = 1:n_comp
        legend_array = {}; legend_loc1 = 1; legend_loc2 = "ne"; 
        fig_handle = make_figure([3,4],{ccc},{ccc},x_array,y_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2);
        set(gcf,'PaperType','usletter','PaperOrientation','landscape','PaperPosition',[0.0 0.0 12.5 8.5],'PaperSize',[12.5 8.5]);
        if ccc == ix_bm
        print(fig_handle,[fig_path, '/fig_16'],'-dpdf', '-r600')
        end
        print(fig_handle,[fig_path, '/addl_figs/z_', num2str(ccc), '_2_fig'],'-dpdf', '-r600')
    end
catch me 
end

  
 % DISASTER SHOCK FIGURES
% ---------------------------------------------------------------------------------

irf_idx = 4;

% appendix figure 1:
series_to_plot = [irf_idxs.log_z;  irf_idxs.Erh; irf_idxs.Erf;  irf_idxs.excrk; irf_idxs.log_qk;       irf_idxs.excrf; ...
                  irf_idxs.log_qx; irf_idxs.log_E;   irf_idxs.thth; irf_idxs.log_kh;    irf_idxs.nfa_ybar; irf_idxs.val_ybar];
y_array = {}; x_array = {}; title_array = {}; ylabel_array = {};
for iii = 1:size(series_to_plot,1)
    y_array{iii}      = squeeze(collected_irfs(start_t:len,series_to_plot(iii),irf_idx,:) );
    x_array{iii}      = repmat((1:len-start_t+1)', [1, n_comp]);
    title_array{iii}  = irf_titles{series_to_plot(iii)};
    ylabel_array{iii} = "bp";
end

try
for ccc = 1:n_comp
    legend_array = {}; legend_loc1 = 1; legend_loc2 = "ne"; 
    fig_handle = make_figure([3,4],{ccc},{ccc},x_array,y_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2);
    set(gcf,'PaperType','usletter','PaperOrientation','landscape','PaperPosition',[0.0 0.0 12.5 8.5],'PaperSize',[12.5 8.5]);
    if ccc == ix_bm
    print(fig_handle,[fig_path, '/fig_13'],'-dpdf', '-r600')
    end
    print(fig_handle,[fig_path, '/addl_figs/dis_', num2str(ccc), '_1_fig'],'-dpdf', '-r600')
end
catch me 
end


    
% appendix figure 2:
series_to_plot = [irf_idxs.nx_ybar;  irf_idxs.log_ch; irf_idxs.log_cf;  irf_idxs.log_x; irf_idxs.log_infl_h;       irf_idxs.log_wh; ...
                  irf_idxs.log_lh; irf_idxs.log_yh;   irf_idxs.log_infl_f; irf_idxs.log_wf;    irf_idxs.log_lf; irf_idxs.log_yf];
y_array = {}; x_array = {}; title_array = {}; ylabel_array = {};
for iii = 1:size(series_to_plot,1)
    y_array{iii}      = squeeze(collected_irfs(start_t:len,series_to_plot(iii),irf_idx,:) );
    x_array{iii}      = repmat((1:len-start_t+1)', [1, n_comp]);
    title_array{iii}  = irf_titles{series_to_plot(iii)};
    ylabel_array{iii} = "bp";
end

try
    for ccc = 1:n_comp
        legend_array = {}; legend_loc1 = 1; legend_loc2 = "ne"; 
        fig_handle = make_figure([3,4],{ccc},{ccc},x_array,y_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2);
        set(gcf,'PaperType','usletter','PaperOrientation','landscape','PaperPosition',[0.0 0.0 12.5 8.5],'PaperSize',[12.5 8.5]);
        if ccc == ix_bm
        print(fig_handle,[fig_path, '/fig_14'],'-dpdf', '-r600')
        end
        print(fig_handle,[fig_path, '/addl_figs/dis_', num2str(ccc), '_2_fig'],'-dpdf', '-r600')
    end
catch me 
end
        
        
% FOREIGN PRODUCTIVITY SHOCK FIGURES
% --------------------------------------------------------------------------------

irf_idx = 5;

% appendix figure 1:
series_to_plot = [irf_idxs.log_zf;  irf_idxs.Erh; irf_idxs.Erf;  irf_idxs.excrk; irf_idxs.log_qk;       irf_idxs.excrf; ...
                   irf_idxs.log_qx; irf_idxs.log_E;   irf_idxs.thth; irf_idxs.log_kh;    irf_idxs.nfa_ybar; irf_idxs.val_ybar];
y_array = {}; x_array = {}; title_array = {}; ylabel_array = {};
for iii = 1:size(series_to_plot,1)
    y_array{iii}      = squeeze(collected_irfs(start_t:len,series_to_plot(iii),irf_idx,:) );
    x_array{iii}      = repmat((1:len-start_t+1)', [1, n_comp]);
    title_array{iii}  = irf_titles{series_to_plot(iii)};
    ylabel_array{iii} = "bp";
end

try
    for ccc = 1:n_comp
        legend_array = {}; legend_loc1 = 1; legend_loc2 = "ne"; 
        fig_handle = make_figure([3,4],{ccc},{ccc},x_array,y_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2);
        set(gcf,'PaperType','usletter','PaperOrientation','landscape','PaperPosition',[0.0 0.0 12.5 8.5],'PaperSize',[12.5 8.5]);
        if ccc == ix_bm
        print(fig_handle,[fig_path, '/fig_17'],'-dpdf', '-r600')
        end
        print(fig_handle,[fig_path, '/addl_figs/zf_', num2str(ccc), '_1_fig'],'-dpdf', '-r600')
    end
catch me 
end


    
% appendix figure 2:
series_to_plot = [irf_idxs.nx_ybar;  irf_idxs.log_ch; irf_idxs.log_cf;  irf_idxs.log_x; irf_idxs.log_infl_h;       irf_idxs.log_wh; ...
                  irf_idxs.log_lh; irf_idxs.log_yh;   irf_idxs.log_infl_f; irf_idxs.log_wf;    irf_idxs.log_lf; irf_idxs.log_yf];
y_array = {}; x_array = {}; title_array = {}; ylabel_array = {};
for iii = 1:size(series_to_plot,1)
    y_array{iii}      = squeeze(collected_irfs(start_t:len,series_to_plot(iii),irf_idx,:) );
    x_array{iii}      = repmat((1:len-start_t+1)', [1, n_comp]);
    title_array{iii}  = irf_titles{series_to_plot(iii)};
    ylabel_array{iii} = "bp";
end

try
    for ccc = 1:n_comp
        legend_array = {}; legend_loc1 = 1; legend_loc2 = "ne"; 
        fig_handle = make_figure([3,4],{ccc},{ccc},x_array,y_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2);
        set(gcf,'PaperType','usletter','PaperOrientation','landscape','PaperPosition',[0.0 0.0 12.5 8.5],'PaperSize',[12.5 8.5]);
        if ccc == ix_bm
        print(fig_handle,[fig_path, '/fig_18'],'-dpdf', '-r600')
        end
        print(fig_handle,[fig_path, '/addl_figs/zf_', num2str(ccc), '_2_fig'],'-dpdf', '-r600')
    end
catch me 
end

