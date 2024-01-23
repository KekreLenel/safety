% -------------------------------------------------------------------------
% make_figure.m: creates formatted matlab figures 
% -------------------------------------------------------------------------
% authors:         Rohan Kekre and Moritz Lenel
% for updates see: https://github.com/KekreLenel/safety
% -------------------------------------------------------------------------
function  fig_handle = make_figure( fig_size, series_set, series_set_2, x_vec_array, y_vec_array, title_array, ylabel_array, legend_array, legend_loc1, legend_loc2, datetick_dmy, thickness, colorscheme)

    fontsize = 14;
    Interpreter = 'latex';

    n_series  = size(series_set,2);
    n_series2 = size(series_set_2,2);

    if exist('thickness','var') && ~isempty(thickness)
        line_thick = thickness;
    else
        line_thick = 2;
    end

    blue        = [0       0.4470   0.7410];
    lightblue   = [107,174,214]./255;
    darkblue    = [8,69,148]./255;
    lightred    = [0.9500  0.4250   0.2];
    darkred     = [0.6350  0.0780   0.1840];
    red         = [0.8500  0.3250   0.0980];
    yellow      = [0.70    0.640    0.1250];
    green       = [0.1     0.75     0.2];
    grey        = [0.4     0.4      0.4];


    if n_series == 1
        color_list = {darkblue};
        style_list = {'-'};
    elseif n_series == 2
        color_list = {lightblue, darkblue};
        if exist('colorscheme','var')
            if strcmp(colorscheme,'modelvsdata') == 1
                color_list = {darkblue, lightblue};
            elseif strcmp(colorscheme,'colorswitch') == 1
                color_list = {darkblue, lightblue};
            end
        end
        style_list = {'-', '-'};
    elseif n_series == 4
        color_list = {lightblue, darkblue, red, green};
        style_list = {'-', '-','-', '-'};
        if exist('colorscheme','var')
            if strcmp(colorscheme,'modelvsdata') == 1
                color_list = {darkblue, red, green, lightblue};
            end
        end
    else
        color_list = {lightblue, blue, darkblue};
        style_list = {'-', '-','-'};
    end

    n_subplots = fig_size(1)*fig_size(2);

    fig_handle = figure;
    set(gcf,'Visible', 'off');
    for nnn = 1:n_subplots 

        subplot(fig_size(1), fig_size(2), nnn)

        x_vec = x_vec_array{nnn};
        y_vec = y_vec_array{nnn};

        length = size(x_vec,1);

        plot(x_vec, zeros(1,length),'-k');
        hold on
        for ppp = 1:n_series
            if ~any(isnan( y_vec(:,series_set{ppp})) )    
                h{ppp} = plot(x_vec(:,series_set{ppp}), y_vec(:,series_set{ppp}) ,'Linewidth',line_thick, 'Color',color_list{ppp},'Linestyle',style_list{ppp});
            end
        end

        title(title_array{nnn},'Fontsize',fontsize,'Interpreter',Interpreter)

        ylabel(ylabel_array{nnn},'Fontsize',fontsize,'Interpreter',Interpreter)

        if exist('legend_array','var')
            if ~isempty(legend_array) && nnn == legend_loc1
            legend([h{:}],legend_array,'Fontsize',fontsize-4,'Interpreter',Interpreter,'Location',legend_loc2)
            end
        end
        if exist('datetick_dmy','var')
            if (datetick_dmy == 1)
            datetick;
            end
        end

        xlim([x_vec(1,1),x_vec(end,1)]);
        if (min(min(y_vec(:,[series_set_2{:}]))) - max(max(y_vec(:,[series_set_2{:}]))) ==0)
            ylim([-0.01 + max(max(y_vec(:,[series_set_2{:}]))), max(max(y_vec(:,[series_set_2{:}]))) + 0.01]);    
        elseif isnan(min(min(y_vec(:,[series_set_2{:}]))))
            ylim([-0.01 , 0.01])
        else
            tmp_spacer = 0.05*abs(max(max(y_vec(:,[series_set_2{:}]))) - min(min(y_vec(:,[series_set_2{:}]))));
            ylim([min(min(y_vec(:,[series_set_2{:}]))) - tmp_spacer, max(max(y_vec(:,[series_set_2{:}])))+ tmp_spacer]);
        end

        ax2 = gca;
        ax2.FontSize = fontsize;
        ax2.TickLabelInterpreter = Interpreter;
    end


end


