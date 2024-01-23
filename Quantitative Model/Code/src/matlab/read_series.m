if strcmp(series_id(1:3),'sim') == 1
    n_periods = n_sim_periods;
    n_series  = n_sims;
elseif strcmp(series_id(1:3),'sam') == 1
    n_periods = n_sample_periods; 
    n_series  = 1;
else 
    n_periods  = n_irf_periods;
    n_series  = 1;
end

fid = fopen([data_path, strcat(series_id, '_state_series.txt')]);
% data_text       = textscan(fid, repmat('%q',1,n_periods*smolyak_d*n_series),'Headerlines',0);
data_text = textscan(fid,'%f');
data_text = data_text{1};
state_series      = reshape(data_text,[smolyak_d,n_periods,n_series]);
fclose(fid);

fid = fopen([data_path, strcat(series_id, '_vars_series.txt')]);
% data_text       = textscan(fid, repmat('%q',1,n_periods*n_interp*n_series),'Headerlines',0);
% data_text       = [data_text{:}];
data_text = textscan(fid,'%f');
data_text = data_text{1};
other_vars_series      = reshape(data_text,[n_interp,n_periods,n_series]);
fclose(fid);

fid = fopen([data_path, strcat(series_id, '_shock_series.txt')]);
% data_text       = textscan(fid, repmat('%q',1,n_periods*n_shocks*n_series),'Headerlines',0);
% data_text       = [data_text{:}];
data_text = textscan(fid,'%f');
data_text = data_text{1};
shock_series      = reshape(data_text,[n_shocks,n_periods,n_series]);
fclose(fid);


if strcmp(series_id(1:3),'sam') == 1 || strcmp(series_id(1:2),'bg') == 1 || ...
        strcmp(series_id(1:3),'fix') == 1 || strcmp(series_id(1:3),'fx0') == 1 
    eval(strcat(series_id, '_series(:,:) = [shock_series; state_series; other_vars_series];'));
else
    fid = fopen([data_path, strcat(series_id, '_value_series.txt')]);
    data_text       = textscan(fid, repmat('%q',1,n_periods*14*n_series),'Headerlines',0);
    data_text       = [data_text{:}];
    value_series      = reshape(str2double(data_text),[14, n_periods, n_series]);
    fclose(fid);
    if strcmp(series_id(1:3),'sim') == 1
        eval(strcat(series_id, '_series(:,:,:) = [shock_series; state_series; other_vars_series; value_series];'));
    else 
        eval(strcat(series_id, '_series(:,:) = [shock_series; state_series; other_vars_series; value_series];'));    
    end
end
