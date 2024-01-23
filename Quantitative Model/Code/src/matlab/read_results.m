%% Read in parameter file

fid = fopen(param_file, 'r');
fmt = repmat('%q ', 1, n_params );
headline = textscan(fid, fmt, 1, 'Delimiter',',','HeaderLines',0);
targ_prm_names = [headline{1:end}];
targ_prm_table = csvread(param_file,1);
n_targ_prm = size(targ_prm_table,2);

trg_prms = cell2struct(cell(1,n_targ_prm),targ_prm_names,2);
for ppp = 1:n_targ_prm
    trg_prms = setfield(trg_prms,targ_prm_names{ppp},targ_prm_table(ppp));
end
trg_prms.debt_to_equity = 1.0;
load('../src/params/additional_params.mat')

%% Read in grid specifications

% read in grid size file
    fid             = fopen([data_path, 'grid.csv']);
    data_text       = textscan(fid, repmat('%q',1,7),'Headerlines',0);
    data_text       = [data_text{:}];
    gridsize_vec    = str2double(data_text);
    fclose(fid);

    fid             = fopen([data_path, 'grid_locs.csv']);
    data_text       = textscan(fid, repmat('%q',1,14),'Headerlines',0);
    data_text       = [data_text{:}];
    gridlocs_vec    = str2double(data_text);
    fclose(fid);
    
    fid             = fopen([data_path, 'extra_data.csv']);
    data_text       = textscan(fid, repmat('%q',1,10),'Headerlines',0);
    data_text       = [data_text{:}];
    extra_data_vec  = str2double(data_text);
    chi0_vec        = extra_data_vec(1:2);
    sig_dis         = extra_data_vec(3);
    sig_omg         = extra_data_vec(4);
    dis_grid_mean   = extra_data_vec(5);
    omg_grid_mean   = extra_data_vec(6);
    yh_ss           = extra_data_vec(7);
    aggr_css        = extra_data_vec(8);


% Assign grid sizes
     n_I           = gridsize_vec(1);
     n_states      = gridsize_vec(2);
     n_active_dims = gridsize_vec(3);
     n_interp      = gridsize_vec(4);
     n_shocks      = gridsize_vec(5);
     n_spread      = gridsize_vec(6);
     smolyak_d     = gridsize_vec(7);
     
     k_grid_mean   = gridlocs_vec(1);
     k_grid_dev    = gridlocs_vec(2);
     tht_h_grid_mean   = gridlocs_vec(3);
     tht_h_grid_dev    = gridlocs_vec(4);
     wh_grid_mean   = gridlocs_vec(5);
     wh_grid_dev    = gridlocs_vec(6);
     wf_grid_mean   = gridlocs_vec(7);
     wf_grid_dev    = gridlocs_vec(8);
     ih_grid_mean   = gridlocs_vec(9);
     ih_grid_dev    = gridlocs_vec(10);
     if_grid_mean   = gridlocs_vec(11);
     if_grid_dev    = gridlocs_vec(12);
     omg_grid_mean   = gridlocs_vec(13);
     omg_grid_dev    = gridlocs_vec(14);


%% Read simulated series and impulse responses

series_id = 'sim';
read_series;

series_id = 'sim_dis';
read_series;

try
series_id = 'sim_omg';
read_series;
catch me
end

try
series_id = 'sim_eps';
read_series;
catch me
end

try
series_id = 'sim_p';
read_series;
catch me
end

try
series_id = 'samp1';
read_series;
series_id = 'samp2';
read_series;
series_id = 'samp3';
read_series;
catch me
end

series_id = 'none_irf';
read_series;

series_id = 'g_irf';
read_series;

series_id = 'p_irf';
read_series;

series_id = 'omg_irf';
read_series;

series_id = 'dis1_irf';
read_series;

series_id = 'zf_irf';
read_series;

try
series_id = 'bg1_irf';
read_series;

series_id = 'fix_irf';
read_series;

series_id = 'fx0_irf';
read_series;

catch me
end




