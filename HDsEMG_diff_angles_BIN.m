function HDsEMG_diff_angles_BIN

[file_aux, path_aux] = uigetfile('*.otb','Importar Arquivo Exportado pelo OT Biolab', 'MultiSelect', 'on');
tic

signal_par = struct;
configuration_par = struct;
abstract = cell(1, size(file_aux,2));
temp_folder = cell(1, size(file_aux,2));
info_signal = cell(1, size(file_aux,2));
info_patient = cell(1, size(file_aux,2));
info_channel = cell(1, size(file_aux,2));

for i = 1:size(file_aux, 2)
    dir_name = strcat('temp_folder', num2str(i));
    temp_folder{i} = [path_aux strcat(dir_name, '\')];
end

for i = 1:size(file_aux, 2)

    mkdir(temp_folder{i});
    compressed_file = [path_aux file_aux{i}];
    unzip(compressed_file, temp_folder{i});
    
    % Preallocating variables
    child_values = cell(1, 2);
    child_names = cell(1, 2);
    
    % Node abstract
    abstract_nodes = xmlread([temp_folder{i} 'abstract.xml']);
    element_abstract = abstract_nodes.getElementsByTagName('abstract');
    abstract_child_nodes = element_abstract.item(0).getChildNodes;
    count = 1;
    for j = 1:2:(abstract_child_nodes.getLength-2)
        child_values{count} = char(abstract_child_nodes.item(j).getTextContent);
        child_names{count} =  char(abstract_child_nodes.item(j).getNodeName);
        count = count + 1;
    end
    abstract{i} = struct(child_names{1}, child_values{1}, child_names{2}, child_values{2},...
        child_names{3}, child_values{3}, child_names{4}, child_values{4},...
        child_names{5}, child_values{5}, child_names{6}, child_values{6},...
        child_names{7}, child_values{7}, child_names{8}, child_values{8},...
        child_names{9}, child_values{9}, child_names{10}, child_values{10});
    
    % Node signal
    element_signal = abstract_nodes.getElementsByTagName('signal');
    signal_child_nodes = element_signal.item(0).getChildNodes;
    count = 1;
    for k = 1:2:(signal_child_nodes.getLength-2)
        child_values{count} = char(signal_child_nodes.item(k).getTextContent);
        child_names{count} =  char(signal_child_nodes.item(k).getNodeName);
        count = count + 1;
    end
    info_signal{i} = struct(child_names{1}, child_values{1}, child_names{2}, child_values{2},...
        child_names{3}, child_values{3}, child_names{4}, child_values{4},...
        child_names{5}, child_values{5}, child_names{6}, child_values{6},...
        child_names{7}, child_values{7}, child_names{8}, child_values{8},...
        child_names{9}, child_values{9}, child_names{10}, child_values{10},...
        child_names{11}, child_values{11}, child_names{12}, child_values{12},...
        child_names{13}, child_values{13}, child_names{14}, child_values{14},...
        child_names{15}, child_values{15}, child_names{16}, child_values{16},...
        child_names{17}, child_values{17}, child_names{18}, child_values{18},...
        child_names{19}, child_values{19}, child_names{20}, child_values{20});
    
    % Node xml patient
    abstract_2_nodes = xmlread([temp_folder{i} info_signal{i}.abstract_path]);
    element_patient = abstract_2_nodes.getElementsByTagName('patient');
    patient_child_nodes = element_patient.item(0).getChildNodes;
    count = 1;
    for l = 1:2:(patient_child_nodes.getLength-2)
        child_values{count} = char(patient_child_nodes.item(l).getTextContent);
        child_names{count} =  char(patient_child_nodes.item(l).getNodeName);
        count = count + 1;
    end
    info_patient{i} = struct(child_names{1}, child_values{1}, child_names{2}, child_values{2},...
        child_names{3}, child_values{3}, child_names{4}, child_values{4},...
        child_names{5}, child_values{5}, child_names{6}, child_values{6},...
        child_names{7}, child_values{7}, child_names{8}, child_values{8},...
        child_names{9}, child_values{9});
    
    % Node xml channel
    element_channel = abstract_2_nodes.getElementsByTagName('channel');
    channel_child_nodes = element_channel.item(0).getChildNodes;
    count = 1;
    for m = 1:2:(channel_child_nodes.getLength-2)
        child_values{count} = char(channel_child_nodes.item(m).getTextContent);
        child_names{count} =  char(channel_child_nodes.item(m).getNodeName);
        count = count + 1;
    end
    info_channel{i} = struct(child_names{1}, child_values{1}, child_names{2}, child_values{2},...
        child_names{3}, child_values{3}, child_names{4}, child_values{4},...
        child_names{5}, child_values{5}, child_names{6}, child_values{6},...
        child_names{7}, child_values{7}, child_names{8}, child_values{8});
    
    if i == 1
        configuration_par.signal_path{i} = info_signal{i}.signal_path;
        configuration_par.n_channels =  str2double(info_signal{i}.channels);
        configuration_par.fsample =  str2double(info_signal{i}.fsample);
        configuration_par.signal_gain =  str2double(info_signal{i}.signal_gain);
        configuration_par.low_pass_filter =  str2double(info_signal{i}.low_pass_filter);
        configuration_par.high_pass_filter =  str2double(info_signal{i}.high_pass_filter);
        configuration_par.ad_bits =  str2double(info_signal{i}.ad_bits);
        configuration_par.patient_id = info_patient{i}.id;
        configuration_par.emg_side = info_channel{i}.side;
        configuration_par.n_conditions = size(file_aux,2);
        
        configuration_par.angle_stim = str2double(info_signal{i}.comments(1:3));
        if isnan(configuration_par.angle_stim)
            configuration_par.angle_stim = str2double(info_signal{i}.comments(1:2));   
        end

        aux_fsample = configuration_par.fsample;
        aux_signal_gain = configuration_par.signal_gain; 
    else
        configuration_par.signal_path{i} = info_signal{i}.signal_path;
        configuration_par.fsample = [configuration_par.fsample str2double(info_signal{i}.fsample)];
        configuration_par.signal_gain = [configuration_par.signal_gain str2double(info_signal{i}.signal_gain)];
        configuration_par.low_pass_filter = [configuration_par.low_pass_filter str2double(info_signal{i}.low_pass_filter)];
        configuration_par.high_pass_filter = [configuration_par.high_pass_filter str2double(info_signal{i}.high_pass_filter)];
        
        configuration_par.angle_stim = [configuration_par.angle_stim str2double(info_signal{i}.comments(1:3))];
        if isnan(configuration_par.angle_stim(i))
            configuration_par.angle_stim = [configuration_par.angle_stim str2double(info_signal{i}.comments(1:2))];
        end
        
        aux_fsample = configuration_par.fsample(i);
        aux_signal_gain = configuration_par.signal_gain(i);
    end
    
    FileID = fopen([temp_folder{i} configuration_par.signal_path{i}]);
      
    % Loading file
    signal_par.raw_data_par{i} = fread(FileID,[configuration_par.n_channels inf],'int16');
    fclose(FileID);
    signal_par.raw_data_par{i} = signal_par.raw_data_par{i}';
    signal_par.emg_data_par{i} = signal_par.raw_data_par{i}(:,...
        [8 9 19 20 27 30 40 41 51 52 61 1 7 10 18 21 28 31 39 42 50 53 60 62 2 6 11 16 22 29 32 38 43 49 54 59 63 3 5 12 15 23 26 34 37 44 47 55 58 64 4 13 14 24 25 35 36 45 46 56 57]);
    % Channel 16 does not have signel - this line is an interpolation
    signal_par.emg_data_par{i}(:, 16) = mean(signal_par.emg_data_par{i}(:,[3 4 5 15 17 27 30 31]),2);
    
    %Coefficients for butterworth filter - passband - 25Hz a 400Hz
    [b,a] = butter(2,[25 400]*2/aux_fsample);
    % Applying butterworth filter
    ad_range = 5.0;
    conv_uv = 1000000.0;
    signal_par.emg_data_par{i} = filtfilt(b,a,(signal_par.emg_data_par{i}*ad_range*conv_uv)/(aux_signal_gain*2^(configuration_par.ad_bits)));
    
%     % Organizing electrodes into an 13x5 array
    signal_par.emg_mono_par{i} = [nan(length(signal_par.emg_data_par{i}), 1) signal_par.emg_data_par{i}(:, 1:11)...
        nan(length(signal_par.emg_data_par{i}), 1) signal_par.emg_data_par{i}(:, 12:50) nan(length(signal_par.emg_data_par{i}), 1)...
        signal_par.emg_data_par{i}(:, 51:end) nan(length(signal_par.emg_data_par{i}), 1)];
    
    signal_par.emg_diff_par{i} = diff(signal_par.emg_data_par{i},1,2);
    signal_par.emg_diff_par{i} = [nan(length(signal_par.emg_diff_par{i}), 1) signal_par.emg_diff_par{i}(:, 1:10)...
        nan(length(signal_par.emg_diff_par{i}), 1) nan(length(signal_par.emg_diff_par{i}), 1)...
        signal_par.emg_diff_par{i}(:, 12:23) nan(length(signal_par.emg_diff_par{i}), 1)...
        signal_par.emg_diff_par{i}(:, 25:36) nan(length(signal_par.emg_diff_par{i}), 1)...
        signal_par.emg_diff_par{i}(:, 38:49) nan(length(signal_par.emg_diff_par{i}), 1)...
        nan(length(signal_par.emg_diff_par{i}), 1) signal_par.emg_diff_par{i}(:, 51:end)...
        nan(length(signal_par.emg_diff_par{i}), 1) nan(length(signal_par.emg_diff_par{i}), 1)];
    
    signal_par.trigger_aux_par{i} = signal_par.raw_data_par{i}(:,65)*5/(2^(configuration_par.ad_bits-1));  
    
%         signal_par.emg_map_par{i}{j} = signal_par.emg_data_diff{i}(:, j);
        
end

% signal_par.emg_map_par = signal_par.emg_data_diff;
% signal_par = rmfield(signal_par, {'raw_data_par'});
signal_par = rmfield(signal_par, 'emg_data_par');

if length(temp_folder) == 8
    configuration_par.angle_stim = [0:45:315];
    
    signal_par.raw_data = {signal_par.raw_data_par{1}, signal_par.raw_data_par{7},...
        signal_par.raw_data_par{8}, signal_par.raw_data_par{2}, signal_par.raw_data_par{4},...
        signal_par.raw_data_par{6},signal_par.raw_data_par{3}, signal_par.raw_data_par{5}};
    signal_par = rmfield(signal_par, 'raw_data_par');
    
%     signal_par.emg_data = {signal_par.emg_data_par{1}, signal_par.emg_data_par{7},...
%         signal_par.emg_data_par{8}, signal_par.emg_data_par{2}, signal_par.emg_data_par{4},...
%         signal_par.emg_data_par{6},signal_par.emg_data_par{3}, signal_par.emg_data_par{5}};
%     signal_par = rmfield(signal_par, 'emg_data_par');
    
    signal_par.emg_mono = {signal_par.emg_mono_par{1}, signal_par.emg_mono_par{7},...
        signal_par.emg_mono_par{8}, signal_par.emg_mono_par{2}, signal_par.emg_mono_par{4},...
        signal_par.emg_mono_par{6},signal_par.emg_mono_par{3}, signal_par.emg_mono_par{5}};
    signal_par = rmfield(signal_par, 'emg_mono_par');    
    
    signal_par.emg_diff = {signal_par.emg_diff_par{1}, signal_par.emg_diff_par{7},...
        signal_par.emg_diff_par{8}, signal_par.emg_diff_par{2}, signal_par.emg_diff_par{4},...
        signal_par.emg_diff_par{6},signal_par.emg_diff_par{3}, signal_par.emg_diff_par{5}};
    signal_par = rmfield(signal_par, 'emg_diff_par');
    
    signal_par.trigger_aux = {signal_par.trigger_aux_par{1}, signal_par.trigger_aux_par{7},...
        signal_par.trigger_aux_par{8}, signal_par.trigger_aux_par{2}, signal_par.trigger_aux_par{4},...
        signal_par.trigger_aux_par{6},signal_par.trigger_aux_par{3}, signal_par.trigger_aux_par{5}};
    signal_par = rmfield(signal_par, 'trigger_aux_par');
end

% i used this because three sujects (11, 12, 13) the angles stim order is different
% from the others

% % if length(temp_folder) == 8
% %     configuration_par.angle_stim = [0:45:315];
% %     
% %     signal_par.raw_data = {signal_par.raw_data_par{1}, signal_par.raw_data_par{7},...
% %         signal_par.raw_data_par{8}, signal_par.raw_data_par{3}, signal_par.raw_data_par{5},...
% %         signal_par.raw_data_par{2},signal_par.raw_data_par{4}, signal_par.raw_data_par{6}};
% %     signal_par = rmfield(signal_par, 'raw_data_par');
% %     
% %     signal_par.emg_mono = {signal_par.emg_mono_par{1}, signal_par.emg_mono_par{7},...
% %         signal_par.emg_mono_par{8}, signal_par.emg_mono_par{3}, signal_par.emg_mono_par{5},...
% %         signal_par.emg_mono_par{2},signal_par.emg_mono_par{4}, signal_par.emg_mono_par{6}};
% %     signal_par = rmfield(signal_par, 'emg_mono_par');    
% %     
% %     signal_par.emg_diff = {signal_par.emg_diff_par{1}, signal_par.emg_diff_par{7},...
% %         signal_par.emg_diff_par{8}, signal_par.emg_diff_par{3}, signal_par.emg_diff_par{5},...
% %         signal_par.emg_diff_par{2},signal_par.emg_diff_par{4}, signal_par.emg_diff_par{6}};
% %     signal_par = rmfield(signal_par, 'emg_diff_par');
% %     
% %     signal_par.trigger_aux = {signal_par.trigger_aux_par{1}, signal_par.trigger_aux_par{7},...
% %         signal_par.trigger_aux_par{8}, signal_par.trigger_aux_par{3}, signal_par.trigger_aux_par{5},...
% %         signal_par.trigger_aux_par{2},signal_par.trigger_aux_par{4}, signal_par.trigger_aux_par{6}};
% %     signal_par = rmfield(signal_par, 'trigger_aux_par');
% % end

signal_name_default = [num2str(configuration_par.patient_id) '_' configuration_par.emg_side '_' 'signal'];
signal_name_default = strcat(signal_name_default,'.mat');

% config_name_default = [num2str(configuration_par.patient_id) '_' configuration_par.emg_side '_' 'config'];
% config_name_default = strcat(config_name_default,'.mat');

[signal_name, signal_path, ~] = uiputfile({'*.mat','MAT-files (*.mat)'},...
    'Save signal and configuration variables', signal_name_default);

% [config_name, config_path, ~] = uiputfile({'*.mat','MAT-files (*.mat)'},...
%     'Save configuration variables', config_name_default);

save([signal_path signal_name], '-struct', 'configuration_par');
save([signal_path signal_name], '-struct', 'signal_par', '-append')

for i = 1:size(file_aux,2)
    rmdir(temp_folder{i}, 's');
end

clear all

toc