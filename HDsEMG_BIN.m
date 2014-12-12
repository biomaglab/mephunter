function output = HDsEMG_BIN()
global file path signal configuration;

[file_aux path_aux] = uigetfile('*.otb','Importar Arquivo Exportado pelo OT Biolab');

path = path_aux;
file = file_aux;
signal.emg_map = {};

temp_folder = [path 'temp_folder\'];
mkdir(temp_folder);
compressed_file = [path file];
unzip(compressed_file, temp_folder);

% Preallocating variables
child_values = cell(2);
child_names = cell(2);

% Node abstract
abstract_nodes = xmlread([temp_folder 'abstract.xml']);
element_abstract = abstract_nodes.getElementsByTagName('abstract');
abstract_child_nodes = element_abstract.item(0).getChildNodes;
count = 1;
for i = 1:2:(abstract_child_nodes.getLength-2)
child_values{count} = char(abstract_child_nodes.item(i).getTextContent);
child_names{count} =  char(abstract_child_nodes.item(i).getNodeName);
count = count + 1;
end
abstract = struct(child_names{1}, child_values{1}, child_names{2}, child_values{2},...
    child_names{3}, child_values{3}, child_names{4}, child_values{4},...
    child_names{5}, child_values{5}, child_names{6}, child_values{6},...
    child_names{7}, child_values{7}, child_names{8}, child_values{8},...
    child_names{9}, child_values{9}, child_names{10}, child_values{10});

% Node signal
clear child_names child_values
element_signal = abstract_nodes.getElementsByTagName('signal');
signal_child_nodes = element_signal.item(0).getChildNodes;
count = 1;
for i = 1:2:(signal_child_nodes.getLength-2)
child_values{count} = char(signal_child_nodes.item(i).getTextContent);
child_names{count} =  char(signal_child_nodes.item(i).getNodeName);
count = count + 1;
end
info_signal = struct(child_names{1}, child_values{1}, child_names{2}, child_values{2},...
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
clear child_names child_values
abstract_2_nodes = xmlread([temp_folder info_signal.abstract_path]);
element_patient = abstract_2_nodes.getElementsByTagName('patient');
patient_child_nodes = element_patient.item(0).getChildNodes;
count = 1;
for i = 1:2:(patient_child_nodes.getLength-2)
child_values{count} = char(patient_child_nodes.item(i).getTextContent);
child_names{count} =  char(patient_child_nodes.item(i).getNodeName);
count = count + 1;
end
info_patient = struct(child_names{1}, child_values{1}, child_names{2}, child_values{2},...
    child_names{3}, child_values{3}, child_names{4}, child_values{4},...
    child_names{5}, child_values{5}, child_names{6}, child_values{6},...
    child_names{7}, child_values{7}, child_names{8}, child_values{8},...
    child_names{9}, child_values{9});

% Node xml channel
clear child_names child_values
element_channel = abstract_2_nodes.getElementsByTagName('channel');
channel_child_nodes = element_channel.item(0).getChildNodes;
count = 1;
for i = 1:2:(channel_child_nodes.getLength-2)
child_values{count} = char(channel_child_nodes.item(i).getTextContent);
child_names{count} =  char(channel_child_nodes.item(i).getNodeName);
count = count + 1;
end
info_channel = struct(child_names{1}, child_values{1}, child_names{2}, child_values{2},...
    child_names{3}, child_values{3}, child_names{4}, child_values{4},...
    child_names{5}, child_values{5}, child_names{6}, child_values{6},...
    child_names{7}, child_values{7}, child_names{8}, child_values{8});

configuration = abstract;
configuration = rmfield(configuration, 'comments');
configuration.abstract_path = info_signal.abstract_path;
configuration.signal_path = info_signal.signal_path;
configuration.channels =  str2double(info_signal.channels);
configuration.fsample =  str2double(info_signal.fsample);
configuration.signal_gain =  str2double(info_signal.signal_gain);
configuration.low_pass_filter =  str2double(info_signal.low_pass_filter);
configuration.high_pass_filter =  str2double(info_signal.high_pass_filter);
configuration.ad_bits =  str2double(info_signal.ad_bits);
configuration.patient_id = str2double(info_patient.id);
configuration.pathology = info_patient.pathology;
configuration.emg_side = info_channel.side;

configuration.angle_stim = str2double(info_signal.comments(1:3));
if isnan(configuration.angle_stim)
    configuration.angle_stim = str2double(info_signal.comments(1:2));
end

%Coefficients for butterworth filter - passband - 50Hz a 400Hz
% [b,a] = butter(2,[25 400]*2/configuration.fsample);

FileID = fopen([temp_folder configuration.signal_path]);

% Loading file
signal.emg_data = fread(FileID,[configuration.channels inf],'int16');
fclose(FileID);
signal.emg_data = signal.emg_data';
signal.raw_data = signal.emg_data;
signal.emg_data = signal.emg_data(:,...
    [8 9 19 20 27 30 40 41 51 52 61 1 7 10 18 21 28 31 39 42 50 53 60 62 2 6 11 16 22 29 32 38 43 49 54 59 63 3 5 12 15 23 26 34 37 44 47 55 58 64 4 13 14 24 25 35 36 45 46 56 57]);
% Channel 16 does not have signel - this line is an interpolation
signal.emg_data(:, 16) = mean(signal.emg_data(:,[3 4 5 15 17 27 30 31]),2);

% Applying butterworth filter
% signal.emg_data = filtfilt(b,a,signal.emg_data*2500/(configuration.signal_gain*2^11));
% Organizing electrodes into an 13x5 array

signal.emg_data_mono = [nan(length(signal.emg_data), 1) signal.emg_data(:, 1:11)...
    nan(length(signal.emg_data), 1) signal.emg_data(:, 12:50) nan(length(signal.emg_data), 1)...
    signal.emg_data(:, 51:end) nan(length(signal.emg_data), 1)];

% for i = 1:size(signal.emg_data, 2)
%      signal.emg_map{i} = signal.emg_data(:, i);
% end

% - Novo
signal.emg_data_diff = diff(signal.emg_data,1,2);
signal.emg_data_diff = [nan(length(signal.emg_data_diff), 1) signal.emg_data_diff(:, 1:10)...
    nan(length(signal.emg_data_diff), 1) nan(length(signal.emg_data_diff), 1) signal.emg_data_diff(:, 12:23)...
    nan(length(signal.emg_data_diff), 1) signal.emg_data_diff(:, 25:36)...
    nan(length(signal.emg_data_diff), 1) signal.emg_data_diff(:, 38:49)...
    nan(length(signal.emg_data_diff), 1) nan(length(signal.emg_data_diff), 1)...
    signal.emg_data_diff(:, 51:end) nan(length(signal.emg_data_diff), 1) nan(length(signal.emg_data_diff), 1)];

signal.emg_data_all = {signal.emg_data_mono signal.emg_data_diff};

for i = 1:65
    for j = 1:2
        signal.emg_map{i}(:, j) = signal.emg_data_all{j}(:, i);
    end     
end

configuration.emg_mode = {'Mono', 'Diff'};

rmdir(temp_folder, 's');

output.signal = signal;
output.configuration = configuration;