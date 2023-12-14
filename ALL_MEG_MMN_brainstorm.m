%% Complete Project_MMN analysis script (modify first two sections for Baseline, 3mo, 6mo and 1yr)

%% (OPTIONAL) Run in server mode: Set Brainstorm nogui on completely headless mode (Baseline) 
% Set a window in 99 first, then run this section and run scripts normally 
% (it won't break despite bad remote connection this way -16/03/2020-)

% Set up the Brainstorm files
clear
% addpath('~/matlab/brainstorm3'); 
addpath('~/matlab/brainstorm3_v20220706');
BrainstormDbDir = '~/brainstorm_db';
% Start Brainstorm
if ~brainstorm('status')
    brainstorm server
end
bst_set('BrainstormDbDir',BrainstormDbDir)
% Select the correct protocol
ProtocolName = 'Project_MMN_baseline'; % Enter the name of your protocol
sProtocol.Comment = ProtocolName;
% sProtocol.SUBJECTS = [home 'anat'];
% /private/path/Project/brainstorm_db/Project_MMN
% or ~/brainstorm_db/Project_MMN_baseline
sProtocol.SUBJECTS = '~/brainstorm_db/Project_FEManat_baseline/anat';
sProtocol.STUDIES = '~/brainstorm_db/Project_MMN/data';
db_edit_protocol('load',sProtocol);
% Get the protocol index
iProtocol = bst_get('Protocol', ProtocolName);
if isempty(iProtocol)
    error(['Unknown protocol: ' ProtocolName]);
end
% Select the current procotol
gui_brainstorm('SetCurrentProtocol', iProtocol);

%% Define variables Project_MMN_Baseline 

% Directories for server in local computer
running_in = 'local'; % 'server' 'local'

% Define paths
if strcmp(running_in,'server')
    root_dir = '/private/path/Project/User/Project_MMN_baseline'; 
    root_dir_bs = '~/brainstorm_db/Project_MMN_baseline'; 
    anat_path = '/private/path/Project/brainstorm_db/Project_FEManat_baseline';
elseif strcmp(running_in,'local')
    root_dir = 'X:/Project/User/Project_MMN_baseline'; 
    root_dir_bs = 'X:/Project/brainstorm_db/Project_MMN'; 
    anat_path = 'X:/Project/brainstorm_db/Project_FEManat_baseline';
end
addpath([root_dir '/Scripts']);

% Get protocol name, in case we don't run it with server mode
pos_last = find(root_dir_bs == '/', 1, 'last');
ProtocolName = root_dir_bs(pos_last+1:end); clear pos_last
load([root_dir '/Final_subject_array_march_22_FEALL.mat'])

% Define participants
participant = {subject_array{:,1}};
% So that it has a different name if we loop through sections
participant_general_list = participant; 
gavr_name = 'GAVR_march_2022'; % Whatever you want the GAVR folder to be named

participant_group = {'FE','C'};
condition = {'1','2','3'};
condition_names = {'Standard','DeviantPitch','DeviantDuration'};
condition_mismatch_names = {'Standard','DeviantPitch','DeviantDuration','DeviantPitch-Standard','DeviantDuration-Standard'};
condition_short_labels = {'STD','PDev','DDev','Pitch_MMN','Dur_MMN'};
modality_data = {'EEG','MEG','BIMODAL'};
epoch_wave = [-0.05, 0.3]; % LLR only
epoch_baseline = [-0.05, 0]; % LLR baseline correction
time_noise_cov = [-0.05, 0]; % time window from which to obtain noise covariance values
reject_EEG_absolute = [0, 50]; % absolute threshold
reject_MEG_GRAD_absolute = [0, 2500];
reject_MEG_MAG_absolute = [0, 2500];
LLR_highpass = 0; 
LLR_lowpass = 20; %
tranband_LLR = 0; % 
crit_sweeps = 30; % minimum number of surviving sweeps to discard EEG, MEG or BIMODAL
crit_percent = 50; % minimum percentage of surviving sweeps to discard EEG, MEG or BIMODAL
reref_option = 1; % 0 = NO 1 = YES Yes or no rereference
ref_EEG = 'M1, M2'; % in case we use it, but we won't for now: Alternative: 'AVERAGE'
sensor_analysis = 4; % 1 = EEG only; 2 = MEG only; 3 = Combined EEG and MEG only; 4 = ALL
source_analysis = 2; % 1 = EEG only; 2 = MEG only; 3 = Combined EEG and MEG only. 4 = ALL
source_space = {'Vol','Surf'}; % Which sources to average (Volume,  Surface or both)
source_noise_option = 'reg'; % For 'NoiseMethod'; 'median'; 'diag'
source_noise_tag = 'Source_Regul'; % 'Source_Eigen' 'Source_diag'
source_inverse_measure = 'dspm2018'; % 'amplitude' or 'dspm2018'
delete_previous_file = 1; % 1 = delete previous head models and source kernels if reran these steps
compute_covariance = 1; % In case we want to repeat sources but not recompute noise covariance matrices
name_volume_source_grid = 'headmodel_grid.mat';
group_default_cortex = 'tess_cortex_pial_02.mat'; % cortical surface where all cortical sources of invididual subjects will be projected
group_default_volume = 'tess_cortex_mixed_02.mat'; % volume surface with brainstem and cerebellum

% If selected more than one, they will be averaged as a cluster
choice_channel_EEG = {'FCz'}; % Watch out because others may be missing for some subjects
choice_channel_MEG = {'MEG2411'};
Peaks_to_extract = {'MMN'}; % Add as many as needed, with time windows below
time_window_scalp_Pitch_MMN_EEG = [70 170]; % (ms) Add as many as needed with different names
time_window_scalp_Dur_MMN_EEG = [110 210]; % (ms) Add as many as needed with different names
time_window_scalp_Pitch_MMN_MEG = [70 170]; % (ms) Add as many as needed with different names
time_window_scalp_Dur_MMN_MEG = [110 210]; % (ms) Add as many as needed with different names
time_window_source_left_Pitch_MMN = [70 170]; % (ms)
time_window_source_right_Pitch_MMN = [70 170]; % (ms)
time_window_source_left_Dur_MMN = [110 210]; % (ms)
time_window_source_right_Dur_MMN = [110 210]; % (ms)
% TRY Project MMN WITH LONGER TIME WINDOWS TOO (+-25 or +-50 around peak) 80-90 to 160

% Original ones
% time_window_scalp_Pitch_MMN_EEG = [110 130]; % (ms) Add as many as needed with different names
% time_window_scalp_Dur_MMN_EEG = [150 170]; % (ms) Add as many as needed with different names
% time_window_scalp_Pitch_MMN_MEG = [110 130]; % (ms) Add as many as needed with different names
% time_window_scalp_Dur_MMN_MEG = [150 170]; % (ms) Add as many as needed with different names
% time_window_source_left_Pitch_MMN = [110 130]; % (ms)
% time_window_source_right_Pitch_MMN = [110 130]; % (ms)
% time_window_source_left_Dur_MMN = [150 170]; % (ms)
% time_window_source_right_Dur_MMN = [150 170]; % (ms)

% +/- 25ms from peak
% time_window_scalp_Pitch_MMN_EEG = [95 145]; % (ms) Add as many as needed with different names
% time_window_scalp_Dur_MMN_EEG = [135 185]; % (ms) Add as many as needed with different names
% time_window_scalp_Pitch_MMN_MEG = [95 145]; % (ms) Add as many as needed with different names
% time_window_scalp_Dur_MMN_MEG = [135 185]; % (ms) Add as many as needed with different names
% time_window_source_left_Pitch_MMN = [95 145]; % (ms)
% time_window_source_right_Pitch_MMN = [95 145]; % (ms)
% time_window_source_left_Dur_MMN = [135 185]; % (ms)
% time_window_source_right_Dur_MMN = [135 185]; % (ms)

% +/- 50s from peak
% time_window_scalp_Pitch_MMN_EEG = [70 170]; % (ms) Add as many as needed with different names
% time_window_scalp_Dur_MMN_EEG = [110 210]; % (ms) Add as many as needed with different names
% time_window_scalp_Pitch_MMN_MEG = [70 170]; % (ms) Add as many as needed with different names
% time_window_scalp_Dur_MMN_MEG = [110 210]; % (ms) Add as many as needed with different names
% time_window_source_left_Pitch_MMN = [70 170]; % (ms)
% time_window_source_right_Pitch_MMN = [70 170]; % (ms)
% time_window_source_left_Dur_MMN = [110 210]; % (ms)
% time_window_source_right_Dur_MMN = [110 210]; % (ms)

header_Project_MMN = {'FE_Pitch_MMN','FE_Dur_MMN','C_Pitch_MMN','C_Dur_MMN'}; % Define based on previous, will determine SPSS matrix

initialVars = who; % variables up until here, which won't be deleted afterwards
initialVars = who; % twice so that InitialVars itself is not deleted

%% Import raw data (ALWAYS IN server)
 
tic
disp(' ');      
disp('-------------------------');  
disp('IMPORTING EEG/MEG DATA FOR Project MMN Baseline');
disp(datetime)
disp('-------------------------');     
disp(' '); 

for p = 1:length(participant)
    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p})); %#ok<*CCAT1>
    if strcmp(subject_array{pos_subj,3},'needs_import')


    folders = dir(['/private/path/Project/analysis/' participant{p} '/']);
    infolder = find(endsWith({folders.name},'tsss_AMICA.fif') & (...
                contains({folders.name},[participant{p} '_MMN'])));
        if isempty(infolder)
            error(['No MMN AMICA files for ' participant{p}]);
        end

        for i = 1:length(infolder)
            line = infolder(i);
            file_name = ['/private/path/Project/analysis/' participant{p} '/' folders(line).name];  %#ok<*SAGROW>             
            
            disp(' ');      
            disp('-------------------------');  
            disp(['Importing data EEG/MEG LLR data for ' participant{p}]);
            disp(datetime)
            disp(' '); 

            sFiles = [];
            % Process: Create link to raw file
            sFiles = bst_process('CallProcess', 'process_import_data_raw', sFiles, [], ...
                'subjectname',    participant{p}, ...
                'datafile',       {file_name, 'FIF'}, ...
                'channelreplace', 0, ...
                'channelalign',   0, ...
                'evtmode',        'value');    
            
            % If successful, update subject_array for this subject
            subject_array{pos_subj,3} = 'needs_events';
            save([root_dir '/subject_array.mat'],'subject_array')  
        end
    end
end

clearvars('-except', initialVars{:});

disp 'DONE WITH IMPORTING DATA EEG/MEG FOR Project MMN Baseline!!!'
disp(datetime)
toc

%% Remove ANY previous projectors and import events (ALWAYS IN server) 

tic
disp(' ');      
disp('-------------------------');  
disp('REMOVING PROJECTORS AND IMPORTING EVENTS (Project MMN)');  
disp(datetime)
disp('-------------------------');     
disp(' '); 

for p = 1:length(participant)
    
    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    if strcmp(subject_array{pos_subj,3},'needs_events')
   
    folders = dir([root_dir_bs '/data/' participant{p} '/']);
    infolder = find(endsWith({folders.name},'tsss_AMICA') & (...
    contains({folders.name},[participant{p} '_MMN'])));
    if isempty(infolder)
        error(['No MMN AMICA folders for ' participant{p}]);
    end

    % Remove active projectors
    for i = 1:length(infolder)
        line = infolder(i);
        sFiles = [root_dir_bs '/data/' participant{p} '/' folders(line).name '/channel_vectorview306_acc1.mat'];  %#ok<*SAGROW>             
        eval(['load ' sFiles])
        variableInfo = who('-file',sFiles);

        disp(' ');      
        disp('-------------------------');  
        disp(['Removing active projectors before importing events: ' participant{p}]);
        disp(datetime)
        disp(' '); 
        
        if ~isempty(Projector)
            num_proj = size(Projector,2);
            for npj = 1:num_proj
                Projector(npj).Status = 0;
            end
        end
        save(sFiles,variableInfo{:});
    end
    
    % Now we are going to import events, so reload the subject
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' '); 

    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);

    % Import events now
    folders = dir([root_dir_bs '/data/' participant{p} '/']);
    infolder =  find(contains({folders.name},'raw') & contains({folders.name},[participant{p}]) & endsWith({folders.name},'tsss_AMICA'));
    if isempty(infolder)
        continue % Go to next subject
    end
    sFiles = {};
    for i= 1:length(infolder)
        line = infolder(i);
        sFiles{i} = [participant{p} '/' folders(line).name '/data_0raw_' folders(line).name(5:end) '.mat'];
    end            
    if isempty(sFiles)
        error(['No AMICA filenames for ' participant{p}])
    end
    
    % Only if this participant has EEG channels
    if ~strcmp(subject_array{pos_subj,4},'no_EEG_chans')
        % Process: CNRL rename EEG Channels
        sFiles = bst_process('CallProcess', 'process_cnrl_rename_eeg', sFiles, [], ...
            'action', 2);  % EasyCap

        % Process: Set channels type
        sFiles = bst_process('CallProcess', 'process_channel_settype', sFiles, [], ...
        'sensortypes', 'HEOG,VEOG', ...
        'newtype',     'EOG');

        % Process: Set channels type
        sFiles = bst_process('CallProcess', 'process_channel_settype', sFiles, [], ...
        'sensortypes', 'ECG', ...
        'newtype',     'ECG');

        if reref_option == 1     
        % Process: Re-reference EEG
        sFiles = bst_process('CallProcess', 'process_eegref', sFiles, [], ...
            'eegref',      ref_EEG, ...
            'sensortypes', 'EEG');     
        end

        % Process: Set channels type: addition on December 2020
        sFiles = bst_process('CallProcess', 'process_channel_settype', sFiles, [], ...
            'sensortypes', 'M1,M2', ...
            'newtype',     'EEG REF');
    end
     
    % Process: Run Matlab command
    sFiles = bst_process('CallProcess', 'process_matlab_eval', sFiles, [], ...
        'matlab',      ['% Available variables: Data, TimeVector' 10 '' 10 'TimeVector = TimeVector - TimeVector(1);' 10 ''], ...
        'sensortypes', '');

    % Process: Import events file
    sFiles = bst_process('CallProcess', 'process_CNRL_evt_import', sFiles, [], ...
        'evtname', '1, 2, 3, boundary');

    % Process: Convert to simple event
    sFiles = bst_process('CallProcess', 'process_evt_simple', sFiles, [], ...
        'eventname', 'boundary', ...
        'method',    1);  % Keep the start of the events

    % Process: Rename event
    sFiles = bst_process('CallProcess', 'process_evt_rename', sFiles, [], ...
        'src',  'boundary', ...
        'dest', 'boundaryStart');

    % Process: Import events file
    sFiles = bst_process('CallProcess', 'process_CNRL_evt_import', sFiles, [], ...
        'evtname', 'boundary');

    % Process: Convert to simple event
    sFiles = bst_process('CallProcess', 'process_evt_simple', sFiles, [], ...
        'eventname', 'boundary', ...
        'method',    3);  % Keep the end of the events

    % Process: Rename event
    sFiles = bst_process('CallProcess', 'process_evt_rename', sFiles, [], ...
        'src',  'boundary', ...
        'dest', 'boundaryEnd');

    for c = 1:length(condition)
    % Process: Remove simultaneous End
    sFiles = bst_process('CallProcess', 'process_evt_remove_simult', sFiles, [], ...
        'remove', condition{c}, ...
        'target', 'boundaryEnd', ...
        'dt',     0.05, ...
        'rename', 0);

    % Process: Remove simultaneous Start
    sFiles = bst_process('CallProcess', 'process_evt_remove_simult', sFiles, [], ...
        'remove', condition{c}, ...
        'target', 'boundaryStart', ...
        'dt',     0.3, ...
        'rename', 0);
    end

    % Process: scale MEG values (only for files needing it)
    if strcmp(subject_array{pos_subj,4},'scale')
        sFiles = bst_process('CallProcess', 'process_scale', sFiles, [], ...
            'factor',      3, ...
            'sensortypes', 'MEG');
    end
    
    % If successful, update subject_array for this subject
    subject_array{pos_subj,3} = 'needs_filter';
    save([root_dir '/subject_array.mat'],'subject_array')
    
    end
end

clearvars('-except', initialVars{:});
disp 'DONE WITH REMOVING PROJECTORS AND IMPORTING EVENTS (Project MMN)!!!'
disp(datetime)
toc

%% Filter

tic
disp(' ');      
disp('-------------------------');  
disp('FILTERING (Project MMN)');  
disp(datetime)
disp('-------------------------');     
disp(' '); 

for p = 1:length(participant)
    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    if strcmp(subject_array{pos_subj,3},'needs_filter')
    
    % Reload subject first
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' '); 

    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
        
    % Find files    
    folders = dir([root_dir_bs '/data/' participant{p} '/']);
    if strcmp(subject_array{pos_subj,4},'scale')
        infolder = find(contains({folders.name},'raw') & contains({folders.name},[participant{p} '_MMN']) & endsWith({folders.name},'matlab_scale'));
    else
        infolder = find(contains({folders.name},'raw') & contains({folders.name},[participant{p} '_MMN']) & endsWith({folders.name},'matlab'));
    end
    if isempty(infolder)
        error(['No MMN matlab or matlab_scale folders for ' participant{p}]);
    end  
        
    sFiles = {};
    for i= 1:length(infolder)
        line = infolder(i);
        sFiles{i} = [participant{p} '/' folders(line).name '/data_0raw_' folders(line).name(5:end) '.mat'];
    end

    disp(' ');      
    disp('-------------------------');  
    disp(['Filtering ' participant{p}]);
    disp(datetime)
    disp(' '); 
    
    try
    % Process: LLR filter (low pass)
    sFiles = bst_process('CallProcess', 'process_bandpass', sFiles, [], ...
        'sensortypes', '', ...
        'highpass',    LLR_highpass, ... 
        'lowpass',     LLR_lowpass, ... 
        'tranband',    tranband_LLR, ...
        'attenuation', 'strict', ...  % 60dB
        'ver',         '2019', ...  % 2019
        'mirror',      0, ...
        'read_all',    0);
    
    % If successful, update subject_array for this subject
    subject_array{pos_subj,3} = 'needs_epoch';
    save([root_dir '/subject_array.mat'],'subject_array')    
    catch
        warning(['No filtering performed for ' participant{p}]);
        continue;
    end
    
    end  
end

clearvars('-except', initialVars{:});
disp 'DONE WITH FILTERING (Project MMN)!!!'
disp(datetime)
toc

%% Epoch and amplitude threshold

tic
disp(' ');      
disp('-------------------------');  
disp('EPOCHING (Project MMN)');  
disp(datetime)
disp('-------------------------');     
disp(' '); 

for p = 1:length(participant)
    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    if strcmp(subject_array{pos_subj,3},'needs_epoch')
    
    % Reload subject first
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' '); 

    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
        
    % Find files    
    folders = dir([root_dir_bs '/data/' participant{p} '/']);
    infolder = find(contains({folders.name},'raw') & contains({folders.name},[participant{p} '_MMN']) & endsWith({folders.name},'low'));

    if isempty(infolder)
        error(['No folder names ending in low for ' participant{p}]);
    end  
        
    sFiles = {};
    for i= 1:length(infolder)
        line = infolder(i);
        sFiles{i} = [participant{p} '/' folders(line).name '/data_0raw_' folders(line).name(5:end) '.mat'];
    end
    
    disp(' ');      
    disp('-------------------------');
    disp(['Renaming events to epoch data: participant ' participant{p}]);
    disp(datetime)
    disp(' ');  

    % Rename events (if they are already renamed, nothing will happen)
    for c =1:length(condition)
    sFiles = bst_process('CallProcess', 'process_evt_rename', sFiles, [], ...
        'src',  condition{c}, ...
        'dest', condition_names{c});
    end    

    disp(' ');      
    disp('-------------------------');  
    disp(['Making epochs for ' participant{p}]);
    disp(datetime)
    disp(' ');  

    % Process: epoch data for normal epochs
    sFiles = bst_process('CallProcess', 'process_import_data_event', sFiles, [], ...
        'subjectname',  participant{p}, ...
        'condition',    '', ...
    ...%    'datafile',     RawFiles, ...
        'eventname',    'Standard, DeviantPitch, DeviantDuration', ...
        'timewindow',   [], ...
        'epochtime',    epoch_wave, ...
        'createcond',   1, ...
        'ignoreshort',  1, ...
        'channelalign', 0, ...
        'usectfcomp',   0, ...
        'usessp',       1, ...
        'freq',         [], ...
        'baseline',     []); 
    
    disp(' ');      
    disp('-------------------------');  
    disp(['Baseline correcting epochs for ' participant{p}]);
    disp(datetime)
    disp(' ');  

    % Process: DC offset correction: [-50ms,1ms]
    sFiles = bst_process('CallProcess', 'process_baseline_norm', sFiles, [], ...
        'baseline',    epoch_baseline, ...
        'sensortypes', '', ...
        'method',      'bl', ...  % DC offset correction:    x_std = x - &mu;
        'overwrite',   1);
    
    disp(' ');      
    disp('-------------------------');  
    disp(['Cleaning epochs for ' participant{p} '(MEG)']);
    disp(datetime)
    disp(' '); 
    
    % Process: Detect bad trials: Absolute threshold MEG
    sFiles_MEG = bst_process('CallProcess', 'process_CNRL_detectbad',sFiles, [], ...
        'timewindow', [], ...
        'meggrad',    reject_MEG_GRAD_absolute, ...
        'megmag',     reject_MEG_MAG_absolute, ...
        'eeg',        [0, 0], ...
        'ieeg',       [0, 0], ...
        'eog',        [0, 0], ...
        'ecg',        [0, 0], ...
        'rejectmode', 2);  % Reject the entire trial
    
    % SENSOR AVERAGE MEG  
    disp(' ');      
    disp('-------------------------');  
    disp(['Averaging epochs for ' participant{p} '(MEG)']);
    disp(datetime)
    disp(' '); 

    sFiles_MEG = bst_process('CallProcess', 'process_average', sFiles_MEG, [], ...
        'avgtype',         5, ...  % By trial group (folder average)
        'avg_func',        1, ...  % Arithmetic average:  mean(x)
        'weighted',        0, ...
        'keepevents', 0);
    
    % Process: Add tag
    sFiles_MEG = bst_process('CallProcess', 'process_add_tag', sFiles_MEG, [], ...
        'tag',           'MEG_average', ...
        'output',        2);  % Add to file name (1 to add a tag)

    % Process: Set name
    sFiles_MEG = bst_process('CallProcess', 'process_set_comment', sFiles_MEG, [], ...
        'tag',           'MEG_average', ...
        'isindex',       1);
    
    % Reset BadTrials variable and reload folder (Std, DevP, DevDur) 
    for c = 1:length(condition_names)
        trials_file = [root_dir_bs '/data/' participant{p} '/' condition_names{c} '/brainstormstudy.mat'];
        if ~exist(trials_file,'file')
            % If there is no file/folder, continue
            warning(['No ' condition_names{c} ' folder for ' participant{p}]);
            continue;
        end  
        eval(['load ' trials_file]); % load them
        variableInfo = who('-file',trials_file);
        % Reset them and save back to original file
        BadTrials = cell([], 1); % exact structure it has when empty
        save(trials_file,variableInfo{:});
        % Reload folder
        [~, iStudies] = bst_get('StudyWithCondition', [participant{p} filesep condition_names{c}]);
        db_reload_studies(iStudies);
    end
    
    % EEG and BIMODAL cleaning and averages will only take place if 
    % that subject doesn't have a 'no_EEG_chans' label
    if ~strcmp(subject_array{pos_subj,4},'no_EEG_chans')
    
    disp(' ');      
    disp('-------------------------');  
    disp(['Cleaning epochs for ' participant{p} '(EEG)']);
    disp(datetime)
    disp(' '); 

    % Process: Detect bad trials: Absolute threshold EEG
    sFiles_EEG = bst_process('CallProcess', 'process_CNRL_detectbad', sFiles, [], ...
        'timewindow', [], ...
        'meggrad',    [0, 0], ...
        'megmag',     [0, 0], ...
        'eeg',        reject_EEG_absolute, ...
        'ieeg',       [0, 0], ...
        'eog',        [0, 0], ...
        'ecg',        [0, 0], ...
        'rejectmode', 2);  % Reject the entire trial
    
    % SENSOR AVERAGE EEG  
    disp(' ');      
    disp('-------------------------');  
    disp(['Averaging epochs for ' participant{p} '(EEG)']);
    disp(datetime)
    disp(' '); 

    sFiles_EEG = bst_process('CallProcess', 'process_average', sFiles_EEG, [], ...
        'avgtype',         5, ...  % By trial group (folder average)
        'avg_func',        1, ...  % Arithmetic average:  mean(x)
        'weighted',        0, ...
        'keepevents', 0);
    
    % Process: Add tag
    sFiles_EEG = bst_process('CallProcess', 'process_add_tag', sFiles_EEG, [], ...
        'tag',           'EEG_average', ...
        'output',        2);  % Add to file name (1 to add a tag)

    % Process: Set name
    sFiles_EEG = bst_process('CallProcess', 'process_set_comment', sFiles_EEG, [], ...
        'tag',           'EEG_average', ...
        'isindex',       1);
    
    % Reset BadTrials variable and reload folder (Std, DevP, DevDur) 
    for c = 1:length(condition_names)
        trials_file = [root_dir_bs '/data/' participant{p} '/' condition_names{c} '/brainstormstudy.mat'];
        if ~exist(trials_file,'file')
            % If there is no file/folder, continue
            warning(['No ' condition_names{c} ' folder for ' participant{p}]);
            continue;
        end  
        eval(['load ' trials_file]); % load them
        variableInfo = who('-file',trials_file);
        % Reset them and save back to original file
        BadTrials = cell([], 1); % exact structure it has when empty
        save(trials_file,variableInfo{:});
        % Reload folder
        [~, iStudies] = bst_get('StudyWithCondition', [participant{p} filesep condition_names{c}]);
        db_reload_studies(iStudies);
    end
    
    disp(' ');      
    disp('-------------------------');  
    disp(['Cleaning epochs for ' participant{p} '(BIMODAL)']);
    disp(datetime)
    disp(' '); 
    
    % Process: Detect bad trials: Abolute threshold BOTH EEG AND MEG
    sFiles_BIMODAL = bst_process('CallProcess', 'process_CNRL_detectbad',sFiles, [], ...
        'timewindow', [], ...
        'meggrad',    reject_MEG_GRAD_absolute, ...
        'megmag',     reject_MEG_MAG_absolute, ...
        'eeg',        reject_EEG_absolute, ...
        'ieeg',       [0, 0], ...
        'eog',        [0, 0], ...
        'ecg',        [0, 0], ...
        'rejectmode', 2);  % Reject the entire trial
    
    % SENSOR AVERAGE BIMODAL  
    disp(' ');      
    disp('-------------------------');  
    disp(['Averaging epochs for ' participant{p} '(BIMODAL)']);
    disp(datetime)
    disp(' '); 

    sFiles_BIMODAL = bst_process('CallProcess', 'process_average', sFiles_BIMODAL, [], ...
        'avgtype',         5, ...  % By trial group (folder average)
        'avg_func',        1, ...  % Arithmetic average:  mean(x)
        'weighted',        0, ...
        'keepevents', 0);
    
    % Process: Add tag
    sFiles_BIMODAL = bst_process('CallProcess', 'process_add_tag', sFiles_BIMODAL, [], ...
        'tag',           'BIMODAL_average', ...
        'output',        2);  % Add to file name (1 to add a tag)

    % Process: Set name
    sFiles_BIMODAL = bst_process('CallProcess', 'process_set_comment', sFiles_BIMODAL, [], ...
        'tag',           'BIMODAL_average', ...
        'isindex',       1);
    
    % Reset BadTrials variable and reload folder (Std, DevP, DevDur) 
    for c = 1:length(condition_names)
        trials_file = [root_dir_bs '/data/' participant{p} '/' condition_names{c} '/brainstormstudy.mat'];
        if ~exist(trials_file,'file')
            % If there is no file/folder, continue
            warning(['No ' condition_names{c} ' folder for ' participant{p}]);
            continue;
        end  
        eval(['load ' trials_file]); % load them
        variableInfo = who('-file',trials_file);
        % Reset them and save back to original file
        BadTrials = cell([], 1); % exact structure it has when empty
        save(trials_file,variableInfo{:});
        % Reload folder
        [~, iStudies] = bst_get('StudyWithCondition', [participant{p} filesep condition_names{c}]);
        db_reload_studies(iStudies);
    end
    
    end
    
    % If successful, update subject_array for this subject
    subject_array{pos_subj,3} = 'needs_mismatch';
    save([root_dir '/subject_array.mat'],'subject_array') 
    
    end  
end

clearvars('-except', initialVars{:});
disp 'DONE WITH EPOCHING (Project MMN)!!!'
disp(datetime)
toc

%% Obtain Mismatch 

tic
disp(' ');      
disp('-------------------------');  
disp('OBTAINING MISMATCH (DEV - STD: Project MMN Baseline)');  
disp(datetime)
disp('-------------------------');     
disp(' '); 

% Obtain MMN
for p = 1:length(participant)
    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    if strcmp(subject_array{pos_subj,3},'needs_mismatch')
    
    % Reload subject epoch folders first
    disp(' ');      
    disp('-------------------------');
    disp(['loading epoch folders for participant ' participant{p}]);
    disp(datetime)
    disp(' '); 
    for c = 1:length(condition_names)
        [~, iStudies] = bst_get('StudyWithCondition', [participant{p} filesep condition_names{c}]);
        db_reload_studies(iStudies);
    end
    
    % Do once for each modality (EEG, MEG, BIMODAL)
    for mode = 1:length(modality_data)
    % Define standard sFile
    files = dir([root_dir_bs '/data/' participant{p} '/' condition_names{1}]);
    if isempty(files)
        error(['No ' condition_names{1} ' files for ' participant{p}]);
    end
    infolder = find(endsWith({files.name},[modality_data{mode} '_average.mat'])); 
    if isempty(infolder)
       % In case, for instance, there is no EEG data for this subject
       warning(['No ' condition_names{1} ' ' modality_data{mode} ' average for ' participant{p}]);
       continue;
    end  
    if length(infolder) > 1
        error(['More than one ' condition_names{1} ' ' modality_data{mode} ' average for ' participant{p}]);
    end
    filename = [participant{p} '/' condition_names{1} '/' files(infolder).name];
    sFiles2 = {}; sFiles2{1} = filename;
    % Do once for DevPitch and again for DevDur
    for c = 2:length(condition_names)
        files = dir([root_dir_bs '/data/' participant{p} '/' condition_names{c}]);
        if isempty(files)
            error(['No ' condition_names{c} ' files for ' participant{p}]);
        end
        infolder = find(endsWith({files.name},[modality_data{mode} '_average.mat'])); 
        if isempty(infolder)
            % In case, for instance, there is no EEG data for this subject
            warning(['No ' condition_names{c} ' ' modality_data{mode} ' average for ' participant{p}]);
            continue;
        end  
        if length(infolder) > 1
            error(['More than one ' condition_names{c} ' ' modality_data{mode} ' average for ' participant{p}]);
        end
        filename = [participant{p} '/' condition_names{c} '/' files(infolder).name];
        sFiles = {}; sFiles{1} = filename;
        
        disp(' ');      
        disp('-------------------------');  
        disp(['Calculating ' condition_names{c} ' MMN for ' participant{p} '(' modality_data{mode} ')']);
        disp(datetime)
        disp(' '); 
        
        % Process: Difference: A-B
        sFiles_mismatch = bst_process('CallProcess', 'process_diff_ab', sFiles, sFiles2);
        
        % Give outcome file a new name
        % Process: Add tag
        sFiles_mismatch = bst_process('CallProcess', 'process_add_tag', sFiles_mismatch, [], ...
        'tag',           ['_' modality_data{mode} '_' condition_names{c} '_MMN'], ...
        'output',        2);  % Add to file name (1 to add a tag)

        % Process: Set name
        sFiles_mismatch = bst_process('CallProcess', 'process_set_comment', sFiles_mismatch, [], ...
        'tag',           [modality_data{mode} '_' condition_names{c} '_MMN'], ...
        'isindex',       1);
    end
    end
    % If successful, update subject_array for this subject
    subject_array{pos_subj,3} = 'needs_forward';
    save([root_dir '/subject_array.mat'],'subject_array')    
    end  
end

clearvars('-except', initialVars{:});
disp 'DONE WITH OBTAINING MISMATCH (DEV - STD: Project MMN Baseline)!!!'
disp(datetime)
toc

%% Store info about number of surviving trials (percentage lost, insuff covariance trials, etc)

tic
disp(' ');      
disp('-------------------------');  
disp('OBTAINING SWEEP COUNT (Project MMN)');  
disp(datetime)
disp('-------------------------');     
disp(' '); 

colNames = {'Group','STD_N','PD_N','DD_N','STD_EEG','PD_EEG','DD_EEG','STD_MEG','PD_MEG','DD_MEG','STD_BI','PD_BI','DD_BI'};
Trial_count = {};
Percentage_count = [];
for p = 1:length(participant)
    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    Subj_type = subject_array{pos_subj,2};
    
    % If subject does not have any of these markers, don't include it in
    % the sweep count (because it may be a bad subject, etc)
    if ~strcmp(subject_array{pos_subj,3},'needs_forward') && ~strcmp(subject_array{pos_subj,3},'needs_sources') && ~strcmp(subject_array{pos_subj,3},'sources_finished')
        continue;
    end
    
    % Do once for DevPitch, DevDur and Std
    for c = 1:length(condition_names)
    % Do once for each modality (EEG, MEG, BIMODAL)
    files = dir([root_dir_bs '/data/' participant{p} '/' condition_names{c}]);
    if isempty(files)
        error(['No ' condition_names{c} ' files for ' participant{p}]);
    end
    number = find(contains({files.name},'_trial')); 
    % Original number of trials found in folder
    eval(['ot_' condition_names{c} ' = length(number);'])
    for mode = 1:length(modality_data)
        infolder = find(endsWith({files.name},[modality_data{mode} '_average.mat']));
        if length(infolder) > 1
            error(['More than one ' condition_names{c} ' ' modality_data{mode} ' average for ' participant{p}]);
        end
        if isempty(infolder)
            warning(['No ' condition_names{c} ' ' modality_data{mode} ' average for ' participant{p}]);
            % Surviving trials
            eval(['st_' modality_data{mode} '_' condition_names{c} ' = 0;'])
            % Percentage of trials
            eval(['pt_' modality_data{mode} '_' condition_names{c} ' = 0%;'])
        else
            filename = [root_dir_bs '/data/' participant{p} '/' condition_names{c} '/' files(infolder).name];
            load(filename);
            % Surviving trials
            eval(['st_' modality_data{mode} '_' condition_names{c} ' = nAvg;'])
            % Percentage of trials
            eval(['pt_' modality_data{mode} '_' condition_names{c} ' = round((nAvg/ot_' condition_names{c} ')*100,0);'])
        end 
    end
    end 
    % Now store values in matrices
    Trial_count{p,1} = Subj_type;
    Trial_count{p,2} = ot_Standard;
    Trial_count{p,3} = ot_DeviantPitch;
    Trial_count{p,4} = ot_DeviantDuration;
    Trial_count{p,5} = st_EEG_Standard;
    Trial_count{p,6} = st_EEG_DeviantPitch;
    Trial_count{p,7} = st_EEG_DeviantDuration;
    Trial_count{p,8} = st_MEG_Standard;
    Trial_count{p,9} = st_MEG_DeviantPitch;
    Trial_count{p,10} = st_MEG_DeviantDuration;
    Trial_count{p,11} = st_BIMODAL_Standard;
    Trial_count{p,12} = st_BIMODAL_DeviantPitch;
    Trial_count{p,13} = st_BIMODAL_DeviantDuration;
    Percentage_count{p,1} = Subj_type;
    Percentage_count{p,2} = ot_Standard;
    Percentage_count{p,3} = ot_DeviantPitch;
    Percentage_count{p,4} = ot_DeviantDuration;
    Percentage_count{p,5} = [num2str(pt_EEG_Standard) ' %'];
    Percentage_count{p,6} = [num2str(pt_EEG_DeviantPitch) ' %'];
    Percentage_count{p,7} = [num2str(pt_EEG_DeviantDuration) ' %'];
    Percentage_count{p,8} = [num2str(pt_MEG_Standard) ' %'];
    Percentage_count{p,9} = [num2str(pt_MEG_DeviantPitch) ' %'];
    Percentage_count{p,10} = [num2str(pt_MEG_DeviantDuration) ' %'];
    Percentage_count{p,11} = [num2str(pt_BIMODAL_Standard) ' %'];
    Percentage_count{p,12} = [num2str(pt_BIMODAL_DeviantPitch) ' %'];
    Percentage_count{p,13} = [num2str(pt_BIMODAL_DeviantDuration) ' %'];
    % Modify the subject_array to specify usability of data and insuff cov trials
    Percentage_average_EEG = round(mean([pt_EEG_Standard,pt_EEG_DeviantPitch,pt_EEG_DeviantDuration]),0); 
    Percentage_average_MEG = round(mean([pt_MEG_Standard,pt_MEG_DeviantPitch,pt_MEG_DeviantDuration]),0); 
    % Determine now if EEG/MEG is bad based on criteria and numbers at hand
    if (Percentage_average_EEG < crit_percent) || (st_EEG_Standard < crit_sweeps) || (st_EEG_DeviantPitch < crit_sweeps) || (st_EEG_DeviantDuration < crit_sweeps)
        if ~strcmp(subject_array{pos_subj,5},'exception_EEG')
            subject_array{pos_subj,5} = 'bad_EEG';
        end
    end
    if (Percentage_average_MEG < crit_percent) || (st_MEG_Standard < crit_sweeps) || (st_MEG_DeviantPitch < crit_sweeps) || (st_MEG_DeviantDuration < crit_sweeps)
        if ~strcmp(subject_array{pos_subj,6},'exception_MEG')
            subject_array{pos_subj,6} = 'bad_MEG';
        end
    end

    % Determine if covariance matrices are ok depending on modality
    % It may be that for MEG alone we do have enough trials to compute a
    % good covariance matrix, but not for bimodal, so this information is
    % useful in case we decide to compute sources only with MEG or BIMODAL
    total_surviving_sweeps_EEG = st_EEG_Standard + st_EEG_DeviantPitch + st_EEG_DeviantDuration; 
    if total_surviving_sweeps_EEG < 37
        % Only if we manually did not consider this an exception
        if ~strcmp(subject_array{pos_subj,7},'exception_EEG')
            subject_array{pos_subj,7} = 'insuf_cov_EEG';
        end
    end
    total_surviving_sweeps_MEG = st_MEG_Standard + st_MEG_DeviantPitch + st_MEG_DeviantDuration; 
    if total_surviving_sweeps_MEG < 940 % minimum with current conditions
        % Only if we manually did not consider this an exception
        if ~strcmp(subject_array{pos_subj,8},'exception_MEG')
            subject_array{pos_subj,8} = 'insuf_cov_MEG';
        end
    end
    total_surviving_sweeps_bimodal = st_BIMODAL_Standard + st_BIMODAL_DeviantPitch + st_BIMODAL_DeviantDuration; 
    if total_surviving_sweeps_bimodal < 1344 % minimum with current conditions
        % Only if we manually did not consider this an exception
        if ~strcmp(subject_array{pos_subj,9},'exception_BIMODAL')
            subject_array{pos_subj,9} = 'insuf_cov_BIMODAL';
        end
    end
    save([root_dir '/subject_array.mat'],'subject_array') 
end
% Convert to tables
Trial_table = array2table(Trial_count,'RowNames',participant,'VariableNames',colNames);
Percentage_table = array2table(Percentage_count,'RowNames',participant,'VariableNames',colNames);
% Save
save([root_dir '/Events/Surviving_sweeps.mat'],'Trial_table');
save([root_dir '/Events/Percentage_sweeps.mat'],'Percentage_table');
% Write table in Excel
writetable(Trial_table, [root_dir '/Events/Surviving_sweeps.xlsx'])
writetable(Percentage_table, [root_dir '/Events/Percentage_sweeps.xlsx'])

clearvars('-except', initialVars{:});
disp 'DONE WITH OBTAINING SWEEP COUNT (Project MMN)!!!'
disp(datetime)
toc

%% GAVR scalp level (Brainstorm) 
% May want to reload the entire protocol first to fix errors

tic
disp(' ');      
disp('-------------------------');  
disp('GAVR SCALP LEVEL (Project MMN)');  
disp(datetime)
disp('-------------------------');     
disp(' '); 

% First, reload all folders for each subject
for p = 1:length(participant)    

    % This step is important as new epoch folders for MMN will be created
    % and won't be read by brainstorm unless we reload the whole subject
    % functional data
    disp(' ');      
    disp('-------------------------');
    disp(['loading participant ' participant{p}]);
    disp(datetime)
    disp(' ');
    
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub);
end

% Delete duplicates/keep the most recent in epoch folders
for p = 1:length(participant)
    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    
    % If subject does not have any of these markers, don't include it in
    % the sweep count (because it may be a bad subject, etc)
    if ~strcmp(subject_array{pos_subj,3},'needs_forward') && ~strcmp(subject_array{pos_subj,3},'needs_sources') && ~strcmp(subject_array{pos_subj,3},'sources_finished')
        continue;
    end
    
    % Only for the mismatch folders
    for c = 1:length(condition_mismatch_names)
    % Find folder
    files = dir([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c}]);
    if isempty(files)
        error(['No ' condition_mismatch_names{c} ' files for ' participant{p}]);
    end
    % Do once for each modality (EEG, MEG, BIMODAL)
    for mode = 1:length(modality_data)
    % data_210528_2057_BIMODAL_DeviantDuration_MMN.mat
    
    if strcmp(condition_mismatch_names{c},'DeviantPitch-Standard') || strcmp(condition_mismatch_names{c},'DeviantDuration-Standard')
        infolder = find(contains({files.name},[modality_data{mode}]) & (...
        endsWith({files.name}, '_MMN.mat')));
    else
        infolder = find(endsWith({files.name},[modality_data{mode} '_average.mat'])); 
    end
    if isempty(infolder)
       % In case, for instance, there is no EEG data for this subject
       warning(['No ' condition_mismatch_names{c} ' ' modality_data{mode} ' average for ' participant{p}]);
       continue; 
    elseif length(infolder) == 1 % there are no duplicates to delete
        continue;
    elseif length(infolder) > 1 % there are duplicates to delete
        for i = 1:length(infolder)-1 % all except last one in list (older)
            line = infolder(i);
            delete([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c} '/' files(line).name]); % always has this format
        end
    end
    end
    end 
end

% GAVR for each modality and condition (separately for FE/C)
for c = 1:length(condition_mismatch_names)
    for mode = 1:length(modality_data)
        for pg = 1:length(participant_group)
        sFiles = {};
        for p = 1:length(participant)
        % Identify patients/controls
        pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
        
        % If subject does not have any of these markers, don't include (because it may be a bad subject, etc)
        if ~strcmp(subject_array{pos_subj,3},'needs_forward') && ~strcmp(subject_array{pos_subj,3},'needs_sources') && ~strcmp(subject_array{pos_subj,3},'sources_finished')
            continue;
        end
        
        if ~strcmp(subject_array{pos_subj,2},participant_group{pg})
            continue; % so only include participants that correspond to the group
        end
        
        % Check if we should include it
        subject_row = {subject_array{pos_subj,:}}; 
        subject_row = subject_row(~cellfun('isempty', subject_row')); % to avoid empty cells in contains
        bad_signal = find(contains(subject_row,['bad_' modality_data{mode}]));
        if strcmp(modality_data{mode},'BIMODAL') % In case of bimodal, check if any EEG or MEG is bad
            bad_signal_EEG = find(contains(subject_row,'bad_EEG'));
            bad_signal_MEG = find(contains(subject_row,'bad_MEG'));
            % If either one of these is bad, mark bimodal as bad
            if ~isempty(bad_signal_EEG) || ~isempty(bad_signal_MEG)
                bad_signal = 1;
            end
        end
        if ~isempty(bad_signal)
            warning([participant{p} ' not included in ' modality_data{mode} ' scalp average']);
            continue; % to next subject
        end     
        
        files = dir([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c}]);
        if isempty(files)
            % If both EEG and MEG were bad there will be no mismatch folder at all
            warning(['No ' condition_mismatch_names{c} ' files for ' participant{p}]);
            continue;
        end
        if strcmp(condition_mismatch_names{c},'DeviantPitch-Standard') || strcmp(condition_mismatch_names{c},'DeviantDuration-Standard')
            infolder = find(contains({files.name},modality_data{mode}) & endsWith({files.name},'_MMN.mat'));
        else
            infolder = find(endsWith({files.name},[modality_data{mode} '_average.mat'])); 
        end
        if isempty(infolder)
           % In case, for instance, there is no EEG data for this subject
           warning(['No ' condition_mismatch_names{c} ' ' modality_data{mode} ' average for ' participant{p}]);
           continue;
        end  
        if length(infolder) > 1
            error(['More than one ' condition_mismatch_names{c} ' ' modality_data{mode} ' average for ' participant{p}]);
        end
        sFiles{p} = [participant{p} '/' condition_mismatch_names{c} '/' files(infolder).name];
        end
        
        sFiles = sFiles(~cellfun('isempty', sFiles')); % to avoid empty cells
        if isempty(sFiles)
            error(['No files to perform GAVR for ' condition_mismatch_names{c} ' ' modality_data{mode}]);
        end
        
        gavr_n = num2str(length(sFiles));
        
        disp(' ');
        disp('-------------------------');  
        disp(['GAVR scalp data for ' condition_mismatch_names{c} ' ' modality_data{mode} ' ' participant_group{pg}]);
        disp(datetime)
        disp(' ');

        % If stated, find and delete any previous GAVR SENSOR data
        if delete_previous_file == 1
            % check if there is already GAVR source in Group analysis folder
            folders_delete = dir([root_dir_bs '/data/Group_analysis/' condition_mismatch_names{c}]);
            % results_average_200520_2029_GAVR_Source_MEG_Regul_MLR_Quietest
            infolder_delete = find(contains({folders_delete.name},'GAVR_SENSOR_')...
            & endsWith({folders_delete.name}, [modality_data{mode} '_' condition_mismatch_names{c} '_' participant_group{pg} '_n' gavr_n '.mat']));
            if ~isempty(infolder_delete) % file exists, therefore delete it
               delete([root_dir_bs '/data/Group_analysis/' condition_mismatch_names{c} '/' folders_delete(infolder_delete).name]);
            end
        end

        % Average using default function (because we are using vertices, not channels)
        sFiles = bst_process('CallProcess', 'process_average', sFiles, [], ...
            'avgtype',         1, ...  % Everything
            'avg_func',        1, ...  % Arithmetic average:  mean(x)
            'weighted',        0, ...
            'scalenormalized', 0);

        % USE OPTION TO NORMALIZE HERE???

        % Process: Add tag
        sFiles = bst_process('CallProcess', 'process_add_tag', sFiles, [], ...
            'tag',           ['GAVR_SENSOR_' modality_data{mode} '_' condition_mismatch_names{c} '_' participant_group{pg} '_n' gavr_n], ...
            'output',        2);  % Add to file name (1 to add a tag)

        % Process: Set name
        sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
            'tag',           ['GAVR_SENSOR_' modality_data{mode} '_' condition_mismatch_names{c} '_' participant_group{pg} '_n' gavr_n], ...
            'isindex',       1);
        end
    end
end

clearvars('-except', initialVars{:});
disp 'DONE WITH GAVR SCALP LEVEL (Project MMN)!!!'
disp(datetime)
toc

%% Extract scalp data out of brainstorm and GAVR outside

tic
disp(' ');      
disp('-------------------------');  
disp('EXTRACTING SCALP VALUES OUT OF BRAINSTORM (Project MMN)');  
disp(datetime)
disp('-------------------------');     
disp(' '); 

% In this case, modality data will always be EEG and MEG only, since it's scalp
modality_data = {'EEG','MEG'};

% Select which channels/clusters are going to be used from matrix
load([root_dir '/Scripts/Areas_channels.mat']);
% EEG channels
if size(choice_channel_EEG,2) > 1 % cluster was selected
    channel_pos_EEG = [];
    for i = 1:size(choice_channel_EEG,2) % Find all the indices
        channel_pos_EEG(i) = find(strcmp({Areas.Name}, choice_channel_EEG{i})==1);
    end
else % single channel
   channel_pos_EEG = find(strcmp({Areas.Name}, choice_channel_EEG)==1); 
end
channel_pos_EEG = channel_pos_EEG - 314; % to adjust for 62 channels
% MEG sensors
if size(choice_channel_MEG,2) > 1 % cluster was selected
    channel_pos_MEG = [];
    for i = 1:size(choice_channel_MEG,2) % Find all the indices
        channel_pos_MEG(i) = find(strcmp({Areas.Name}, choice_channel_MEG{i})==1);
    end
else % single sensor
   channel_pos_MEG = find(strcmp({Areas.Name}, choice_channel_MEG)==1); 
end

% Extract individual responses of each participant and condition
for p = 1:length(participant)
    for mode = 1:length(modality_data)
        pos_subj = find(strcmp({subject_array{:,1}},participant{p}));

        % If subject does not have any of these markers, don't include it in
        % the sweep count (because it may be a bad subject, etc)
        if ~strcmp(subject_array{pos_subj,3},'needs_forward') && ~strcmp(subject_array{pos_subj,3},'needs_sources') && ~strcmp(subject_array{pos_subj,3},'sources_finished')
            continue;
        end
        
        % Check if we should include it
        subject_row = {subject_array{pos_subj,:}}; 
        subject_row = subject_row(~cellfun('isempty', subject_row')); % to avoid empty cells in contains
        bad_signal = find(contains(subject_row,['bad_' modality_data{mode}]));
        if ~isempty(bad_signal)
            warning([participant{p} ' scalp data not extracted for ' modality_data{mode}]);
            continue; % to next subject
        end    

        if ~exist([root_dir '/Scalp/' participant{p}], 'dir')
            mkdir([root_dir '/Scalp/'], [participant{p}]);
        end
        
        for c = 1:length(condition_mismatch_names)  
            
        files = dir([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c}]);
        if isempty(files)
            % If both EEG and MEG were bad there will be no mismatch folder at all
            warning(['No ' condition_mismatch_names{c} ' files for ' participant{p}]);
            continue;
        end
        if strcmp(condition_mismatch_names{c},'DeviantPitch-Standard') || strcmp(condition_mismatch_names{c},'DeviantDuration-Standard')
            infolder_average = find(contains({files.name},modality_data{mode}) & endsWith({files.name},'_MMN.mat'));
        else
            infolder = find(endsWith({files.name},[modality_data{mode} '_average.mat'])); 
        end
        if isempty(infolder)
           % In case, for instance, there is no EEG data for this subject
           warning(['No ' condition_mismatch_names{c} ' ' modality_data{mode} ' average for ' participant{p}]);
           continue;
        end  
        if length(infolder) > 1
            error(['More than one ' condition_mismatch_names{c} ' ' modality_data{mode} ' average for ' participant{p}]);
        end
        file_name = [root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c} '/' files(infolder).name];
        channel_name = [root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c} '/channel_vectorview306_acc1.mat'];
        
        % This way so that we can preserve the structure for the plotting script
        eval(['load ' file_name])
        % Load also channel file to check indices (since some subjects only have MEG channels
        eval(['load ' channel_name])
        chan_index = find(contains({Channel.Type},modality_data{mode}));
        F = F(chan_index,:);
        save([root_dir '/Scalp/' participant{p} '/' modality_data{mode} '_' condition_mismatch_names{c} '.mat'], 'F');     
        end
    end
end

% GAVR outside brainstorm (and obtain Matrix)
for mode = 1:length(modality_data)
    for c = 1:length(condition_mismatch_names)
        for pg = 1:length(participant_group)   
        file_matrix = {};
        for p = 1:length(participant)
            pos_subj = find(strcmp({subject_array{:,1}},participant{p}));

            % If subject does not have any of these markers, don't include it in
            % the sweep count (because it may be a bad subject, etc)
            if ~strcmp(subject_array{pos_subj,3},'needs_forward') && ~strcmp(subject_array{pos_subj,3},'needs_sources') && ~strcmp(subject_array{pos_subj,3},'sources_finished')
                continue;
            end
            
            % Identify patients/controls
            if ~strcmp(subject_array{pos_subj,2},participant_group{pg})
                continue; % so only include participants that correspond to the group
            end
            
            % Again, check if we should include it
            subject_row = {subject_array{pos_subj,:}}; 
            subject_row = subject_row(~cellfun('isempty', subject_row')); % to avoid empty cells in contains
            bad_signal = find(contains(subject_row,['bad_' modality_data{mode}]));
            if ~isempty(bad_signal)
                warning([participant{p} ' scalp data not averaged for ' modality_data{mode}]);
                continue; % to next subject
            end  
            
            load([root_dir '/Scalp/' participant{p} '/' modality_data{mode} '_' condition_mismatch_names{c} '.mat']);
            if strcmp(modality_data{mode},'EEG')
                F = F*1e6; % Transform from Volts to microVolts
            elseif strcmp(modality_data{mode},'MEG') % MEG
                F = F*1e15; % Transform from Tesla to Femtotesla
            end
            if size(F,1) > 306 && strcmp(modality_data{mode},'MEG')
                disp(['WARNING!! MEG DATA HAS MORE THAN 306 CHANNELS in ' participant{p}]);
                F(307:end,:)=[];
            elseif size(F,1) > 62 && strcmp(modality_data{mode},'EEG')
                disp(['WARNING!! EEG DATA HAS MORE THAN 62 CHANNELS in ' participant{p}]);                   
                F(63:end,:)=[];
            end
            file_matrix{p} = F;
        end
        
        % Delete empty subjects
        file_matrix = file_matrix(~cellfun('isempty', file_matrix')); % to avoid empty cells in contains
        % Generate files
        Matrix = cat(3,file_matrix{:});
        gavr = mean(Matrix,3);
        STD = squeeze(std(Matrix,0,3));
        num_valid_subjects = size(Matrix,3);
        STD_ERR = STD/sqrt(num_valid_subjects);

        % Save files
        if ~exist([root_dir '/Scalp/' gavr_name], 'dir')
            mkdir([root_dir '/Scalp/' gavr_name], 'gavr');
            mkdir([root_dir '/Scalp/' gavr_name], 'std_dev');
            mkdir([root_dir '/Scalp/' gavr_name], 'std_err');
        end
        F = gavr;
        eval(['save(''' root_dir '/Scalp/' gavr_name '/gavr/' modality_data{mode} '_' condition_mismatch_names{c} '_' participant_group{pg} ''', ''F'');'])
        F = STD;
        eval(['save(''' root_dir '/Scalp/' gavr_name '/std_dev/' modality_data{mode} '_' condition_mismatch_names{c} '_' participant_group{pg} ''', ''F'');'])
        F = STD_ERR;
        eval(['save(''' root_dir '/Scalp/' gavr_name  '/std_err/' modality_data{mode} '_' condition_mismatch_names{c} '_' participant_group{pg} ''', ''F'');'])
        % Save the matrix for peak average extraction too
        eval(['save(''' root_dir '/Scalp/' gavr_name '/Matrix_' modality_data{mode} '_' condition_mismatch_names{c} '_' participant_group{pg} ''', ''Matrix'');'])  
       end
    end      
end

% Obtain peak matrices
for mode = 1:length(modality_data)
    eval(['chan_num = channel_pos_' modality_data{mode} ';'])
for pg = 1:length(participant_group)   

% Define time samples, which we are gonna need
plot_baseline = -50;
plot_post = 300;
s_r = 1000;
time_samples=linspace(plot_baseline,plot_post,((((plot_baseline*(-1)) + plot_post)/1000)*s_r) +1);

% Obtain amplitude values from Matrix, selecting channel (conditions)
for c = 4:length(condition_mismatch_names) % ONLY MMN CONDITIONS FOR NOW
    % Load matrix of data
    load ([root_dir '/Scalp/' gavr_name '/Matrix_' modality_data{mode} '_' condition_mismatch_names{c} '_' participant_group{pg} '.mat'])
    % Ensure correct number of channels
    if strcmp(modality_data{mode}, 'EEG')
        if size(Matrix,1) ~= 62
            error('watch out, matrix does not have 62 channels')
        end
    elseif strcmp(modality_data{mode}, 'MEG')
        if size(Matrix,1) ~= 306
            error('watch out, matrix does not have 306 channels')
        end
    end

   % Retrieve values from channels
   Amplitudes = Matrix(chan_num,:,:);
   % if it's a cluster of channels, average first
   if size(Amplitudes,1) ~= 1
       Amplitudes = squeeze(mean(Amplitudes,1));
   else
       Amplitudes = squeeze(Amplitudes);
   end
  
   % Now before stats, we have to cut this matrix (MLR or LLR) into one
   % with 2 or 5 values instead of e.g. 276 for each subject (which will
   % correspond to N1 and P2, or to P0, Na, Pa, Nb, P50)
   for pe = 1:length(Peaks_to_extract)
       if strcmp(modality_data{mode},'EEG')
           eval(['time_window_scalp_Pitch_' Peaks_to_extract{pe} ' = time_window_scalp_Pitch_' Peaks_to_extract{pe} '_EEG;'])
           eval(['time_window_scalp_Dur_' Peaks_to_extract{pe} ' = time_window_scalp_Dur_' Peaks_to_extract{pe} '_EEG;'])
       elseif strcmp(modality_data{mode},'MEG')
           eval(['time_window_scalp_Pitch_' Peaks_to_extract{pe} ' = time_window_scalp_Pitch_' Peaks_to_extract{pe} '_MEG;'])
           eval(['time_window_scalp_Dur_' Peaks_to_extract{pe} ' = time_window_scalp_Dur_' Peaks_to_extract{pe} '_MEG;'])
       end
   % Find the indices of the start and end of time window from time_samples variable
   if strcmp(condition_mismatch_names{c},'DeviantPitch-Standard')
       eval(['time_window = time_window_scalp_Pitch_' Peaks_to_extract{pe} ';'])
   elseif strcmp(condition_mismatch_names{c},'DeviantDuration-Standard')
       eval(['time_window = time_window_scalp_Dur_' Peaks_to_extract{pe} ';'])
   end
   [~,closestIndex] = min(abs(time_samples-time_window(1)));
   init_time = closestIndex;
   [~,closestIndex] = min(abs(time_samples-time_window(2)));
   end_time = closestIndex;
   F = mean(Amplitudes(init_time:end_time,:),1);
   
   if ~exist([root_dir '/Scalp/' gavr_name '/Peaks'],'dir')
     mkdir([root_dir '/Scalp/' gavr_name '/'], 'Peaks');
   end
   
   % size(choice_channel_EEG,2) > 1
   if length(chan_num) >1 % Group of electrodes/sensors
        save([root_dir '/Scalp/' gavr_name '/Peaks/' modality_data{mode} '_' condition_mismatch_names{c} '_' Peaks_to_extract{pe} '_' participant_group{pg} '_cluster.mat'],'F');
   elseif length(chan_num) == 1
       if strcmp(modality_data{mode},'EEG')
           save([root_dir '/Scalp/' gavr_name '/Peaks/' modality_data{mode} '_' condition_mismatch_names{c} '_' Peaks_to_extract{pe} '_' participant_group{pg} '_' choice_channel_EEG{:} '.mat'],'F');
       else
           save([root_dir '/Scalp/' gavr_name '/Peaks/' modality_data{mode} '_' condition_mismatch_names{c} '_' Peaks_to_extract{pe} '_' participant_group{pg} '_' choice_channel_MEG{:} '.mat'],'F');
       end
   end
   end
end
end
end

% SPSS Matrix
for mode = 1:length(modality_data)
eval(['chan_num = choice_channel_' modality_data{mode} ';'])
temporal = []; iter = 1;
for pg = 1:length(participant_group)
for c = 4:length(condition_mismatch_names) % ONLY MMN CONDITIONS FOR NOW
for pe = 1:length(Peaks_to_extract)
    if size(chan_num,2) > 1 % Cluster
        load ([root_dir '/Scalp/' gavr_name '/Peaks/' modality_data{mode} '_' condition_mismatch_names{c} '_' Peaks_to_extract{pe} '_' participant_group{pg} '_cluster.mat']);
    else % Single channel
        load ([root_dir '/Scalp/' gavr_name '/Peaks/' modality_data{mode} '_' condition_mismatch_names{c} '_' Peaks_to_extract{pe} '_' participant_group{pg} '_' chan_num{:} '.mat']);
    end
    temporal(1:length(F),iter) = F;
    iter = iter + 1;
end
end
end

 if ~exist([root_dir '/Scalp/' gavr_name '/Statistics'],'dir')
     mkdir([root_dir '/Scalp/' gavr_name '/'], 'Statistics');
 end
 
SPSS_matrix = [header_Project_MMN; num2cell(temporal)];
if size(chan_num,2) > 1 % Cluster
    save([root_dir '/Scalp/' gavr_name '/Statistics/SPSS_Matrix_' modality_data{mode} '_cluster.mat'],'SPSS_matrix');
    xlswrite([root_dir '/Scalp/' gavr_name '/Statistics/SPSS_Matrix_' modality_data{mode} '_cluster.xlsx'],SPSS_matrix);
else % Single channel
    save([root_dir '/Scalp/' gavr_name '/Statistics/SPSS_Matrix_' modality_data{mode} '_' chan_num{:} '.mat'],'SPSS_matrix');
    xlswrite([root_dir '/Scalp/' gavr_name '/Statistics/SPSS_Matrix_' modality_data{mode} '_' chan_num{:} '.xlsx'],SPSS_matrix);
end
end

% Reset modality data to its original self
modality_data = {'EEG','MEG','BIMODAL'};

clearvars('-except', initialVars{:});
disp 'DONE WITH EXTRACTING SCALP VALUES OUT OF BRAINSTORM (Project MMN)!!!'
disp(datetime)
toc

%% (!) Ensure anatomy (FEM) and corregistration/EEG projection is correct before moving forward (!)

% Project EEG electrodes to surface (better do it manually upon checking corregistrations)

tic
disp(' ');      
disp('-------------------------');  
disp('PROJECTING EEG ELECTRODES TO SURFACE (Project MMN)');  
disp(datetime)
disp('-------------------------');     
disp(' '); 

% To ensure we keep track of which ones are projected
no_projection = 0;

for p = 1:length(participant)
    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    % Only if all of these conditions are met we can carry on
    if strcmp(subject_array{pos_subj,10},'corregistration_done') && ~strcmp(subject_array{pos_subj,13},'EEG_projected')
    
    for c = 1:length(condition_mismatch_names)
        [~, iStudies] = bst_get('StudyWithCondition', [participant{p} filesep condition_mismatch_names{c}]);
        
        disp(' ');      
        disp('-------------------------');
        disp(['loading epoch folders for participant ' participant{p} ' ' condition_mismatch_names{c}]);
        disp(datetime)
        disp(' ');
        
        db_reload_studies(iStudies);

        % Get name of any brainstorm file in the folder (from which channel file will be retrieved automatically by function)
        dir_sweeps = dir([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c}]);
        sweeps_list = find(startsWith({dir_sweeps.name},'data_'));
        position = sweeps_list(1); 
        sFiles = dir([participant{p} '/' condition_mismatch_names{c} '/' dir_sweeps(position).name]);
        
        disp(' ');      
        disp('-------------------------');
        disp(['projecting EEG electrodes to surface for ' participant{p} ' ' condition_mismatch_names{c}]);
        disp(datetime)
        disp(' ');
        
        try
            % Project electrodes to surface
            sFiles = bst_process('CallProcess', 'process_channel_project', sFiles, []);
        catch
            no_projection = 1;
        end
        
    end   
       
    if no_projection == 0
        % If successful, update subject_array for this subject
        subject_array{pos_subj,13} = 'EEG_projected';
        save([root_dir '/subject_array.mat'],'subject_array')  
    else % at least one condition was not projected, reset for next subj.
        no_projection = 0;
    end
    end
end

clearvars('-except', initialVars{:});
disp 'DONE WITH PROJECTING EEG ELECTRODES TO SURFACE (Project MMN)!!!'
disp(datetime)
toc

%% Compute forward models SURFACE and/or VOLUME (ONLY CNRL 98)
% (since this is what takes longest and requires QC)

tic
disp(' ');      
disp('-------------------------');  
disp('COMPUTING FORWARD MODELS (Project MMN)');  
disp(datetime)
disp('-------------------------');     
disp(' '); 

% Change modality_data list depending on the choice of sources
if source_analysis == 1
    modality_data = {'EEG'};
elseif source_analysis == 2
    modality_data = {'MEG'};
elseif source_analysis == 3
    modality_data = {'BIMODAL'};
elseif source_analysis == 4
    modality_data = {'EEG','MEG','BIMODAL'};
end

No_head_model_created = {};
count_no_head = 1;

% Reload Group analysis folder (important for volume forward model)
prot_subs = bst_get('ProtocolSubjects');
current_sub = find(strcmp({prot_subs.Subject.Name}, 'Group_analysis'));
db_reload_conditions(current_sub);

for p = 1:length(participant)
    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    % Only if all of these conditions are met we can carry on with sources for that subject
    if strcmp(subject_array{pos_subj,3},'needs_forward') && strcmp(subject_array{pos_subj,10},'corregistration_done')...
          && strcmp(subject_array{pos_subj,11},'FEM_anatomy_ready') && ~strcmp(subject_array{pos_subj,12},'No_HCP')...
          && strcmp(subject_array{pos_subj,13},'EEG_projected')
        
    % To check later if we can perform sources
    subject_row = {subject_array{pos_subj,:}}; 
    subject_row = subject_row(~cellfun('isempty', subject_row')); % to avoid empty cells in contains
    
    % Reload subject epoch folders first
    disp(' ');      
    disp('-------------------------');
    disp(['loading epoch folders for participant ' participant{p}]);
    disp(datetime)
    disp(' ');
    
    for c = 1:length(condition_mismatch_names)
        [~, iStudies] = bst_get('StudyWithCondition', [participant{p} filesep condition_mismatch_names{c}]);
        db_reload_studies(iStudies);
    end   
    
    % Reload anatomy folder of that subject too
    disp(' ');      
    disp('-------------------------');
    disp(['loading anatomy folder for participant ' participant{p}]);
    disp(datetime)
    disp(' ');
    
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_subjects(current_sub);
    
    % Do once for each modality (EEG, MEG, BIMODAL)
    for mode = 1:length(modality_data)
    % Check if we can perform cov or sources based on signal quality
    % If not, we don't need to bother computing the forward model
    bad_signal = find(contains(subject_row,['bad_' modality_data{mode}]));
    if strcmp(modality_data{mode},'BIMODAL') % In case of bimodal, check if any EEG or MEG is bad
        bad_signal_EEG = find(contains(subject_row,'bad_EEG'));
        bad_signal_MEG = find(contains(subject_row,'bad_MEG'));
        % If either one of these is bad, mark bimodal as bad
        if ~isempty(bad_signal_EEG) || ~isempty(bad_signal_MEG)
            bad_signal = 1;
        end
    end
    if ~isempty(bad_signal)
        warning([modality_data{mode} ' forward model not computed for ' participant{p}]);
        continue; % to next modality
    end
    % Check if we can perform cov or sources based on covariance sweeps
    % If not, we don't need to bother computing the forward model
    no_cov_sweeps = find(contains(subject_row,['insuf_cov_' modality_data{mode}]));
    if ~isempty(no_cov_sweeps)
        warning([modality_data{mode} ' forward model not computed for ' participant{p}]);
        continue; % to next modality
    end
    
    % Now, if we reach this point, we can compute forward models overwriting any previous one
    % Determine files to compute sources from (standard file)
    files = dir([root_dir_bs '/data/' participant{p} '/' condition_names{1}]);
    if isempty(files)
        error(['No ' condition_names{1} ' files for ' participant{p}]);
    end
    % It doesn't really matter if it's EEG or MEG for the purpose of the
    % forward, but since we are looping through that let's keep this
    infolder = find(endsWith({files.name},[modality_data{mode} '_average.mat'])); 
    if isempty(infolder)
       % In case, for instance, there is no EEG data for this subject
       error(['No ' condition_names{1} ' ' modality_data{mode} ' average for ' participant{p}]);
    end  
    if length(infolder) > 1
        error(['More than one ' condition_names{1} ' ' modality_data{mode} ' average for ' participant{p}]);
    end
    sFiles = [participant{p} '/' condition_names{1} '/' files(infolder).name];
    for se = 1:length(source_space)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EEG SURF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(source_space{se},'Surf') && strcmp(modality_data{mode},'EEG')
        
        % FEM DUNNEURO SIMNIBS     
        disp(' ');
        disp('-------------------------');  
        disp(['Computing forward model (' source_space{se} ' ' modality_data{mode} ') for ' participant{p}]);
        disp(datetime)
        disp(' ');   

        % If stated, find and delete any previous head model (overwrite it)    
        if delete_previous_file == 1
            % check if there is already a headmodel OF THIS TYPE in the folder
            if isfile([root_dir_bs '/data/' participant{p} '/' condition_names{1} '/headmodel_surf_duneuro.mat'])
                delete([root_dir_bs '/data/' participant{p} '/' condition_names{1} '/headmodel_surf_duneuro.mat'])
            end            
        end
        
        try
        % FEM FROM DUNNEURO SIMNIBS
        sFiles_EEG_Surf = bst_process('CallProcess', 'process_headmodel', sFiles, [], ...
            'Comment',     '', ...
            'sourcespace', 1, ...  % Cortex surface
            'meg',         1, ...  % 
            'eeg',         4, ...  % DUNEuro FEM
            'ecog',        1, ...  % 
            'seeg',        1, ...  % 
            'duneuro',     struct(...
                 'FemCond',             [0.14, 0.33, 1.79, 0.008, 0.43], ...
                 'FemSelect',           [1, 1, 1, 1, 1], ...
                 'UseTensor',           0, ...
                 'Isotropic',           1, ...
                 'SrcShrink',           0, ...
                 'SrcForceInGM',        1, ...
                 'FemType',             'fitted', ...
                 'SolverType',          'cg', ...
                 'GeometryAdapted',     0, ...
                 'Tolerance',           1e-08, ...
                 'ElecType',            'normal', ...
                 'MegIntorderadd',      0, ...
                 'MegType',             'physical', ...
                 'SolvSolverType',      'cg', ...
                 'SolvPrecond',         'amg', ...
                 'SolvSmootherType',    'ssor', ...
                 'SolvIntorderadd',     0, ...
                 'DgSmootherType',      'ssor', ...
                 'DgScheme',            'sipg', ...
                 'DgPenalty',           20, ...
                 'DgEdgeNormType',      'houston', ...
                 'DgWeights',           1, ...
                 'DgReduction',         1, ...
                 'SolPostProcess',      1, ...
                 'SolSubstractMean',    0, ...
                 'SolSolverReduction',  1e-10, ...
                 'SrcModel',            'venant', ...
                 'SrcIntorderadd',      0, ...
                 'SrcIntorderadd_lb',   2, ...
                 'SrcNbMoments',        3, ...
                 'SrcRefLen',           20, ...
                 'SrcWeightExp',        1, ...
                 'SrcRelaxFactor',      6, ...
                 'SrcMixedMoments',     1, ...
                 'SrcRestrict',         1, ...
                 'SrcInit',             'closest_vertex', ...
                 'BstSaveTransfer',     0, ...
                 'BstEegTransferFile',  'eeg_transfer.dat', ...
                 'BstMegTransferFile',  'meg_transfer.dat', ...
                 'BstEegLfFile',        'eeg_lf.dat', ...
                 'BstMegLfFile',        'meg_lf.dat', ...
                 'UseIntegrationPoint', 1, ...
                 'EnableCacheMemory',   0, ...
                 'MegPerBlockOfSensor', 0), ...
            'channelfile', '');
        catch
            No_head_model_created{count_no_head,1} = [participant{p} ' has no' source_space{se} ' ' modality_data{mode} ' forward model'];
            count_no_head = count_no_head +1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MEG SURF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif strcmp(source_space{se},'Surf') && strcmp(modality_data{mode},'MEG')
            
        % FEM DUNNEURO SIMNIBS     
        disp(' ');
        disp('-------------------------');  
        disp(['Computing forward model (' source_space{se} ' ' modality_data{mode} ') for ' participant{p}]);
        disp(datetime)
        disp(' ');   

        % If stated, find and delete any previous head model (overwrite it)  
        if delete_previous_file == 1
            % check if there is already a headmodel OF THIS TYPE in the folder
            if isfile([root_dir_bs '/data/' participant{p} '/' condition_names{1} '/headmodel_surf_os_meg.mat'])
                delete([root_dir_bs '/data/' participant{p} '/' condition_names{1} '/headmodel_surf_os_meg.mat'])
            end            
        end
        
        try
        % Process: Compute head model
        sFiles_MEG_Surf = bst_process('CallProcess', 'process_headmodel', sFiles, [], ...
            'Comment',     '', ...
            'sourcespace', 1, ...  % Cortex surface
            'volumegrid',  struct(...
                 'Method',        'isotropic', ...
                 'nLayers',       17, ...
                 'Reduction',     3, ...
                 'nVerticesInit', 4000, ...
                 'Resolution',    0.005, ...
                 'FileName',      []), ...
            'meg',         3, ...  % Overlapping spheres
            'eeg',         1, ...  % 
            'ecog',        1, ...  % 
            'seeg',        1, ...  % 
            'openmeeg',    struct(...
                 'BemFiles',     {{}}, ...
                 'BemNames',     {{'Scalp', 'Skull', 'Brain'}}, ...
                 'BemCond',      [1, 0.0125, 1], ...
                 'BemSelect',    [1, 1, 1], ...
                 'isAdjoint',    0, ...
                 'isAdaptative', 1, ...
                 'isSplit',      0, ...
                 'SplitLength',  4000));
        catch
            No_head_model_created{count_no_head,1} = [participant{p} ' has no' source_space{se} ' ' modality_data{mode} ' forward model'];
            count_no_head = count_no_head +1;
        end
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BIMODAL SURF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif strcmp(source_space{se},'Surf') && strcmp(modality_data{mode},'BIMODAL')
        % We won't get here unless requirements for EEG and MEG are met
        % So if we are here it means we SHOULD have EEG and MEG forward
        % models for that subject. If we don't it must be because of a
        % problem computing them (using try catch).
        
        EEG_leadfield = [root_dir_bs '/data/' participant{p} '/' condition_names{1} '/headmodel_surf_duneuro.mat'];
        MEG_leadfield = [root_dir_bs '/data/' participant{p} '/' condition_names{1} '/headmodel_surf_os_meg.mat'];

        if ~exist(EEG_leadfield,'file') || ~exist(MEG_leadfield,'file')
            No_head_model_created{count_no_head,1} = [participant{p} ' has no' source_space{se} ' ' modality_data{mode} ' forward model'];
            count_no_head = count_no_head +1;
        else
            
            disp(' ');
            disp('-------------------------');  
            disp(['Merging EEG/MEG leadfields(' source_space{se} ') for ' participant{p}]);
            disp(datetime)
            disp(' '); 
            
            if delete_previous_file == 1
                % check if there is already a headmodel OF THIS TYPE in the folder
                if isfile([root_dir_bs '/data/' participant{p} '/' condition_names{1} '/headmodel_surf_duneuro_os_meg.mat'])
                    delete([root_dir_bs '/data/' participant{p} '/' condition_names{1} '/headmodel_surf_duneuro_os_meg.mat'])
                end            
            end
            
            % FUNCTION TO COMBINE LEADFIELDS
            merge_leadfields(EEG_leadfield, MEG_leadfield)
            % Will combine and delete the original ones
            % Outcome will be saved with same name as if it was
            % computed manually: headmodel_surf_duneuro_os_meg.mat
            
        end
  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EEG VOL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif strcmp(source_space{se},'Vol') && strcmp(modality_data{mode},'EEG')
          
        disp(' ');
        disp('-------------------------');  
        disp(['Computing forward model (' source_space{se} ' ' modality_data{mode} ') for ' participant{p}]);
        disp(datetime)
        disp(' ');   

        % If stated, find and delete any previous head model (overwrite it)    
        if delete_previous_file == 1
            % check if there is already a headmodel OF THIS TYPE in the folder
            if isfile([root_dir_bs '/data/' participant{p} '/' condition_names{1} '/headmodel_vol_duneuro.mat'])
                delete([root_dir_bs '/data/' participant{p} '/' condition_names{1} '/headmodel_vol_duneuro.mat'])
            end            
        end
        
        try
        % FEM FROM DUNNEURO SIMNIBS
        sFiles_EEG_Vol = bst_process('CallProcess', 'process_headmodel', sFiles, [], ...
            'Comment',     '', ...
            'sourcespace', 2, ...  % MRI volume
            'volumegrid',  struct(...
                 'Method',        'group', ...
                 'nLayers',       17, ...
                 'Reduction',     3, ...
                 'nVerticesInit', 4000, ...
                 'Resolution',    0.003, ...
                 'FileName',      ['Group_analysis/@default_study/' name_volume_source_grid]), ...
            'meg',         1, ...  % 
            'eeg',         4, ...  % DUNEuro FEM
            'ecog',        1, ...  % 
            'seeg',        1, ...  % 
            'duneuro',     struct(...
                 'FemCond',             [0.14, 0.33, 1.79, 0.008, 0.43], ...
                 'FemSelect',           [1, 1, 1, 1, 1], ...
                 'UseTensor',           0, ...
                 'Isotropic',           1, ...
                 'SrcShrink',           0, ...
                 'SrcForceInGM',        1, ...
                 'FemType',             'fitted', ...
                 'SolverType',          'cg', ...
                 'GeometryAdapted',     0, ...
                 'Tolerance',           1e-08, ...
                 'ElecType',            'normal', ...
                 'MegIntorderadd',      0, ...
                 'MegType',             'physical', ...
                 'SolvSolverType',      'cg', ...
                 'SolvPrecond',         'amg', ...
                 'SolvSmootherType',    'ssor', ...
                 'SolvIntorderadd',     0, ...
                 'DgSmootherType',      'ssor', ...
                 'DgScheme',            'sipg', ...
                 'DgPenalty',           20, ...
                 'DgEdgeNormType',      'houston', ...
                 'DgWeights',           1, ...
                 'DgReduction',         1, ...
                 'SolPostProcess',      1, ...
                 'SolSubstractMean',    0, ...
                 'SolSolverReduction',  1e-10, ...
                 'SrcModel',            'venant', ...
                 'SrcIntorderadd',      0, ...
                 'SrcIntorderadd_lb',   2, ...
                 'SrcNbMoments',        3, ...
                 'SrcRefLen',           20, ...
                 'SrcWeightExp',        1, ...
                 'SrcRelaxFactor',      6, ...
                 'SrcMixedMoments',     1, ...
                 'SrcRestrict',         1, ...
                 'SrcInit',             'closest_vertex', ...
                 'BstSaveTransfer',     0, ...
                 'BstEegTransferFile',  'eeg_transfer.dat', ...
                 'BstMegTransferFile',  'meg_transfer.dat', ...
                 'BstEegLfFile',        'eeg_lf.dat', ...
                 'BstMegLfFile',        'meg_lf.dat', ...
                 'UseIntegrationPoint', 1, ...
                 'EnableCacheMemory',   0, ...
                 'MegPerBlockOfSensor', 0), ...
            'channelfile', '');
        catch
            No_head_model_created{count_no_head,1} = [participant{p} ' has no' source_space{se} ' ' modality_data{mode} ' forward model'];
            count_no_head = count_no_head +1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MEG VOL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        elseif strcmp(source_space{se},'Vol') && strcmp(modality_data{mode},'MEG')
        
        disp(' ');
        disp('-------------------------');  
        disp(['Computing forward model (' source_space{se} ' ' modality_data{mode} ') for ' participant{p}]);
        disp(datetime)
        disp(' ');   

        % If stated, find and delete any previous head model (overwrite it)    
        if delete_previous_file == 1
            % check if there is already a headmodel OF THIS TYPE in the folder
            if isfile([root_dir_bs '/data/' participant{p} '/' condition_names{1} '/headmodel_vol_os_meg.mat'])
                delete([root_dir_bs '/data/' participant{p} '/' condition_names{1} '/headmodel_vol_os_meg.mat'])
            end            
        end
            
        try
        % Process: Compute head model
        sFiles_MEG_Vol = bst_process('CallProcess', 'process_headmodel', sFiles, [], ...
            'Comment',     '', ...
            'sourcespace', 2, ...  % MRI volume
            'volumegrid',  struct(...
                 'Method',        'group', ...
                 'nLayers',       17, ...
                 'Reduction',     3, ...
                 'nVerticesInit', 4000, ...
                 'Resolution',    0.003, ...
                 'FileName',      ['Group_analysis/@default_study/' name_volume_source_grid]), ...
            'meg',         3, ...  % Overlapping spheres
            'eeg',         1, ...  % 
            'ecog',        1, ...  % 
            'seeg',        1, ...  % 
            'channelfile', '');
        catch
            No_head_model_created{count_no_head,1} = [participant{p} ' has no' source_space{se} ' ' modality_data{mode} ' forward model'];
            count_no_head = count_no_head +1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BIMODAL VOL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        elseif strcmp(source_space{se},'Vol') && strcmp(modality_data{mode},'BIMODAL')
            
        EEG_leadfield = [root_dir_bs '/data/' participant{p} '/' condition_names{1} '/headmodel_vol_duneuro.mat'];
        MEG_leadfield = [root_dir_bs '/data/' participant{p} '/' condition_names{1} '/headmodel_vol_os_meg.mat'];

        if ~exist(EEG_leadfield,'file') || ~exist(MEG_leadfield,'file')
            No_head_model_created{count_no_head,1} = [participant{p} ' has no' source_space{se} ' ' modality_data{mode} ' forward model'];
            count_no_head = count_no_head +1;
        else
            
            disp(' ');
            disp('-------------------------');  
            disp(['Merging EEG/MEG leadfields(' source_space{se} ') for ' participant{p}]);
            disp(datetime)
            disp(' '); 
            
            if delete_previous_file == 1
                % check if there is already a headmodel OF THIS TYPE in the folder
                if isfile([root_dir_bs '/data/' participant{p} '/' condition_names{1} '/headmodel_vol_duneuro_os_meg.mat'])
                    delete([root_dir_bs '/data/' participant{p} '/' condition_names{1} '/headmodel_vol_duneuro_os_meg.mat'])
                end            
            end
            
            % FUNCTION TO COMBINE LEADFIELDS
            merge_leadfields(EEG_leadfield, MEG_leadfield)
            % Will combine and delete the original ones
            % Outcome will be saved with same name as if it was
            % computed manually: headmodel_surf_duneuro_os_meg.mat
            
        end
        end   
    end
    end
    % Now, copy all forward models to other conditions
    files = dir([root_dir_bs '/data/' participant{p} '/' condition_names{1}]);
    infolder = find(contains({files.name},'headmodel'));
    % If there is something to copy
    if ~isempty(infolder)
        for c = 2:length(condition_mismatch_names)

            disp(' ');      
            disp('-------------------------');  
            disp(['Copying forward models to ' condition_mismatch_names{c} ' for ' participant{p}]);
            disp(datetime)
            disp(' '); 
            
            % Theoretically it will never be more than one. If EEG or MEG
            % cannot be computed there will be no BIMODAL (therefore only
            % one computed). If they both can be computed, they will be
            % merged and only the BIMODAL will remain. But just in case.
            for i= 1:length(infolder)
                line = infolder(i);
                file_name_leadfield = [root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{1} '/' files(line).name];
                % If stated, find and delete any previous head model (overwrite it)
                if delete_previous_file == 1
                    % check if there is already a headmodel OF THIS
                    % TYPE in the folder where we are going to copy the file
                    if isfile([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c} '/' files(line).name])
                        delete([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c} '/' files(line).name])
                    end             
                end
                % If destiny folder exists
                if exist([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c} '/'],'dir')
                    copyfile(file_name_leadfield,[root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c} '/' files(line).name]);
                end
            end
        end
    end
    % If successful, update subject_array for this subject
    subject_array{pos_subj,3} = 'needs_sources';
    save([root_dir '/subject_array.mat'],'subject_array')    
    end  
end

% Reset modality_data variable to its original self
modality_data = {'EEG','MEG','BIMODAL'};
save([root_dir '/QC/No_head_model_created.mat'],'No_head_model_created');

clearvars('-except', initialVars{:});
disp 'DONE WITH COMPUTING FORWARD MODELS (Project MMN)!!!'
disp(datetime)
toc

%% Compute noise covariance and inverse solutions

tic
disp(' ');      
disp('-------------------------');  
disp('COMPUTING NOISE COVARIANCE AND SOURCES (Project MMN)');  
disp(datetime)
disp('-------------------------');     
disp(' '); 

% Change modality_data list depending on the choice of sources
if source_analysis == 1
    modality_data = {'EEG'};
elseif source_analysis == 2
    modality_data = {'MEG'};
elseif source_analysis == 3
    modality_data = {'BIMODAL'};
elseif source_analysis == 4
    modality_data = {'EEG','MEG','BIMODAL'};
end

No_inverse_solution = {};
count_no_inverse = 1;

% Reload Group analysis folder (important for volume forward model)
prot_subs = bst_get('ProtocolSubjects');
current_sub = find(strcmp({prot_subs.Subject.Name}, 'Group_analysis'));
db_reload_conditions(current_sub);

for p = 1:length(participant)
    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    % Only if all of these conditions are met, carry on
    if strcmp(subject_array{pos_subj,3},'needs_sources') && strcmp(subject_array{pos_subj,10},'corregistration_done')...
          && strcmp(subject_array{pos_subj,11},'FEM_anatomy_ready') && ~strcmp(subject_array{pos_subj,12},'No_HCP')...
          && strcmp(subject_array{pos_subj,13},'EEG_projected')
    % To check later if we can perform cov and sources for this modality
    subject_row = {subject_array{pos_subj,:}}; 
    subject_row = subject_row(~cellfun('isempty', subject_row')); % to avoid empty cells in contains
    % For bimodal vs individual covariances later
    bimodal_cov = 0;
    
    % Reload subject epoch folders first
    disp(' ');      
    disp('-------------------------');
    disp(['loading epoch folders for participant ' participant{p}]);
    disp(datetime)
    disp(' ');
    
    for c = 1:length(condition_mismatch_names)
        [~, iStudies] = bst_get('StudyWithCondition', [participant{p} filesep condition_mismatch_names{c}]);
        db_reload_studies(iStudies);
    end   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%% COVARIANCE MATRIX BIMODAL %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if compute_covariance == 1 % In case we don't want to repeat this
    % Check which covariance matrix to compute (EEG/MEG only, or BIMODAL)
    bad_signal_EEG = find(contains(subject_row,'bad_EEG'));
    bad_signal_MEG = find(contains(subject_row,'bad_MEG'));
    no_cov_sweeps_EEG = find(contains(subject_row,'insuf_cov_EEG'));
    no_cov_sweeps_MEG = find(contains(subject_row,'insuf_cov_MEG'));
    no_cov_sweeps_BIMODAL = find(contains(subject_row,'insuf_cov_BIMODAL'));
    % If both EEG and MEG are good cov and data-wise, compute BIMODAL covariance matrix
    % WILL BE USEFUL FOR THE REST OF THE CONDITIONS/SURFACES/MODALITIES
    if isempty(bad_signal_EEG) && isempty(bad_signal_MEG) && isempty(no_cov_sweeps_EEG)...
        && isempty(no_cov_sweeps_MEG) && isempty(no_cov_sweeps_BIMODAL)
    
    % COMPUTE BIMODAL NOISE COVARIANCE
    % Find sweeps used to compute this average as input
    cov_names = {}; count = 1;
    for c = 1:length(condition_names) % where real trials are
        files = dir([root_dir_bs '/data/' participant{p} '/' condition_names{c}]);
        if isempty(files)
            error(['No ' condition_names{c} ' files for ' participant{p}]);
        end
        infolder = find(endsWith({files.name},'_BIMODAL_average.mat')); 
        if isempty(infolder)
            % If it didn't exist (although it should if enough trials were present)
            error(['No ' condition_names{c} ' BIMODAL average for ' participant{p}]);
        end 
        if length(infolder) > 1
            error(['More than one BIMODAL average for ' participant{p}]);
        end
        % Load average to get the sweep names from file History
        load([root_dir_bs '/data/' participant{p} '/' condition_names{c} '/' files(infolder).name]);
        % Structure of trials
        trial_indices = find(contains(History(:,3),['- ' participant{p} '/' condition_names{c} '/data_']));
        for i = 1:length(trial_indices)
            position = trial_indices(i);
            cov_name = History{position,3};
            cov_name = cov_name(4:end); % to remove - symbol
            cov_names{count} = cov_name;
            count = count +1;
        end
    end
    cov_names = cov_names(~cellfun('isempty', cov_names')); % to avoid empty cells
    
    disp(' ');
    disp('-------------------------');  
    disp(['Computing BIMODAL covariance for ' participant{p}]);
    disp(datetime)
    disp(' '); 
    
    % Process: Compute noise covariance
    cov_names = bst_process('CallProcess', 'process_noisecov', cov_names, [], ...
        'baseline',       time_noise_cov, ...
        'datatimewindow', [0, 0], ...
        'sensortypes',    'MEG, EEG', ... % both in this case (BIMODAL)
        'target',         1, ...  % Noise covariance (covariance over baseline time window)
        'dcoffset',       1, ...  % Block by block, to avoid effects of slow shifts in data
        'identity',       0, ...
        'copycond',       0, ...
        'copysubj',       0, ...
        'copymatch',      0, ...
        'replacefile',    1);  % Replace (will always replace the previous file)
    
    % Now copy cov matrices to MMN folders
    origin_file = [root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{1} '/noisecov_full.mat'];
    % Any of the three first can do, but let's use standard one
    for c = 4:length(condition_mismatch_names) % last two folders of list 
        copyfile(origin_file,[root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c} '/noisecov_full.mat']);
        % Reload destiny folder in bs
        file_to_search = [participant{p} '/' condition_mismatch_names{c} '/channel_vectorview306_acc1.mat'];
        [sStudy, iStudy] = bst_get('AnyFile', file_to_search);
        db_reload_studies(iStudy);
    end
    
    % And indicate that bimodal was computed for this subject
    bimodal_cov = 1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%% SOURCES COMPUTATION %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for mode = 1:length(modality_data) % Do sources once for each modality (EEG, MEG, BIMODAL). If possible
    % Double check if we can perform cov or sources based on signal quality
    bad_signal = find(contains(subject_row,['bad_' modality_data{mode}]));
    if strcmp(modality_data{mode},'BIMODAL') % In case of bimodal, check if any EEG or MEG is bad
        bad_signal_EEG = find(contains(subject_row,'bad_EEG'));
        bad_signal_MEG = find(contains(subject_row,'bad_MEG'));
        % If either one of these is bad, mark bimodal as bad
        if ~isempty(bad_signal_EEG) || ~isempty(bad_signal_MEG)
            bad_signal = 1;
        end
    end
    if ~isempty(bad_signal)
        warning([modality_data{mode} ' covariance and sources not computed for ' participant{p}]);
        continue; % to next modality
    end
    % Check if we can perform cov or sources based on covariance sweeps
    no_cov_sweeps = find(contains(subject_row,['insuf_cov_' modality_data{mode}]));
    if ~isempty(no_cov_sweeps)
        warning([modality_data{mode} ' covariance and sources not computed for ' participant{p}]);
        continue; % to next modality
    end
    
    
    % IF BIMODAL IS OK WE WILL USE THAT COV MATRIX WITH EEG AND MEG
    % BUT IF IT'S NOT, COV MATRIX WILL HAVE TO BE COMPUTED INDIVIDUALLY
    % FOR EACH MODALITY HERE BEFORE SOMPUTING SOURCES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%% COVARIANCE MATRIX EEG/MEG %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if compute_covariance == 1 % In case we don't want to repeat this
    if bimodal_cov == 0 % so do only if bimodal cov matrix was not computed
    % Try to compute EEG or MEG only cov matrix (if getting here means we can)
  
    % COMPUTE EEG OR MEG ONLY NOISE COVARIANCE
    % Find sweeps used to compute this average as input
    cov_names = {}; count = 1;
    for c = 1:length(condition_names) % where real trials are
        files = dir([root_dir_bs '/data/' participant{p} '/' condition_names{c}]);
        if isempty(files)
            error(['No ' condition_names{c} ' files for ' participant{p}]);
        end
        infolder = find(endsWith({files.name},['_' modality_data{mode} '_average.mat'])); 
        if isempty(infolder)
            % If it didn't exist (although it should if enough trials were present)
            error(['No ' condition_names{c} ' ' modality_data{mode} ' average for ' participant{p}]);
        end 
        if length(infolder) > 1
            error(['More than one ' modality_data{mode} ' average for ' participant{p}]);
        end
        % Load average to get the sweep names from file History
        load([root_dir_bs '/data/' participant{p} '/' condition_names{c} '/' files(infolder).name]);
        % Structure of trials
        trial_indices = find(contains(History(:,3),['- ' participant{p} '/' condition_names{c} '/data_']));
        for i = 1:length(trial_indices)
            position = trial_indices(i);
            cov_name = History{position,3};
            cov_name = cov_name(4:end); % to remove - symbol
            cov_names{count} = cov_name;
            count = count +1;
        end
    end
    cov_names = cov_names(~cellfun('isempty', cov_names')); % to avoid empty cells
    
    disp(' ');
    disp('-------------------------');  
    disp(['Computing ' modality_data{mode} ' covariance for ' participant{p}]);
    disp(datetime)
    disp(' '); 
    
    % Process: Compute noise covariance
    cov_names = bst_process('CallProcess', 'process_noisecov', cov_names, [], ...
        'baseline',       time_noise_cov, ...
        'datatimewindow', [0, 0], ...
        'sensortypes',    modality_data{mode}, ... % Only EEG or MEG in this case
        'target',         1, ...  % Noise covariance (covariance over baseline time window)
        'dcoffset',       1, ...  % Block by block, to avoid effects of slow shifts in data
        'identity',       0, ...
        'copycond',       0, ...
        'copysubj',       0, ...
        'copymatch',      0, ...
        'replacefile',    1);  % Replace (will always replace the previous file)
    
    % Now copy cov matrices to MMN folders
    origin_file = [root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{1} '/noisecov_full.mat'];
    % Any of the three first can do, but let's use standard one
    for c = 4:length(condition_mismatch_names) % last two folders of list 
        copyfile(origin_file,[root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c} '/noisecov_full.mat']);
        % Reload destiny folder in bs
        file_to_search = [participant{p} '/' condition_mismatch_names{c} '/channel_vectorview306_acc1.mat'];
        [sStudy, iStudy] = bst_get('AnyFile', file_to_search);
        db_reload_studies(iStudy);
    end
    
    end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Now, if we reach this point, we may compute inverse solutions...
    for c = 1:length(condition_mismatch_names)
    % Determine files to compute sources from
    files = dir([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c}]);
    if isempty(files)
        error(['No ' condition_mismatch_names{c} ' files for ' participant{p}]);
    end
    % Select average accordingly
    if strcmp(condition_mismatch_names{c},'DeviantPitch-Standard') || strcmp(condition_mismatch_names{c},'DeviantDuration-Standard')
        infolder = find(contains({files.name},modality_data{mode}) & endsWith({files.name},'_MMN.mat'));
    else % normal conditions
        infolder = find(endsWith({files.name},[modality_data{mode} '_average.mat'])); 
    end
    if isempty(infolder)
       % Just in case
       warning(['No ' condition_mismatch_names{c} ' ' modality_data{mode} ' average for ' participant{p}]);
       No_inverse_solution{count_no_inverse,1} = [participant{p} ' has no ' modality_data{mode} ' inverse solutions (no average)'];
       count_no_inverse = count_no_inverse +1;
       continue; % to next condition or next modality
    end  
    if length(infolder) > 1
        warning(['More than one ' condition_mismatch_names{c} ' ' modality_data{mode} ' average for ' participant{p}]);
        No_inverse_solution{count_no_inverse,1} = [participant{p} ' has no ' modality_data{mode} ' inverse solutions (more than one average)'];
        count_no_inverse = count_no_inverse +1;
        continue; % to next condition or next modality
    end
    
    % IMPORTANT: keep the original sFiles variable alive within this loop
    sFiles = [participant{p} '/' condition_mismatch_names{c} '/' files(infolder).name];
    
    % Check if covariance matrix is computed before moving forward 
    if ~exist([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c} '/noisecov_full.mat'], 'file')
        No_inverse_solution{count_no_inverse,1} = [participant{p} ' ' condition_mismatch_names{c} ' has no ' modality_data{mode} ' inverse solutions (no cov matrix)'];
        count_no_inverse = count_no_inverse +1;
        continue; % Go to next condition/modality
    end
    
    for se = 1:length(source_space)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EEG SURF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(source_space{se},'Surf') && strcmp(modality_data{mode},'EEG')
        
        % Check if any forward model for this modality and space exists before moving forward
        if ~exist([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c} '/headmodel_surf_duneuro.mat'], 'file')...
            && ~exist([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c} '/headmodel_surf_duneuro_os_meg.mat'], 'file')
            % If BIMODAL one exists there should be no individual EEG one (merged with MEG) and that's still ok!
            No_inverse_solution{count_no_inverse,1} = [participant{p} ' ' condition_mismatch_names{c} ' has no ' modality_data{mode} ' ' source_space{se} ' inverse solutions (no forward model)'];
            count_no_inverse = count_no_inverse +1;
            continue; % Go to next condition/modality
        end    
                 
        disp(' ');
        disp('-------------------------');  
        disp(['Computing cortical solutions (' source_space{se} ' ' modality_data{mode} ') for ' participant{p}]);
        disp(datetime)
        disp(' ');   
        
        % If stated, find and delete any previous source kernel (overwrite it) before copying 
        if delete_previous_file == 1
            % check if there is already a source of THIS TYPE in the folder where we are gonna calculate it
            folders_delete = dir([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c}]);
            results_delete =  contains({folders_delete.name},'KERNEL') & endsWith({folders_delete.name}, ['Surf_' source_noise_tag '_' modality_data{mode} '.mat']); 
            infolder_delete = find(results_delete);
            if ~isempty(infolder_delete) % file exists, therefore delete it
               delete([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c} '/' folders_delete(infolder_delete).name]);
            end            
        end

        % NOW WE NEED TO SPECIFY WHICH FORWARD MODEL IS THE ONE SELECTED
        file_to_search = [participant{p} '/' condition_mismatch_names{c} '/channel_vectorview306_acc1.mat'];
        [sStudy, iStudy] = bst_get('AnyFile', file_to_search);
        % Possible ones here: surf_duneuro.mat, surf_os_meg.mat, surf_duneuro_os_meg.mat
        % Try to find the bimodal first for this source space
        surface_index = find(endsWith({sStudy.HeadModel.FileName},'surf_duneuro_os_meg.mat'), 1);
        if isempty(surface_index)
            % Try the one specific for EEG then
            surface_index = find(endsWith({sStudy.HeadModel.FileName},'surf_duneuro.mat'), 1);
        end
        % We now at this point that one of the two exists, but if more than one copy is present, show error
        if length(surface_index) > 1
            error('more than one surface head model found');
        end
        sStudy.iHeadModel = surface_index;
        bst_set('Study', iStudy, sStudy);
        % According to brainstorm developers, better not to reload
        % db_reload_studies(iStudy);

        % Finnally, compute the inverse solutions or keep track if not
        try
        % Process: Compute sources [2018]
        sFiles_source = bst_process('CallProcess', 'process_inverse_2018', sFiles, [], ...
            'output',  1, ...  % Kernel only: shared
            'inverse', struct(...
                 'Comment',        ['dSPM-unscaled: ' modality_data{mode} ' ALL'], ...
                 'InverseMethod',  'minnorm', ...
                 'InverseMeasure', 'dspm2018', ...
                 'SourceOrient',   {{'loose'}}, ...
                 'Loose',          0.4, ...
                 'UseDepth',       1, ...
                 'WeightExp',      0.5, ...
                 'WeightLimit',    10, ...
                 'NoiseMethod',    source_noise_option, ...
                 'NoiseReg',       0.1, ...
                 'SnrMethod',      'fixed', ...
                 'SnrRms',         1e-06, ...
                 'SnrFixed',       3, ...
                 'ComputeKernel',  1, ...
                 'DataTypes',      {{'EEG'}}));
        catch
            No_inverse_solution{count_no_inverse,1} = [participant{p} ' ' condition_mismatch_names{c} ' has no ' modality_data{mode} ' ' source_space{se} ' inverse solutions (computing failed)'];
            count_no_inverse = count_no_inverse +1;
            continue; % Go to next condition/modality
        end

        % Process: Change name in folder (for code)
        sFiles_name = bst_process('CallProcess', 'process_add_tag', sFiles_source, [], ...
            'tag',           ['Surf_' source_noise_tag '_' modality_data{mode}], ...
            'output',        2);  % Add to file path

        % Process: Set name to see in GUI
        sFiles_name = bst_process('CallProcess', 'process_set_comment', sFiles_name, [], ...
            'tag',           ['Surf_' source_noise_tag '_' modality_data{mode}], ...
            'isindex',       0);
        clear sFiles_source sFiles_name
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MEG SURF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif strcmp(source_space{se},'Surf') && strcmp(modality_data{mode},'MEG')  
        % Check if any forward model for this modality and space exists before moving forward
        if ~exist([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c} '/headmodel_surf_os_meg.mat'], 'file')...
            && ~exist([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c} '/headmodel_surf_duneuro_os_meg.mat'], 'file')
            % If BIMODAL one exists there should be no individual MEG one (merged with MEG) and that's still ok!
            No_inverse_solution{count_no_inverse,1} = [participant{p} ' ' condition_mismatch_names{c} ' has no ' modality_data{mode} ' ' source_space{se} ' inverse solutions (no forward model)'];
            count_no_inverse = count_no_inverse +1;
            continue; % Go to next condition/modality
        end    
                
        disp(' ');
        disp('-------------------------');  
        disp(['Computing cortical solutions (' source_space{se} ' ' modality_data{mode} ') for ' participant{p}]);
        disp(datetime)
        disp(' ');   
        
        % If stated, find and delete any previous source kernel (overwrite it) before copying 
        if delete_previous_file == 1
            % check if there is already a source of THIS TYPE in the folder where we are gonna calculate it
            folders_delete = dir([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c}]);
            results_delete =  contains({folders_delete.name},'KERNEL') & endsWith({folders_delete.name}, ['Surf_' source_noise_tag '_' modality_data{mode} '.mat']); 
            infolder_delete = find(results_delete);
            if ~isempty(infolder_delete) % file exists, therefore delete it
               delete([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c} '/' folders_delete(infolder_delete).name]);
            end            
        end

        % NOW WE NEED TO SPECIFY WHICH FORWARD MODEL IS THE ONE SELECTED
        file_to_search = [participant{p} '/' condition_mismatch_names{c} '/channel_vectorview306_acc1.mat'];
        [sStudy, iStudy] = bst_get('AnyFile', file_to_search);
        % Possible ones here: surf_duneuro.mat, surf_os_meg.mat, surf_duneuro_os_meg.mat
        % Try to find the bimodal first for this source space
        surface_index = find(endsWith({sStudy.HeadModel.FileName},'surf_duneuro_os_meg.mat'), 1);
        if isempty(surface_index)
            % Try the one specific for MEG then
            surface_index = find(endsWith({sStudy.HeadModel.FileName},'surf_os_meg.mat'), 1);
        end
        % We now at this point that one of the two exists, but if more than one copy is present, show error
        if length(surface_index) > 1
            error('more than one surface head model found');
        end
        sStudy.iHeadModel = surface_index;
        bst_set('Study', iStudy, sStudy);
        % According to brainstorm developers, better not to reload
        % db_reload_studies(iStudy);

        % Finnally, compute the inverse solutions or keep track if not
        try
        % Process: Compute sources [2018]
        sFiles_source = bst_process('CallProcess', 'process_inverse_2018', sFiles, [], ...
            'output',  1, ...  % Kernel only: shared
            'inverse', struct(...
                 'Comment',        ['dSPM-unscaled: ' modality_data{mode} ' ALL'], ...
                 'InverseMethod',  'minnorm', ...
                 'InverseMeasure', 'dspm2018', ...
                 'SourceOrient',   {{'loose'}}, ...
                 'Loose',          0.4, ...
                 'UseDepth',       1, ...
                 'WeightExp',      0.5, ...
                 'WeightLimit',    10, ...
                 'NoiseMethod',    source_noise_option, ...
                 'NoiseReg',       0.1, ...
                 'SnrMethod',      'fixed', ...
                 'SnrRms',         1e-06, ...
                 'SnrFixed',       3, ...
                 'ComputeKernel',  1, ...
                 'DataTypes',      {{'MEG GRAD', 'MEG MAG'}}));
        catch
            No_inverse_solution{count_no_inverse,1} = [participant{p} ' ' condition_mismatch_names{c} ' has no ' modality_data{mode} ' ' source_space{se} ' inverse solutions (computing failed)'];
            count_no_inverse = count_no_inverse +1;
            continue; % Go to next condition/modality
        end

        % Process: Change name in folder (for code)
        sFiles_name = bst_process('CallProcess', 'process_add_tag', sFiles_source, [], ...
            'tag',           ['Surf_' source_noise_tag '_' modality_data{mode}], ...
            'output',        2);  % Add to file path

        % Process: Set name to see in GUI
        sFiles_name = bst_process('CallProcess', 'process_set_comment', sFiles_name, [], ...
            'tag',           ['Surf_' source_noise_tag '_' modality_data{mode}], ...
            'isindex',       0);
        clear sFiles_source sFiles_name
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BIMODAL SURF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif strcmp(source_space{se},'Surf') && strcmp(modality_data{mode},'BIMODAL')
        
        % Check if the forward model for this modality and space exists before moving forward
        if ~exist([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c} '/headmodel_surf_duneuro_os_meg.mat'], 'file')
            No_inverse_solution{count_no_inverse,1} = [participant{p} ' ' condition_mismatch_names{c} ' has no ' modality_data{mode} ' ' source_space{se} ' inverse solutions (no forward model)'];
            count_no_inverse = count_no_inverse +1;
            continue; % Go to next condition/modality
        end    
                 
        disp(' ');
        disp('-------------------------');  
        disp(['Computing inverse solutions (' source_space{se} ' ' modality_data{mode} ') for ' participant{p}]);
        disp(datetime)
        disp(' ');   
        
        % If stated, find and delete any previous source kernel (overwrite it) before copying 
        if delete_previous_file == 1
            % check if there is already a source of THIS TYPE in the folder where we are gonna calculate it
            folders_delete = dir([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c}]);
            results_delete =  contains({folders_delete.name},'KERNEL') & endsWith({folders_delete.name}, ['Surf_' source_noise_tag '_' modality_data{mode} '.mat']); 
            infolder_delete = find(results_delete);
            if ~isempty(infolder_delete) % file exists, therefore delete it
               delete([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c} '/' folders_delete(infolder_delete).name]);
            end            
        end

        % NOW WE NEED TO SPECIFY WHICH FORWARD MODEL IS THE ONE SELECTED
        file_to_search = [participant{p} '/' condition_mismatch_names{c} '/channel_vectorview306_acc1.mat'];
        [sStudy, iStudy] = bst_get('AnyFile', file_to_search);
        % Possible ones here: surf_duneuro.mat, surf_os_meg.mat, surf_duneuro_os_meg.mat
        % Try to find the bimodal (only one useful here)
        surface_index = find(endsWith({sStudy.HeadModel.FileName},'surf_duneuro_os_meg.mat'), 1);
        % We now at this point that one of the two exists, but if more than one copy is present, show error
        if length(surface_index) > 1
            error('more than one surface head model found');
        end
        sStudy.iHeadModel = surface_index;
        bst_set('Study', iStudy, sStudy);
        % According to brainstorm developers, better not to reload
        % db_reload_studies(iStudy);

        % Finnally, compute the inverse solutions or keep track if not
        try
        % Process: Compute sources [2018]
        sFiles_source = bst_process('CallProcess', 'process_inverse_2018', sFiles, [], ...
            'output',  1, ...  % Kernel only: shared
            'inverse', struct(...
                 'Comment',        ['dSPM-unscaled: ' modality_data{mode} ' ALL'], ...
                 'InverseMethod',  'minnorm', ...
                 'InverseMeasure', 'dspm2018', ...
                 'SourceOrient',   {{'loose'}}, ...
                 'Loose',          0.4, ...
                 'UseDepth',       1, ...
                 'WeightExp',      0.5, ...
                 'WeightLimit',    10, ...
                 'NoiseMethod',    source_noise_option, ...
                 'NoiseReg',       0.1, ...
                 'SnrMethod',      'fixed', ...
                 'SnrRms',         1e-06, ...
                 'SnrFixed',       3, ...
                 'ComputeKernel',  1, ...
                 'DataTypes',      {{'MEG GRAD', 'MEG MAG', 'EEG'}}));
        catch
            No_inverse_solution{count_no_inverse,1} = [participant{p} ' ' condition_mismatch_names{c} ' has no ' modality_data{mode} ' ' source_space{se} ' inverse solutions (computing failed)'];
            count_no_inverse = count_no_inverse +1;
            continue; % Go to next condition/modality
        end

        % Process: Change name in folder (for code)
        sFiles_name = bst_process('CallProcess', 'process_add_tag', sFiles_source, [], ...
            'tag',           ['Surf_' source_noise_tag '_' modality_data{mode}], ...
            'output',        2);  % Add to file path

        % Process: Set name to see in GUI
        sFiles_name = bst_process('CallProcess', 'process_set_comment', sFiles_name, [], ...
            'tag',           ['Surf_' source_noise_tag '_' modality_data{mode}], ...
            'isindex',       0);
        clear sFiles_source sFiles_name
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EEG VOL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif strcmp(source_space{se},'Vol') && strcmp(modality_data{mode},'EEG')
            
        % Check if any forward model for this modality and space exists before moving forward
        if ~exist([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c} '/headmodel_vol_duneuro.mat'], 'file')...
            && ~exist([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c} '/headmodel_vol_duneuro_os_meg.mat'], 'file')
            % If BIMODAL one exists there should be no individual EEG one (merged with MEG) and that's still ok!
            No_inverse_solution{count_no_inverse,1} = [participant{p} ' ' condition_mismatch_names{c} ' has no ' modality_data{mode} ' ' source_space{se} ' inverse solutions (no forward model)'];
            count_no_inverse = count_no_inverse +1;
            continue; % Go to next condition/modality
        end    
                 
        disp(' ');
        disp('-------------------------');  
        disp(['Computing volume solutions (' source_space{se} ' ' modality_data{mode} ') for ' participant{p}]);
        disp(datetime)
        disp(' ');   
        
        % If stated, find and delete any previous source kernel (overwrite it) before copying 
        if delete_previous_file == 1
            % check if there is already a source of THIS TYPE in the folder where we are gonna calculate it
            folders_delete = dir([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c}]);
            results_delete =  contains({folders_delete.name},'KERNEL') & endsWith({folders_delete.name}, ['Vol_' source_noise_tag '_' modality_data{mode} '.mat']); 
            infolder_delete = find(results_delete);
            if ~isempty(infolder_delete) % file exists, therefore delete it
               delete([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c} '/' folders_delete(infolder_delete).name]);
            end            
        end

        % NOW WE NEED TO SPECIFY WHICH FORWARD MODEL IS THE ONE SELECTED
        file_to_search = [participant{p} '/' condition_mismatch_names{c} '/channel_vectorview306_acc1.mat'];
        [sStudy, iStudy] = bst_get('AnyFile', file_to_search);
        % Possible ones here: surf_duneuro.mat, surf_os_meg.mat, surf_duneuro_os_meg.mat
        % Try to find the bimodal first for this source space
        surface_index = find(endsWith({sStudy.HeadModel.FileName},'vol_duneuro_os_meg.mat'), 1);
        if isempty(surface_index)
            % Try the one specific for EEG then
            surface_index = find(endsWith({sStudy.HeadModel.FileName},'vol_duneuro.mat'), 1);
        end
        % We now at this point that one of the two exists, but if more than one copy is present, show error
        if length(surface_index) > 1
            error('more than one surface head model found');
        end
        sStudy.iHeadModel = surface_index;
        bst_set('Study', iStudy, sStudy);
        % According to brainstorm developers, better not to reload
        % db_reload_studies(iStudy);

        % Finnally, compute the inverse solutions or keep track if not
        try
        % Process: Compute sources [2018]
        sFiles_source = bst_process('CallProcess', 'process_inverse_2018', sFiles, [], ...
            'output',  1, ...  % Kernel only: shared
            'inverse', struct(...
                 'Comment',        ['dSPM-unscaled: ' modality_data{mode} ' ALL'], ...
                 'InverseMethod',  'minnorm', ...
                 'InverseMeasure', 'dspm2018', ...
                 'SourceOrient',   {{'free'}}, ...
                 'Loose',          0.4, ...
                 'UseDepth',       1, ...
                 'WeightExp',      0.5, ...
                 'WeightLimit',    10, ...
                 'NoiseMethod',    source_noise_option, ...
                 'NoiseReg',       0.1, ...
                 'SnrMethod',      'fixed', ...
                 'SnrRms',         1e-06, ...
                 'SnrFixed',       3, ...
                 'ComputeKernel',  1, ...
                 'DataTypes',      {{'EEG'}}));
        catch
            No_inverse_solution{count_no_inverse,1} = [participant{p} ' ' condition_mismatch_names{c} ' has no ' modality_data{mode} ' ' source_space{se} ' inverse solutions (computing failed)'];
            count_no_inverse = count_no_inverse +1;
            continue; % Go to next condition/modality
        end

        % Process: Change name in folder (for code)
        sFiles_name = bst_process('CallProcess', 'process_add_tag', sFiles_source, [], ...
            'tag',           ['Vol_' source_noise_tag '_' modality_data{mode}], ...
            'output',        2);  % Add to file path

        % Process: Set name to see in GUI
        sFiles_name = bst_process('CallProcess', 'process_set_comment', sFiles_name, [], ...
            'tag',           ['Vol_' source_noise_tag '_' modality_data{mode}], ...
            'isindex',       0);
        clear sFiles_source sFiles_name
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MEG VOL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif strcmp(source_space{se},'Vol') && strcmp(modality_data{mode},'MEG')
            
        % Check if any forward model for this modality and space exists before moving forward
        if ~exist([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c} '/headmodel_vol_os_meg.mat'], 'file')...
            && ~exist([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c} '/headmodel_vol_duneuro_os_meg.mat'], 'file')
            % If BIMODAL one exists there should be no individual EEG one (merged with MEG) and that's still ok!
            No_inverse_solution{count_no_inverse,1} = [participant{p} ' ' condition_mismatch_names{c} ' has no ' modality_data{mode} ' ' source_space{se} ' inverse solutions (no forward model)'];
            count_no_inverse = count_no_inverse +1;
            continue; % Go to next condition/modality
        end    
                 
        disp(' ');
        disp('-------------------------');  
        disp(['Computing volume solutions (' source_space{se} ' ' modality_data{mode} ') for ' participant{p}]);
        disp(datetime)
        disp(' ');   
        
        % If stated, find and delete any previous source kernel (overwrite it) before copying 
        if delete_previous_file == 1
            % check if there is already a source of THIS TYPE in the folder where we are gonna calculate it
            folders_delete = dir([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c}]);
            results_delete =  contains({folders_delete.name},'KERNEL') & endsWith({folders_delete.name}, ['Vol_' source_noise_tag '_' modality_data{mode} '.mat']); 
            infolder_delete = find(results_delete);
            if ~isempty(infolder_delete) % file exists, therefore delete it
               delete([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c} '/' folders_delete(infolder_delete).name]);
            end            
        end

        % NOW WE NEED TO SPECIFY WHICH FORWARD MODEL IS THE ONE SELECTED
        file_to_search = [participant{p} '/' condition_mismatch_names{c} '/channel_vectorview306_acc1.mat'];
        [sStudy, iStudy] = bst_get('AnyFile', file_to_search);
        % Possible ones here: surf_duneuro.mat, surf_os_meg.mat, surf_duneuro_os_meg.mat
        % Try to find the bimodal first for this source space
        surface_index = find(endsWith({sStudy.HeadModel.FileName},'vol_duneuro_os_meg.mat'), 1);
        if isempty(surface_index)
            % Try the one specific for MEG then
            surface_index = find(endsWith({sStudy.HeadModel.FileName},'vol_os_meg.mat'), 1);
        end
        % We know at this point that one of the two exists, but if more than one copy is present, show error
        if length(surface_index) > 1
            error('more than one surface head model found');
        end
        sStudy.iHeadModel = surface_index;
        bst_set('Study', iStudy, sStudy);
        % According to brainstorm developers, better not to reload
        % db_reload_studies(iStudy);

        % Finnally, compute the inverse solutions or keep track if not
        try
        % Process: Compute sources [2018]
        sFiles_source = bst_process('CallProcess', 'process_inverse_2018', sFiles, [], ...
            'output',  1, ...  % Kernel only: shared
            'inverse', struct(...
                 'Comment',        ['dSPM-unscaled: ' modality_data{mode} ' ALL'], ...
                 'InverseMethod',  'minnorm', ...
                 'InverseMeasure', 'dspm2018', ...
                 'SourceOrient',   {{'free'}}, ...
                 'Loose',          0.4, ...
                 'UseDepth',       1, ...
                 'WeightExp',      0.5, ...
                 'WeightLimit',    10, ...
                 'NoiseMethod',    source_noise_option, ...
                 'NoiseReg',       0.1, ...
                 'SnrMethod',      'fixed', ...
                 'SnrRms',         1e-06, ...
                 'SnrFixed',       3, ...
                 'ComputeKernel',  1, ...
                 'DataTypes',      {{'MEG GRAD', 'MEG MAG'}}));
        catch
            No_inverse_solution{count_no_inverse,1} = [participant{p} ' ' condition_mismatch_names{c} ' has no ' modality_data{mode} ' ' source_space{se} ' inverse solutions (computing failed)'];
            count_no_inverse = count_no_inverse +1;
            continue; % Go to next condition/modality
        end

        % Process: Change name in folder (for code)
        sFiles_name = bst_process('CallProcess', 'process_add_tag', sFiles_source, [], ...
            'tag',           ['Vol_' source_noise_tag '_' modality_data{mode}], ...
            'output',        2);  % Add to file path

        % Process: Set name to see in GUI
        sFiles_name = bst_process('CallProcess', 'process_set_comment', sFiles_name, [], ...
            'tag',           ['Vol_' source_noise_tag '_' modality_data{mode}], ...
            'isindex',       0);
        clear sFiles_source sFiles_name
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BIMODAL VOL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif strcmp(source_space{se},'Vol') && strcmp(modality_data{mode},'BIMODAL')
        
        % Check if any forward model for this modality and space exists before moving forward
        if ~exist([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c} '/headmodel_vol_duneuro_os_meg.mat'], 'file')
            No_inverse_solution{count_no_inverse,1} = [participant{p} ' ' condition_mismatch_names{c} ' has no ' modality_data{mode} ' ' source_space{se} ' inverse solutions (no forward model)'];
            count_no_inverse = count_no_inverse +1;
            continue; % Go to next condition/modality
        end    
                 
        disp(' ');
        disp('-------------------------');  
        disp(['Computing volume solutions (' source_space{se} ' ' modality_data{mode} ') for ' participant{p}]);
        disp(datetime)
        disp(' ');   
        
        % If stated, find and delete any previous source kernel (overwrite it) before copying 
        if delete_previous_file == 1
            % check if there is already a source of THIS TYPE in the folder where we are gonna calculate it
            folders_delete = dir([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c}]);
            results_delete =  contains({folders_delete.name},'KERNEL') & endsWith({folders_delete.name}, ['Vol_' source_noise_tag '_' modality_data{mode} '.mat']); 
            infolder_delete = find(results_delete);
            if ~isempty(infolder_delete) % file exists, therefore delete it
               delete([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c} '/' folders_delete(infolder_delete).name]);
            end            
        end

        % NOW WE NEED TO SPECIFY WHICH FORWARD MODEL IS THE ONE SELECTED
        file_to_search = [participant{p} '/' condition_mismatch_names{c} '/channel_vectorview306_acc1.mat'];
        [sStudy, iStudy] = bst_get('AnyFile', file_to_search);
        % Possible ones here: surf_duneuro.mat, surf_os_meg.mat, surf_duneuro_os_meg.mat
        % Try to find the bimodal first for this source space
        surface_index = find(endsWith({sStudy.HeadModel.FileName},'vol_duneuro_os_meg.mat'), 1);
        % We now at this point that one of the two exists, but if more than one copy is present, show error
        if length(surface_index) > 1
            error('more than one surface head model found');
        end
        sStudy.iHeadModel = surface_index;
        bst_set('Study', iStudy, sStudy);
        % According to brainstorm developers, better not to reload
        % db_reload_studies(iStudy);

        % Finnally, compute the inverse solutions or keep track if not
        try
        % Process: Compute sources [2018]
        sFiles_source = bst_process('CallProcess', 'process_inverse_2018', sFiles, [], ...
            'output',  1, ...  % Kernel only: shared
            'inverse', struct(...
                 'Comment',        ['dSPM-unscaled: ' modality_data{mode} ' ALL'], ...
                 'InverseMethod',  'minnorm', ...
                 'InverseMeasure', 'dspm2018', ...
                 'SourceOrient',   {{'free'}}, ...
                 'Loose',          0.4, ...
                 'UseDepth',       1, ...
                 'WeightExp',      0.5, ...
                 'WeightLimit',    10, ...
                 'NoiseMethod',    source_noise_option, ...
                 'NoiseReg',       0.1, ...
                 'SnrMethod',      'fixed', ...
                 'SnrRms',         1e-06, ...
                 'SnrFixed',       3, ...
                 'ComputeKernel',  1, ...
                 'DataTypes',      {{'MEG GRAD', 'MEG MAG', 'EEG'}}));
        catch
            No_inverse_solution{count_no_inverse,1} = [participant{p} ' ' condition_mismatch_names{c} ' has no ' modality_data{mode} ' ' source_space{se} ' inverse solutions (computing failed)'];
            count_no_inverse = count_no_inverse +1;
            continue; % Go to next condition/modality
        end

        % Process: Change name in folder (for code)
        sFiles_name = bst_process('CallProcess', 'process_add_tag', sFiles_source, [], ...
            'tag',           ['Vol_' source_noise_tag '_' modality_data{mode}], ...
            'output',        2);  % Add to file path

        % Process: Set name to see in GUI
        sFiles_name = bst_process('CallProcess', 'process_set_comment', sFiles_name, [], ...
            'tag',           ['Vol_' source_noise_tag '_' modality_data{mode}], ...
            'isindex',       0);
        clear sFiles_source sFiles_name
        end
    end
    end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % If successful, update subject_array for this subject
    subject_array{pos_subj,3} = 'sources_finished';
    save([root_dir '/subject_array.mat'],'subject_array')  
    end
end

% Reset modality_data variable to its original self
modality_data = {'EEG','MEG','BIMODAL'};

save([root_dir '/QC/No_inverse_solution_created.mat'],'No_inverse_solution');
clearvars('-except', initialVars{:});
disp 'DONE WITH COMPUTING NOISE COVARIANCE AND SOURCES (Project MMN)!!!'
disp(datetime)
toc

%% Project to common anatomy 

tic
disp(' ');      
disp('-------------------------');  
disp('PROJECT ABSOLUTE TO COMMON ANATOMY (Project MMN)'); 
disp(datetime)
disp('-------------------------');     
disp(' '); 

if source_analysis == 1
    modality_data = {'EEG'};
elseif source_analysis == 2
    modality_data = {'MEG'};
elseif source_analysis == 3
    modality_data = {'BIMODAL'};
elseif source_analysis == 4
    modality_data = {'EEG','MEG','BIMODAL'};
end

% Project to common anatomy
for p = 1:length(participant)  
    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    % Only if sources are computed we continue 
    if strcmp(subject_array{pos_subj,3},'sources_finished')
        
    % If sources-or-anatomy-specific problems, go to next subject
    if ~strcmp(subject_array{pos_subj,10},'corregistration_done') || ~strcmp(subject_array{pos_subj,11},'FEM_anatomy_ready')...
      || strcmp(subject_array{pos_subj,12},'No_HCP') || ~strcmp(subject_array{pos_subj,13},'EEG_projected')
        continue; % on to next subject
    end
    
    % Reload subject epoch folders first
    disp(' ');      
    disp('-------------------------');
    disp(['loading epoch folders for participant ' participant{p}]);
    disp(datetime)
    disp(' ');
    
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_conditions(current_sub); 
    
    % Reload anatomy folder of that subject too
    disp(' ');      
    disp('-------------------------');
    disp(['loading anatomy folder for participant ' participant{p}]);
    disp(datetime)
    disp(' ');
    
    prot_subs = bst_get('ProtocolSubjects');
    current_sub = find(strcmp({prot_subs.Subject.Name}, participant{p}));
    db_reload_subjects(current_sub);
    
    % Now for every condition, modality and source space, project to common anatomy
    for c = 1:length(condition_mismatch_names)
    for mode = 1:length(modality_data)
        for se = 1:length(source_space)
        % Find average and kernel files
        files = dir([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c}]);
        if isempty(files)
            error(['No ' condition_mismatch_names{c} ' files for ' participant{p}]);
        end
        % Select average accordingly
        if strcmp(condition_mismatch_names{c},'DeviantPitch-Standard') || strcmp(condition_mismatch_names{c},'DeviantDuration-Standard')
            infolder_average = find(contains({files.name},modality_data{mode}) & endsWith({files.name},'_MMN.mat'));
        else % normal conditions
            infolder_average = find(endsWith({files.name},[modality_data{mode} '_average.mat'])); 
        end
        % Select kernel
        infolder_kernel = find(contains({files.name},'KERNEL')...
            & endsWith({files.name},[source_space{se} '_' source_noise_tag '_' modality_data{mode} '.mat']));
        
        % Now it can happen that either the average or the kernel are not present
        % If the folder exists but not the Kernel or the average
        if isempty(infolder_kernel) || isempty(infolder_average); continue; end
        % If there are more than one kernels or averages (shouln't be) stop
        if length(infolder_average) > 1; error('more than one average of the same type'); end
        if length(infolder_kernel) > 1; error('more than one kernel of the same type'); end
        
        % Name the links, not the kernel or the avg file separately
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sFiles = {['link|' participant{p} '/' condition_mismatch_names{c} '/' ...
            files(infolder_kernel).name '|' ...
            participant{p} '/' condition_mismatch_names{c} '/' files(infolder_average).name]};
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        
        % If stated, find and delete any previous individual
        % sources in group average with projected sources
        if delete_previous_file == 1
            % check if there is already an individual source in
            % group folder (which would be the one with projected
            % sources into common anatomy)
            folders_delete = dir([root_dir_bs '/data/Group_analysis/' condition_mismatch_names{c} '/']);
            infolder_delete = find(endsWith({folders_delete.name}, [participant{p} '_' modality_data{mode} '_' source_space{se} '_' source_noise_tag '_' condition_mismatch_names{c} '.mat'])); 
            if ~isempty(infolder_delete) % file exists, therefore delete it
               delete([root_dir_bs '/data/Group_analysis/' condition_mismatch_names{c} '/' folders_delete(infolder_delete).name]);
            end
        end

        % Project to surface (it will automatically set them as norm)
        if strcmp(source_space{se},'Vol')
            % Ensure appropiate cortex/volume is loaded for VOLUME sources
            load([anat_path '/anat/@default_subject/brainstormsubject.mat']);
            if ~strcmp(Cortex,['@default_subject/' group_default_volume])
                Cortex = ['@default_subject/' group_default_volume];
                variableInfo = who('-file',[anat_path '/anat/@default_subject/brainstormsubject.mat']);
                save([anat_path '/anat/@default_subject/brainstormsubject.mat'],variableInfo{:});
            end
            % Then, reload group anatomy, as that folder is the destiny of the abs
            db_reload_subjects(0); % As 0 is the default anatomy
            % Also, reload Group analysis folder (important for volume grid)
            prot_subs = bst_get('ProtocolSubjects');
            current_sub = find(strcmp({prot_subs.Subject.Name}, 'Group_analysis'));
            db_reload_conditions(current_sub);
            % Will project on grid (same than bst_project_grid)
            sFiles = bst_process('CallProcess', 'process_project_sources', sFiles, [], ...
                'headmodeltype', 'volume');  % Volume
        elseif strcmp(source_space{se},'Surf')
            % Ensure appropiate cortex is loaded for cortical sources
            load([anat_path '/anat/@default_subject/brainstormsubject.mat']);
            if ~strcmp(Cortex,['@default_subject/' group_default_cortex])
                Cortex = ['@default_subject/' group_default_cortex];
                variableInfo = who('-file',[anat_path '/anat/@default_subject/brainstormsubject.mat']);
                save([anat_path '/anat/@default_subject/brainstormsubject.mat'],variableInfo{:});
            end
            % Then, reload group anatomy, as that folder is the destiny of the abs
            db_reload_subjects(0); % As 0 is the default anatomy
            % Process: Project on default anatomy: surface
            sFiles = bst_process('CallProcess', 'process_project_sources', sFiles, [], ...
                'headmodeltype', 'surface');  % Cortex surface
        end

        % Process: Add tag
        sFiles = bst_process('CallProcess', 'process_add_tag', sFiles, [], ...
            'tag',           [participant{p} '_' modality_data{mode} '_' source_space{se} '_' source_noise_tag '_' condition_mismatch_names{c}], ...
            'output',        2);  % Add to file name (1 to add a tag)

        % Process: Set name
        sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
            'tag',           [participant{p} '_' modality_data{mode} '_' source_space{se} '_' source_noise_tag '_' condition_mismatch_names{c}], ...
            'isindex',       1);
        
        end
    end
    end
    end
end

% Reset modality_data variable to its original self
modality_data = {'EEG','MEG','BIMODAL'};

clearvars('-except', initialVars{:});
disp 'DONE WITH PROJECT ABSOLUTE TO COMMON ANATOMY (Project MMN)!!!'
disp(datetime)
toc

%% GAVR inverse solutions

tic
disp(' ');      
disp('-------------------------');  
disp('GAVR INVERSE SOLUTIONS (Project MMN)'); 
disp(datetime)
disp('-------------------------');     
disp(' '); 

if source_analysis == 1
    modality_data = {'EEG'};
elseif source_analysis == 2
    modality_data = {'MEG'};
elseif source_analysis == 3
    modality_data = {'BIMODAL'};
elseif source_analysis == 4
    modality_data = {'EEG','MEG','BIMODAL'};
end

% Before GAVR, reload Group analysis folder again(since projected subject sources are there)
prot_subs = bst_get('ProtocolSubjects');
current_sub = find(strcmp({prot_subs.Subject.Name}, 'Group_analysis'));
db_reload_conditions(current_sub);

% GAVR
for c = 1:length(condition_mismatch_names)
    for mode = 1:length(modality_data)
    for se = 1:length(source_space)
        
        % Ensure the right cortex/volume is present in default anatomy
        % (Just in case)
        if strcmp(source_space{se},'Vol')
            load([anat_path '/anat/@default_subject/brainstormsubject.mat']);
            if ~strcmp(Cortex,['@default_subject/' group_default_volume])
                Cortex = ['@default_subject/' group_default_volume];
                variableInfo = who('-file',[anat_path '/anat/@default_subject/brainstormsubject.mat']);
                save([anat_path '/anat/@default_subject/brainstormsubject.mat'],variableInfo{:});
            end
            % Then, reload group anatomy, as that folder is the destiny of the abs
            db_reload_subjects(0); % As 0 is the default anatomy
            % Also, reload Group analysis folder (important for volume grid)
            prot_subs = bst_get('ProtocolSubjects');
            current_sub = find(strcmp({prot_subs.Subject.Name}, 'Group_analysis'));
            db_reload_conditions(current_sub);
        elseif strcmp(source_space{se},'Surf')
            load([anat_path '/anat/@default_subject/brainstormsubject.mat']);
            if ~strcmp(Cortex,['@default_subject/' group_default_cortex])
                Cortex = ['@default_subject/' group_default_cortex];
                variableInfo = who('-file',[anat_path '/anat/@default_subject/brainstormsubject.mat']);
                save([anat_path '/anat/@default_subject/brainstormsubject.mat'],variableInfo{:});
            end
            % Then, reload group anatomy, as that folder is the destiny of the abs
            db_reload_subjects(0); % As 0 is the default anatomy
        end
        
        for pg = 1:length(participant_group)
        sFiles = {};
        for p = 1:length(participant)
        % Identify patients/controls
        pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
        % If sources are not finished, go to next subject
        if ~strcmp(subject_array{pos_subj,3},'sources_finished')
            continue;
        end
        if ~strcmp(subject_array{pos_subj,2},participant_group{pg})
            continue; % so only include participants that correspond to the group
        end
            
        % Check if we should include it
        % General bad EEG/MEG
        subject_row = {subject_array{pos_subj,:}}; 
        subject_row = subject_row(~cellfun('isempty', subject_row')); % to avoid empty cells in contains
        bad_signal = find(contains(subject_row,['bad_' modality_data{mode}]));
        if strcmp(modality_data{mode},'BIMODAL') % In case of bimodal, check if any EEG or MEG is bad
            bad_signal_EEG = find(contains(subject_row,'bad_EEG'));
            bad_signal_MEG = find(contains(subject_row,'bad_MEG'));
            % If either one of these is bad, mark bimodal as bad
            if ~isempty(bad_signal_EEG) || ~isempty(bad_signal_MEG)
                bad_signal = 1;
            end
        end
        if ~isempty(bad_signal)
            warning([participant{p} ' not included in ' modality_data{mode} ' source average']);
            continue; % to next subject
        end     
        % Insuf covariance problems
        no_cov_sweeps = find(contains(subject_row,['insuf_cov_' modality_data{mode}]));
        if ~isempty(no_cov_sweeps)
            continue; % to next modality
        end
        % Sources-or-anatomy-specific problems
        if ~strcmp(subject_array{pos_subj,10},'corregistration_done') || ~strcmp(subject_array{pos_subj,11},'FEM_anatomy_ready')...
          || strcmp(subject_array{pos_subj,12},'No_HCP') || ~strcmp(subject_array{pos_subj,13},'EEG_projected')
            continue; % on to next subject
        end
        
        files = dir([root_dir_bs '/data/Group_analysis/' condition_mismatch_names{c}]);
        if isempty(files)
            error(['No ' condition_mismatch_names{c} ' Group analysis folder for ']);
        end
        infolder = find(endsWith({files.name},[participant{p} '_' modality_data{mode} '_' source_space{se} '_' source_noise_tag '_' condition_mismatch_names{c} '.mat']));
        if isempty(infolder)
           % In case, for instance, there is no EEG data for this subject
           warning(['No ' condition_mismatch_names{c} ' ' modality_data{mode} ' ' source_space{se} ' average for ' participant{p}]);
           continue;
        end  
        if length(infolder) > 1
            error(['More than one ' condition_mismatch_names{c} ' ' modality_data{mode} ' ' source_space{se} ' average for ' participant{p}]);
        end
        sFiles{p} = ['Group_analysis/' condition_mismatch_names{c} '/' files(infolder).name];
        end
        
        sFiles = sFiles(~cellfun('isempty', sFiles')); % to avoid empty cells
        if isempty(sFiles)
            % It can actually happen (BIMODAL Standard)
            warning(['No files to perform GAVR for ' condition_mismatch_names{c} ' ' modality_data{mode} ' ' source_space{se}]);
            continue;
        end
        
        gavr_n = num2str(length(sFiles));
        
        disp(' ');
        disp('-------------------------');  
        disp(['GAVR source data for ' condition_mismatch_names{c} ' ' modality_data{mode} ' ' source_space{se} ' ' participant_group{pg}]);
        disp(datetime)
        disp(' ');

        % If stated, find and delete any previous GAVR SENSOR data
        if delete_previous_file == 1
            % check if there is already GAVR source in Group analysis folder
            folders_delete = dir([root_dir_bs '/data/Group_analysis/' condition_mismatch_names{c}]);
            infolder_delete = find(contains({folders_delete.name},'GAVR_SOURCE_')...
            & endsWith({folders_delete.name}, [modality_data{mode} '_' condition_mismatch_names{c} '_' source_space{se} '_' participant_group{pg} '_n' gavr_n '_' gavr_name '.mat']));
            if ~isempty(infolder_delete) % file exists, therefore delete it
                for id = 1:length(infolder_delete)
                    inf_delete = infolder_delete(id);
                    delete([root_dir_bs '/data/Group_analysis/' condition_mismatch_names{c} '/' folders_delete(inf_delete).name]);
                end
            end
        end

        % Average using default function (because we are using vertices, not channels)
        sFiles = bst_process('CallProcess', 'process_average', sFiles, [], ...
            'avgtype',         1, ...  % Everything
            'avg_func',        1, ...  % Arithmetic average:  mean(x)
            'weighted',        0, ...  % We don't need weighted
            'scalenormalized', 0);

        % USE OPTION TO NORMALIZE HERE???

        % Process: Add tag
        sFiles = bst_process('CallProcess', 'process_add_tag', sFiles, [], ...
            'tag',           ['GAVR_SOURCE_' modality_data{mode} '_' condition_mismatch_names{c} '_' source_space{se} '_' participant_group{pg} '_n' gavr_n '_' gavr_name], ...
            'output',        2);  % Add to file name (1 to add a tag)

        % Process: Set name
        sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
            'tag',           ['GAVR_SOURCE_' modality_data{mode} '_' condition_mismatch_names{c} '_' source_space{se} '_' participant_group{pg} '_n' gavr_n '_' gavr_name], ...
            'isindex',       1);
        
        
        % Create a folder to save new GAVRs
        % [root_dir_bs '/data/Group_analysis/' condition_mismatch_names{c}]
        if ~exist([root_dir_bs '/data/Group_analysis/' gavr_name],'dir')
            iStudy = db_add_condition('Group_analysis', gavr_name, 1);
        end
        
        % Now copy the new gavr recently created into this folder (names of
        % variables are inherited from a previous script, no worries)
        folders_delete = dir([root_dir_bs '/data/Group_analysis/' condition_mismatch_names{c}]); % just to have the same variable name if it was created
        results_delete =  find(endsWith({folders_delete.name}, ['GAVR_SOURCE_' modality_data{mode} '_' condition_mismatch_names{c} '_' source_space{se} '_' participant_group{pg} '_n' gavr_n '_' gavr_name '.mat'])); 
        if ~isempty(results_delete)
            movefile([root_dir_bs '/data/Group_analysis/' condition_mismatch_names{c} '/' folders_delete(results_delete).name],[root_dir_bs '/data/Group_analysis/' gavr_name '/']);
        end
        end
    end
    end
end

% Reset modality_data variable to its original self
modality_data = {'EEG','MEG','BIMODAL'};

clearvars('-except', initialVars{:});
disp 'DONE WITH GAVR INVERSE SOLUTIONS (Project MMN)!!!'
disp(datetime)
toc

%% Extract source scouts (cortical MEG only for now)

tic
disp(' ');      
disp('-------------------------');  
disp('EXTRACTING SOURCE SCOUTS (Project MMN)');  
disp(datetime)
disp('-------------------------');     
disp(' '); 

% Change modality_data list depending on the choice of sources
if source_analysis == 1
    modality_data = {'EEG'};
elseif source_analysis == 2
    modality_data = {'MEG'};
elseif source_analysis == 3
    modality_data = {'BIMODAL'};
elseif source_analysis == 4
    modality_data = {'EEG','MEG','BIMODAL'};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% SPECIFIC PARAMETERS OF THIS SECTION %%%%%%%%%%%%%%%%%%%%%%
time_window_scout = [-0.05, 0.3];
scout_function = [1]; %#ok<*NBRAK> % [1,2,3,4] % 'mean','PCA','Max','Std'
scout_function_tag = {'mean'}; % 'mean','PCA','Max','Std'
normalization_option = 1; % 0 = NO; 1 = DO normalize across orientations % not sure how this works
if normalization_option == 0
    normalization_option_tag = '3D';
else
    normalization_option_tag = 'NORM';
end
atlas_to_look = 'Auditory_and_frontal'; % 'HCPMMP1' % 'Plots_left_right'
scouts = {'52_L','A1_L','A4_L','A5_L','LBelt_L', ...
    'MBelt_L','OP4_L','PBelt_L','RI_L','STSdp_L','STSvp_L', ...
    '52_R','A1_R','A4_R','A5_R','LBelt_R','MBelt_R', ...
    'OP4_R','PBelt_R','RI_R','STSdp_R','STSvp_R', ...
    '45_L','IFSa_L','IFSp_L','45_R','IFSa_R','IFSp_R' ...
    'OFC_L','pOFC_L','OFC_R','pOFC_R','IFG_L','IFG_R','OFC_L_merged','OFC_R_merged',...
    'AUDCORTEX_L','AUDCORTEX_R'};
scouts_peak = {'52_L','52_R','A1_L','A1_R','A4_L','A4_R','A5_L','A5_R','LBelt_L','LBelt_R','MBelt_L','MBelt_R',...
    'OP4_L','OP4_R','PBelt_L','PBelt_R','RI_L','RI_R','STSdp_L','STSdp_R','STSvp_L','STSvp_R',...
    '45_L','45_R','IFSa_L','IFSa_R','IFSp_L','IFSp_R','OFC_L','OFC_R','pOFC_L','pOFC_R',...
    'IFG_L','IFG_R','OFC_L_merged','OFC_R_merged','AUDCORTEX_L','AUDCORTEX_R'};
% Always in pairs of left and right
scouts_peak = {'A1_L','A1_R'}; % In case we want to extract the MMN peak average only for one ROI
scouts_label = scouts; % inherited from when we did it with HCPMMP1(see below)
extract_scouts_problems = {}; probl = 1;
% Define time samples first, which we are gonna need for mean peak amplitudes
plot_baseline = -50;
plot_post = 300;
s_r = 1000;
time_samples=linspace(plot_baseline,plot_post,((((plot_baseline*(-1)) + plot_post)/1000)*s_r) +1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reload group anatomy, just in case
db_reload_subjects(0); % As 0 is the default anatomy

% Also, reload Group analysis folder (where all norm single subjects are)
prot_subs = bst_get('ProtocolSubjects');
current_sub = find(strcmp({prot_subs.Subject.Name}, 'Group_analysis'));
db_reload_conditions(current_sub);

% Normalize sources before extracting (brainstorm)
for mode = 1:length(modality_data) 
for c = 1:length(condition_mismatch_names)
    sFiles = {};
    for p = 1:length(participant)  
        % Check if we should extract this one
        % General bad EEG/MEG
        pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
        subject_row = {subject_array{pos_subj,:}}; 
        subject_row = subject_row(~cellfun('isempty', subject_row')); % to avoid empty cells in contains
        bad_signal = find(contains(subject_row,['bad_' modality_data{mode}]));
        if strcmp(modality_data{mode},'BIMODAL') % In case of bimodal, check if any EEG or MEG is bad
            bad_signal_EEG = find(contains(subject_row,['bad_EEG']));
            bad_signal_MEG = find(contains(subject_row,['bad_MEG']));
            % If either one of these is bad, mark bimodal as bad
            if ~isempty(bad_signal_EEG) || ~isempty(bad_signal_MEG)
                bad_signal = 1;
            end
        end
        if ~isempty(bad_signal)
            warning([participant{p} ' source waveforms not extraced in ' modality_data{mode}]);
            continue; % to next subject
        end    
        % Insuf covariance problems
        no_cov_sweeps = find(contains(subject_row,['insuf_cov_' modality_data{mode}]));
        if ~isempty(no_cov_sweeps)
            continue; % to next modality
        end
        % Sources-or-anatomy-specific problems
        if ~strcmp(subject_array{pos_subj,10},'corregistration_done') || ~strcmp(subject_array{pos_subj,11},'FEM_anatomy_ready')...
          || strcmp(subject_array{pos_subj,12},'No_HCP') || ~strcmp(subject_array{pos_subj,13},'EEG_projected')
            continue; % on to next subject
        end

        % If sources are not finished, go to next subject
        if ~strcmp(subject_array{pos_subj,3},'sources_finished')
            continue;
        end
        
        % Reload subject epoch folders first
        disp(' ');      
        disp('-------------------------');
        disp(['loading epoch folders for participant ' participant{p}]);
        disp(datetime)
        disp(' ');
    
        [~, iStudies] = bst_get('StudyWithCondition', [participant{p} filesep condition_mismatch_names{c}]);
        db_reload_studies(iStudies);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Find average and kernel files
        files = dir([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c}]);
        if isempty(files)
            error(['No ' condition_mismatch_names{c} ' files for ' participant{p}]);
        end
        % Select average accordingly
        if strcmp(condition_mismatch_names{c},'DeviantPitch-Standard') || strcmp(condition_mismatch_names{c},'DeviantDuration-Standard')
            infolder_average = find(contains({files.name},modality_data{mode}) & endsWith({files.name},'_MMN.mat'));
        else % normal conditions
            infolder_average = find(endsWith({files.name},[modality_data{mode} '_average.mat'])); 
        end
        % Select kernel
        infolder_kernel = find(contains({files.name},'KERNEL')...
            & endsWith({files.name},['Surf_' source_noise_tag '_' modality_data{mode} '.mat']));
        % Now it can happen that either the average or the kernel are not present
        % If the folder exists but not the Kernel or the average
        if isempty(infolder_kernel) || isempty(infolder_average); continue; end
        % If there are more than one kernels or averages (shouln't be) stop
        if length(infolder_average) > 1; error('more than one average of the same type'); end
        if length(infolder_kernel) > 1; error('more than one kernel of the same type'); end

        % Name the links, not the kernel or the avg file separately
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        source_file = ['link|' participant{p} '/' condition_mismatch_names{c} '/' ...
            files(infolder_kernel).name '|' ...
            participant{p} '/' condition_mismatch_names{c} '/' files(infolder_average).name];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if isempty(source_file)
            warning(['subject ' participant{p} ' ' modality_data{mode} ' has no ' condition_mismatch_names{c} ' source']);
            continue;
        end
        if size(source_file) > 1
            warning(['subject ' participant{p} ' ' modality_data{mode} ' has more than one ' condition_mismatch_names{c} ' source']);
            continue;
        end

        sFiles = source_file;

        disp(' ');      
        disp('-------------------------');  
        disp(['Normalizing ' participant{p} ' ' condition_mismatch_names{c} ' Surf_' source_noise_tag '_' modality_data{mode} ' before extracting scouts']);
        disp(datetime)
        disp(' ');

        %%%%%%%%%%%%%%%%%%%% Delete previous %%%%%%%%%%%%%%%%%%%%%%%%%
        folders_delete = dir([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c} '/']);
        results_delete = contains({folders_delete.name}, 'results_norm') & endsWith({folders_delete.name},['Surf_' source_noise_tag '_' modality_data{mode} '_norm.mat']);
        infolder_delete = find(results_delete);
        if ~isempty(infolder_delete) % file exists, therefore delete it (any present there)
           for i = 1:length(infolder_delete)
               position = infolder_delete(i);
               delete([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c} '/' folders_delete(position).name]);
           end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Process: Unconstrained to flat map
        sFiles = bst_process('CallProcess', 'process_source_flat', sFiles, [], ...
            'method', 1);  % Norm: sqrt(x^2+y^2+z^2)

        % Process: Add tag: % It will already have results_norm by
        % default but make it end in that
        sFiles = bst_process('CallProcess', 'process_add_tag', sFiles, [], ...
            'tag',           ['Surf_' source_noise_tag '_' modality_data{mode} '_norm'], ... 
            'output',        2);  % Add to file name (1 to add a tag)

        % Process: Set name
        sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
            'tag',           ['Surf_' source_noise_tag '_' modality_data{mode} '_norm'], ...
            'isindex',       1);

    end
end
end

% Extract source waveforms for every condition (Brainstorm)
for so = 1:length(scout_function)
    for sc = 1:length(scouts)
        for mode = 1:length(modality_data) 
        for c = 1:length(condition_mismatch_names)
            sFiles = {};
            for p = 1:length(participant)  
                % Check if we should extract this one
                % General bad EEG/MEG
                pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
                subject_row = {subject_array{pos_subj,:}}; 
                subject_row = subject_row(~cellfun('isempty', subject_row')); % to avoid empty cells in contains
                bad_signal = find(contains(subject_row,['bad_' modality_data{mode}]));
                if strcmp(modality_data{mode},'BIMODAL') % In case of bimodal, check if any EEG or MEG is bad
                    bad_signal_EEG = find(contains(subject_row,['bad_EEG']));
                    bad_signal_MEG = find(contains(subject_row,['bad_MEG']));
                    % If either one of these is bad, mark bimodal as bad
                    if ~isempty(bad_signal_EEG) || ~isempty(bad_signal_MEG)
                        bad_signal = 1;
                    end
                end
                if ~isempty(bad_signal)
                    warning([participant{p} ' source waveforms not extraced in ' modality_data{mode}]);
                    continue; % to next subject
                end    
                % Insuf covariance problems
                no_cov_sweeps = find(contains(subject_row,['insuf_cov_' modality_data{mode}]));
                if ~isempty(no_cov_sweeps)
                    continue; % to next modality
                end
                % Sources-or-anatomy-specific problems
                if ~strcmp(subject_array{pos_subj,10},'corregistration_done') || ~strcmp(subject_array{pos_subj,11},'FEM_anatomy_ready')...
                  || strcmp(subject_array{pos_subj,12},'No_HCP') || ~strcmp(subject_array{pos_subj,13},'EEG_projected')
                    continue; % on to next subject
                end
                
                % If sources are not finished, go to next subject
                if ~strcmp(subject_array{pos_subj,3},'sources_finished')
                    continue;
                end
                
                % Reload subject epoch folders first
                disp(' ');      
                disp('-------------------------');
                disp(['loading epoch folders for participant ' participant{p}]);
                disp(datetime)
                disp(' ');

                [~, iStudies] = bst_get('StudyWithCondition', [participant{p} filesep condition_mismatch_names{c}]);
                db_reload_studies(iStudies);
                
                %%%%%%%%%%%%%%%%%%%% Find normalized file %%%%%%%%%%%%%%%%%%%%%%%%%
                folders = dir([root_dir_bs '/data/' participant{p} '/' condition_mismatch_names{c} '/']);
                results = contains({folders.name}, 'results_norm') & endsWith({folders.name},['Surf_' source_noise_tag '_' modality_data{mode} '_norm.mat']);
                infolder = find(results);
                if isempty(infolder) % file doesn't exist, so continue
                    continue;
                end
                if length(infolder) > 1
                    infolder = infolder(1); % In case there are more than one, just pick whichever, they should all be the same
                end
                
                norm_file = [participant{p} '/' condition_mismatch_names{c} '/' folders(infolder).name];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                sFiles{p} = norm_file; %#ok<*SAGROW>
            end
            % to avoid empty cells
            sFiles = sFiles(~cellfun('isempty', sFiles')); 

            disp(' ');      
            disp('-------------------------');  
            disp(['Extracting ' scout_function_tag{so} '_' modality_data{mode} '_' scouts{sc} '_' condition_mismatch_names{c}]);
            disp(datetime)
            disp(' ');  
            
            %%%%%%%%%%%%%%%%%%%% Delete previous %%%%%%%%%%%%%%%%%%%%%%%%%
            folders_delete = dir([root_dir_bs '/data/Group_analysis/' condition_mismatch_names{c} '/']);
            results_delete = contains({folders_delete.name}, 'matrix_scout') & endsWith({folders_delete.name},[scout_function_tag{so} '_' modality_data{mode} '_' scouts{sc} '_' condition_mismatch_names{c} '_' normalization_option_tag '.mat']);
            infolder_delete = find(results_delete);
            if ~isempty(infolder_delete) % file exists, therefore delete it
               delete([root_dir_bs '/data/Group_analysis/' condition_mismatch_names{c} '/' folders_delete(infolder_delete).name]);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Start a new report (as next step gives error that are
            % undetectable without this)
            bst_report('Start', sFiles);
            
            scout_array = scouts{sc};
            % Process: Scouts time series: 
            sFiles = bst_process('CallProcess', 'process_extract_scout', sFiles, [], ...
                'timewindow',     time_window_scout, ...
                'scouts',         {atlas_to_look, {scout_array}}, ...
                'scoutfunc',      scout_function(so), ...  % In case of trying different ones (mean, pca, etc)
                'isflip',         1, ...
                'isnorm',         normalization_option, ...
                'concatenate',    1, ...
                'save',           1, ...
                'addrowcomment',  1, ...
                'addfilecomment', 1);

            if isempty(sFiles) % means it was not executed
                % Save and display report
                ReportFile = bst_report('Save', sFiles);
                bst_report('Open', ReportFile);
                error(['problems for ' scout_function_tag{so} '_' scouts{sc} '_' condition_mismatch_names{c}])
            end
            
            % Process: Add tag
            sFiles = bst_process('CallProcess', 'process_add_tag', sFiles, [], ...
                'tag',           [scout_function_tag{so} '_' modality_data{mode} '_' scouts{sc} '_' condition_mismatch_names{c} '_' normalization_option_tag], ...
                'output',        2);  % Add to file name (1 to add a tag)

            % Process: Set name
            sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
                'tag',           ['Matrix_' scout_function_tag{so} '_' modality_data{mode} '_' scouts{sc} '_' condition_mismatch_names{c} '_' normalization_option_tag], ...
                'isindex',       1);
        end
        end
    end
end

% Extract source waveforms out of brainstorm for every condition
for so = 1:length(scout_function_tag)
    for mode = 1:length(modality_data)
    for c = 1:length(condition_mismatch_names)
        for p = 1:length(participant)
            % Check if we should extract this one
            % General bad EEG/MEG
            pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
            subject_row = {subject_array{pos_subj,:}}; 
            subject_row = subject_row(~cellfun('isempty', subject_row')); % to avoid empty cells in contains
            bad_signal = find(contains(subject_row,['bad_' modality_data{mode}]));
            if strcmp(modality_data{mode},'BIMODAL') % In case of bimodal, check if any EEG or MEG is bad
                bad_signal_EEG = find(contains(subject_row,['bad_EEG']));
                bad_signal_MEG = find(contains(subject_row,['bad_MEG']));
                % If either one of these is bad, mark bimodal as bad
                if ~isempty(bad_signal_EEG) || ~isempty(bad_signal_MEG)
                    bad_signal = 1;
                end
            end
            if ~isempty(bad_signal)
                warning([participant{p} ' source waveforms not extraced in ' modality_data{mode}]);
                continue; % to next subject
            end    
            % Insuf covariance problems
            no_cov_sweeps = find(contains(subject_row,['insuf_cov_' modality_data{mode}]));
            if ~isempty(no_cov_sweeps)
                continue; % to next modality
            end
            % Sources-or-anatomy-specific problems
            if ~strcmp(subject_array{pos_subj,10},'corregistration_done') || ~strcmp(subject_array{pos_subj,11},'FEM_anatomy_ready')...
              || strcmp(subject_array{pos_subj,12},'No_HCP') || ~strcmp(subject_array{pos_subj,13},'EEG_projected')
                continue; % on to next subject
            end
            
            % If sources are not finished, go to next subject
            if ~strcmp(subject_array{pos_subj,3},'sources_finished')
                continue;
            end
            
            if ~exist([root_dir '/Sources/Waveforms/' modality_data{mode} '/' participant{p}], 'dir')
                mkdir([root_dir '/Sources/Waveforms/' modality_data{mode} '/'], [participant{p} '_' normalization_option_tag]);
            end
            sources = [];
            for sc = 1:length(scouts_label)
                scout_file = dir([root_dir_bs '/data/Group_analysis/' condition_mismatch_names{c} '/matrix_scout_*_' scout_function_tag{so} '_' modality_data{mode} '_' scouts_label{sc} '_' condition_mismatch_names{c} '_' normalization_option_tag '.mat']);
                if isempty(scout_file)
                    warning(['No ' scout_function_tag{so} '_' modality_data{mode} '_' scouts_label{sc} '_' condition_mismatch_names{c} '_' normalization_option_tag ' scouts']);
                    continue;
                elseif size(scout_file) > 1
                    warning(['More than one ' scout_function_tag{so} '_' modality_data{mode} '_' scouts_label{sc} '_' condition_mismatch_names{c} '_' normalization_option_tag ' scouts']);
                    continue;
                end

                scout_file = [root_dir_bs '/data/Group_analysis/' condition_mismatch_names{c} '/' scout_file.name];
                eval(['load ' scout_file])
                subject_position = find(contains(Description,[participant{p} '/'])); % something unique about their name
                if isempty(subject_position)
                    warning(['subject ' participant{p} 'is not in the scout file'])
                    continue; % With next participant
                end
                
                disp(' ');      
                disp('-------------------------');  
                disp(['Extracting source waveforms out of bs for ' participant{p} ' ' scouts_label{sc} '_' modality_data{mode} '_' scout_function_tag{so} '_' condition_mismatch_names{c} '_' normalization_option_tag]);
                disp(datetime)
                disp(' '); 
                
                if normalization_option == 0
                    for i = 1:length(subject_position)% for each of the three orientations
                    F(sc,:,i) = Value(subject_position(i),:); %#ok<*SAGROW>
                    % So these files will have three dimensions, the third one being orientation
                    end
                else
                    F(sc,:) =  Value(subject_position,:); %#ok<*SAGROW>
                end 
            end
            % Once all scouts are stored in variable sources, save
            save([root_dir '/Sources/Waveforms/' modality_data{mode} '/' participant{p} '_' normalization_option_tag '/' scout_function_tag{so} '_' condition_mismatch_names{c} '.mat'],'F');   
        end
    end
    end
end

% Prepare GAVR folders outside brainstorm
for mode = 1:length(modality_data)
% root_dir '/Sources/Waveforms/' modality_data{mode} '/'], [participant{p} '_' normalization_option_tag]
if ~exist([root_dir '/Sources/Waveforms/' modality_data{mode} '/' gavr_name], 'dir')
    mkdir([root_dir '/Sources/Waveforms/' modality_data{mode} '/'], [gavr_name '_' normalization_option_tag]);
    mkdir([root_dir '/Sources/Waveforms/' modality_data{mode} '/'], [gavr_name '_' normalization_option_tag '/gavr']);
    mkdir([root_dir '/Sources/Waveforms/' modality_data{mode} '/'], [gavr_name '_' normalization_option_tag '/std_dev']);
    mkdir([root_dir '/Sources/Waveforms/' modality_data{mode} '/'], [gavr_name '_' normalization_option_tag '/std_err']);
end
end

% GAVR source waveforms for every condition outside brainstorm 
% (and extract Matrix)
for mode = 1:length(modality_data)
    % Create directory to save if not present
    if ~exist([root_dir '/Sources/Waveforms/' modality_data{mode} '/' gavr_name '_' normalization_option_tag],'dir')
        mkdir([root_dir '/Sources/Waveforms/' modality_data{mode} '/'], [gavr_name '_' normalization_option_tag]);
    end
for so = 1:length(scout_function_tag)
    for c = 1:length(condition_mismatch_names)
    for pg = 1:length(participant_group)    
        sources = [];
        for p = 1:length(participant)
        % Identify patients/controls
        pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
        if ~strcmp(subject_array{pos_subj,2},participant_group{pg})
            continue; % so only include participants that correspond to the group
        end
            
        % Check if we should extract this one
        % General bad EEG/MEG
        subject_row = {subject_array{pos_subj,:}}; 
        subject_row = subject_row(~cellfun('isempty', subject_row')); % to avoid empty cells in contains
        bad_signal = find(contains(subject_row,['bad_' modality_data{mode}]));
        if strcmp(modality_data{mode},'BIMODAL') % In case of bimodal, check if any EEG or MEG is bad
            bad_signal_EEG = find(contains(subject_row,['bad_EEG']));
            bad_signal_MEG = find(contains(subject_row,['bad_MEG']));
            % If either one of these is bad, mark bimodal as bad
            if ~isempty(bad_signal_EEG) || ~isempty(bad_signal_MEG)
                bad_signal = 1;
            end
        end
        if ~isempty(bad_signal)
            warning([participant{p} ' source waveforms not extraced in ' modality_data{mode}]);
            continue; % to next subject
        end     
        % Sources-or-anatomy-specific problems
        if ~strcmp(subject_array{pos_subj,10},'corregistration_done') || ~strcmp(subject_array{pos_subj,11},'FEM_anatomy_ready')...
          || strcmp(subject_array{pos_subj,12},'No_HCP') || ~strcmp(subject_array{pos_subj,13},'EEG_projected')
            warning([participant{p} ' source waveforms not extraced in ' modality_data{mode}]);      
            continue; % on to next subject
        end

        % If sources are not finished, go to next subject
        if ~strcmp(subject_array{pos_subj,3},'sources_finished')
            continue;
        end
        
        for sc = 1:length(scouts_label)
            scout_file = dir([root_dir_bs '/data/Group_analysis/' condition_mismatch_names{c} '/matrix_scout_*_' scout_function_tag{so} '_' modality_data{mode} '_' scouts_label{sc} '_' condition_mismatch_names{c} '_' normalization_option_tag '.mat']);
            if isempty(scout_file)
                warning(['No ' scout_function_tag{so} '_' modality_data{mode} '_' scouts_label{sc} '_' condition_mismatch_names{c} '_' normalization_option_tag ' scouts']);
                continue;
            elseif size(scout_file) > 1
                warning(['More than one ' scout_function_tag{so} '_' modality_data{mode} '_' scouts_label{sc} '_' condition_mismatch_names{c} '_' normalization_option_tag ' scouts']);
                continue;
            end

            scout_file = [root_dir_bs '/data/Group_analysis/' condition_mismatch_names{c} '/' scout_file.name];
            eval(['load ' scout_file])
            subject_position = find(contains(Description,[participant{p} '/'])); % something unique about their name
            if isempty(subject_position)
                warning(['subject ' participant{p} ' is not in the scout file'])
                continue; % on to next subject
            end
            if normalization_option == 0
                for i = 1:length(subject_position)% for each of the three orientations
                sources(sc,:,i,p) = Value(subject_position(i),:); %#ok<*SAGROW>
                % So these files will have four dimensions, the third one being orientation
                end
            else
                sources(sc,:,p) =  Value(subject_position,:); %#ok<*SAGROW>
            end 
        end
        end
        
        disp(' ');      
        disp('-------------------------');  
        disp(['GAVR source waveforms for ' scout_function_tag{so} '_' modality_data{mode} '_' condition_mismatch_names{c} '_' participant_group{pg} '_' normalization_option_tag]);
        disp(datetime)
        disp(' '); 
        
        % Redefine source matrix eliminating empty subjects
        non_empty_subj_idx = []; pos_nes = 1;
        if normalization_option == 1 % only one orientation, 3D variable
            [rows, columns, slices] = size(sources);
            for sli = 1:slices % subjects, third dimension
                thisslice = sources(:,:,sli);
                % Check if it's all zeros (means all ROIs and time points for that subject are empty)
                 if all(thisslice(:) == 0)
                    disp(['disregarding empty participant']);
                 else
                     % save indices of subjects that are non-zeros
                     non_empty_subj_idx(pos_nes) = sli;
                     pos_nes = pos_nes +1;
                 end
            end
            sources2 = sources(:,:,non_empty_subj_idx);
        elseif normalization_option == 0 % 4D variable
            [rows, columns, third_dim, slices] = size(sources);
            for sli = 1:slices % subjects, third dimension
                thisslice = sources(:,:,:,sli);
                % Check if it's all zeros (means all ROIs and time points for that subject are empty)
                 if all(thisslice(:) == 0)
                    disp(['disregarding empty participant']);
                 else
                     % save indices of subjects that are non-zeros
                     non_empty_subj_idx(pos_nes) = sli;
                     pos_nes = pos_nes +1;
                 end
            end
            sources2 = sources(:,:,:,non_empty_subj_idx);
        end
        
        % Non empty matrix
        sources = sources2;
        
        % Once all scouts and subjects are stored in variables, average and save
        if normalization_option == 0
            % So these files will have four dimensions, the third one being orientation
            gavr = mean(sources,4);
            std_dev = std(sources,0,4);
            std_err = std_dev/sqrt(size(sources,4));   
        else
            gavr = mean(sources,3);
            std_dev = std(sources,0,3);
            std_err = std_dev/sqrt(size(sources,3));
        end 
    
        
        
        % Save Matrix (easier to access data this way later)
        save([root_dir '/Sources/Waveforms/' modality_data{mode} '/' gavr_name '_' normalization_option_tag '/Matrix_' scout_function_tag{so} '_' condition_mismatch_names{c} '_' participant_group{pg} '.mat'],'sources');   
        
        F = gavr; %#ok<*NASGU>
        save([root_dir '/Sources/Waveforms/' modality_data{mode} '/' gavr_name '_' normalization_option_tag '/gavr/' scout_function_tag{so} '_' condition_mismatch_names{c} '_' participant_group{pg} '.mat'],'F');   
        F = std_dev;
        save([root_dir '/Sources/Waveforms/' modality_data{mode} '/' gavr_name '_' normalization_option_tag '/std_dev/' scout_function_tag{so} '_' condition_mismatch_names{c} '_' participant_group{pg} '.mat'],'F');   
        F = std_err;
        save([root_dir '/Sources/Waveforms/' modality_data{mode} '/' gavr_name '_' normalization_option_tag '/std_err/' scout_function_tag{so} '_' condition_mismatch_names{c} '_' participant_group{pg} '.mat'],'F');   
    end
    end
end
end

% Extract MMN for every ROI (MMN conditions only) and create SPSS matrix
for mode = 1:length(modality_data)
    % To use for SPSS matrix later
    temporal = []; iter = 1;
    % Create directory to save if not present
    if ~exist([root_dir '/Sources/Waveforms/' modality_data{mode} '/' gavr_name '_' normalization_option_tag],'dir')
        mkdir([root_dir '/Sources/Waveforms/' modality_data{mode} '/'], [gavr_name '_' normalization_option_tag]);
    end
for so = 1:length(scout_function_tag)
    header_Project_MMN_source = {};
    for pg = 1:length(participant_group) 
    for c = 4:length(condition_mismatch_names) % ONLY MISMATCH CONDITIONS FOR PEAK EXTRACTION!
        sources = [];
        for p = 1:length(participant)
        % Identify patients/controls
        pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
        if ~strcmp(subject_array{pos_subj,2},participant_group{pg})
            continue; % so only include participants that correspond to the group
        end
            
        % Check if we should extract this one
        % General bad EEG/MEG
        subject_row = {subject_array{pos_subj,:}}; 
        subject_row = subject_row(~cellfun('isempty', subject_row')); % to avoid empty cells in contains
        bad_signal = find(contains(subject_row,['bad_' modality_data{mode}]));
        if strcmp(modality_data{mode},'BIMODAL') % In case of bimodal, check if any EEG or MEG is bad
            bad_signal_EEG = find(contains(subject_row,['bad_EEG']));
            bad_signal_MEG = find(contains(subject_row,['bad_MEG']));
            % If either one of these is bad, mark bimodal as bad
            if ~isempty(bad_signal_EEG) || ~isempty(bad_signal_MEG)
                bad_signal = 1;
            end
        end
        if ~isempty(bad_signal)
            warning([participant{p} ' source waveforms not extraced in ' modality_data{mode}]);
            continue; % to next subject
        end     
        % Sources-or-anatomy-specific problems
        if ~strcmp(subject_array{pos_subj,10},'corregistration_done') || ~strcmp(subject_array{pos_subj,11},'FEM_anatomy_ready')...
          || strcmp(subject_array{pos_subj,12},'No_HCP') || ~strcmp(subject_array{pos_subj,13},'EEG_projected')
            warning([participant{p} ' source waveforms not extraced in ' modality_data{mode}]);      
            continue; % on to next subject
        end

        % If sources are not finished, go to next subject
        if ~strcmp(subject_array{pos_subj,3},'sources_finished')
            continue;
        end
        
        for sc = 1:length(scouts_label)
            scout_file = dir([root_dir_bs '/data/Group_analysis/' condition_mismatch_names{c} '/matrix_scout_*_' scout_function_tag{so} '_' modality_data{mode} '_' scouts_label{sc} '_' condition_mismatch_names{c} '_' normalization_option_tag '.mat']);
            if isempty(scout_file)
                warning(['No ' scout_function_tag{so} '_' modality_data{mode} '_' scouts_label{sc} '_' condition_mismatch_names{c} '_' normalization_option_tag ' scouts']);
                continue;
            elseif size(scout_file) > 1
                warning(['More than one ' scout_function_tag{so} '_' modality_data{mode} '_' scouts_label{sc} '_' condition_mismatch_names{c} '_' normalization_option_tag ' scouts']);
                continue;
            end

            scout_file = [root_dir_bs '/data/Group_analysis/' condition_mismatch_names{c} '/' scout_file.name];
            eval(['load ' scout_file])
            subject_position = find(contains(Description,[participant{p} '/'])); % something unique about their name
            if isempty(subject_position)
                warning(['subject ' participant{p} ' is not in the scout file'])
                continue; % on to next subject
            end
            if normalization_option == 0
                for i = 1:length(subject_position)% for each of the three orientations
                sources(sc,:,i,p) = Value(subject_position(i),:); %#ok<*SAGROW>
                % So these files will have four dimensions, the third one being orientation
                end
            else
                sources(sc,:,p) =  Value(subject_position,:); %#ok<*SAGROW>
            end 
        end
        end
        
        % Redefine source matrix eliminating empty subjects
        non_empty_subj_idx = []; pos_nes = 1;
        if normalization_option == 1 % only one orientation, 3D variable
            [rows, columns, slices] = size(sources);
            for sli = 1:slices % subjects, third dimension
                thisslice = sources(:,:,sli);
                % Check if it's all zeros (means all ROIs and time points for that subject are empty)
                 if all(thisslice(:) == 0)
                    disp(['disregarding empty participant']);
                 else
                     % save indices of subjects that are non-zeros
                     non_empty_subj_idx(pos_nes) = sli;
                     pos_nes = pos_nes +1;
                 end
            end
            sources2 = sources(:,:,non_empty_subj_idx);
        elseif normalization_option == 0 % 4D variable
            [rows, columns, third_dim, slices] = size(sources);
            for sli = 1:slices % subjects, third dimension
                thisslice = sources(:,:,:,sli);
                % Check if it's all zeros (means all ROIs and time points for that subject are empty)
                 if all(thisslice(:) == 0)
                    disp(['disregarding empty participant']);
                 else
                     % save indices of subjects that are non-zeros
                     non_empty_subj_idx(pos_nes) = sli;
                     pos_nes = pos_nes +1;
                 end
            end
            sources2 = sources(:,:,:,non_empty_subj_idx);
        end
        
        % Non empty matrix
        sources = sources2;

        % Ensure correct number of scouts
        if size(sources,1) ~= length(scouts)
            error('watch out, matrix size does not match with scout labels')
        end

        % Retrieve values from each of the scouts from which we want peaks
        for sco_p = 1:length(scouts_peak)
            
        disp(' ');      
        disp('-------------------------');  
        disp(['Extracting ' scouts_peak{sco_p} ' mean peak amplitudes for ' scout_function_tag{so} '_' modality_data{mode} '_' condition_mismatch_names{c} '_' participant_group{pg} '_' normalization_option_tag]);
        disp(datetime)
        disp(' '); 
            
        % Find position of that scout in scout list
        sco_num = find(strcmp(scouts,scouts_peak{sco_p}));
        if normalization_option == 0
            Scout_amplitudes = sources(sco_num,:,:,:);
        else
            
            Scout_amplitudes = squeeze(sources(sco_num,:,:));
               for pe = 1:length(Peaks_to_extract)
               % Find the indices of the start and end of time window from time_samples variable
               if strcmp(condition_mismatch_names{c},'DeviantPitch-Standard') && endsWith(scouts_peak{sco_p},'_L')
                   eval(['time_window = time_window_source_left_Pitch_' Peaks_to_extract{pe} ';'])
               elseif strcmp(condition_mismatch_names{c},'DeviantDuration-Standard') && endsWith(scouts_peak{sco_p},'_L')
                   eval(['time_window = time_window_source_left_Dur_' Peaks_to_extract{pe} ';'])
               elseif strcmp(condition_mismatch_names{c},'DeviantPitch-Standard') && endsWith(scouts_peak{sco_p},'_R')
                   eval(['time_window = time_window_source_right_Pitch_' Peaks_to_extract{pe} ';'])
               elseif strcmp(condition_mismatch_names{c},'DeviantDuration-Standard') && endsWith(scouts_peak{sco_p},'_R')
                   eval(['time_window = time_window_source_right_Dur_' Peaks_to_extract{pe} ';'])
               elseif strcmp(condition_mismatch_names{c},'DeviantPitch-Standard') && endsWith(scouts_peak{sco_p},'_L_merged')
                   eval(['time_window = time_window_source_left_Pitch_' Peaks_to_extract{pe} ';'])
               elseif strcmp(condition_mismatch_names{c},'DeviantDuration-Standard') && endsWith(scouts_peak{sco_p},'_L_merged')
                   eval(['time_window = time_window_source_left_Dur_' Peaks_to_extract{pe} ';'])
               elseif strcmp(condition_mismatch_names{c},'DeviantPitch-Standard') && endsWith(scouts_peak{sco_p},'_R_merged')
                   eval(['time_window = time_window_source_right_Pitch_' Peaks_to_extract{pe} ';'])
               elseif strcmp(condition_mismatch_names{c},'DeviantDuration-Standard') && endsWith(scouts_peak{sco_p},'_R_merged')
                   eval(['time_window = time_window_source_right_Dur_' Peaks_to_extract{pe} ';'])
               else
                   error('Peaks have to be extracted at either Pitch or Dur MMN')
               end
               [~,closestIndex] = min(abs(time_samples-time_window(1)));
               init_time = closestIndex;
               [~,closestIndex] = min(abs(time_samples-time_window(2)));
               end_time = closestIndex;
               
               F = mean(Scout_amplitudes(init_time:end_time,:),1);

               % Create directory to save if not present
               if ~exist([root_dir '/Sources/Waveforms/' modality_data{mode} '/' gavr_name '_' normalization_option_tag],'dir')
                    mkdir([root_dir '/Sources/Waveforms/' modality_data{mode} '/'], [gavr_name '_' normalization_option_tag]);
               end

               % Now save the Peak averages for all subjects for this condition and ROI
               save([root_dir '/Sources/Waveforms/' modality_data{mode} '/' gavr_name '_' normalization_option_tag '/' Peaks_to_extract{pe} '_' scout_function_tag{so} '_' condition_mismatch_names{c} '_' scouts_peak{sco_p} '_' participant_group{pg} '.mat'],'F');  
               
               % Before moving forward, store that F in the SPSS matrix
               temporal(1:length(F),iter) = F;
               % Add the corresponding title to the header for later
               header_Project_MMN_source{iter} =  [participant_group{pg} '_' condition_short_labels{c} '_' scouts_peak{sco_p}];
               iter = iter + 1;
               end
        end
        end
    end
    end
    
    % Now create a Matrix with the peaks for later use in SPSS 
    % (FE > Cond > ROI > hemisphere)
    SPSS_matrix = [header_Project_MMN_source; num2cell(temporal)];
    if ~exist([root_dir '/Sources/Waveforms/' modality_data{mode} '/' gavr_name '_' normalization_option_tag '/Statistics'],'dir')
        mkdir([root_dir '/Sources/Waveforms/' modality_data{mode} '/' gavr_name '_' normalization_option_tag '/'], 'Statistics');
    end
    save([root_dir '/Sources/Waveforms/' modality_data{mode} '/' gavr_name '_' normalization_option_tag '/Statistics/SPSS_matrix_' scout_function_tag{so} '.mat'],'SPSS_matrix');
    xlswrite([root_dir '/Sources/Waveforms/' modality_data{mode} '/' gavr_name '_' normalization_option_tag '/Statistics/SPSS_matrix_' scout_function_tag{so} '.xlsx'],SPSS_matrix);
end
end

% Reset modality_data variable to its original self
modality_data = {'EEG','MEG','BIMODAL'};
save([root_dir '/QC/extract_scouts_problems_.mat'],'extract_scouts_problems');

clearvars('-except', initialVars{:});
disp 'DONE WITH EXTRACTING SOURCE SCOUTS (Project MMN)!!!'
disp(datetime)
toc

%% Define Plot-specific variables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% SPECIFIC PARAMETERS OF THIS SECTION %%%%%%%%%%%%%%%%%%%%%%
set(0,'defaultfigurecolor',[1 1 1]); % I want white backgrounds
abreviated_cond_titles = {'STD','PDEV','DDEV','PitchMMN','DurMMN'};
hemisphere = {'L','R'};
% ROIs = {'MBelt','A1','LBelt','PBelt','STSdp','STSvp'}; % The ones you want to plot...
ROIs = {'A1','AUDCORTEX','pOFC','IFG'}; % The ones you want to plot...
%... based on scouts extracted (Has to be the same than 'Extract source scouts' section)
scout_positions = {'52_L','A1_L','A4_L','A5_L','LBelt_L', ...
    'MBelt_L','OP4_L','PBelt_L','RI_L','STSdp_L','STSvp_L', ...
    '52_R','A1_R','A4_R','A5_R','LBelt_R','MBelt_R', ...
    'OP4_R','PBelt_R','RI_R','STSdp_R','STSvp_R', ...
    '45_L','IFSa_L','IFSp_L','45_R','IFSa_R','IFSp_R' ...
    'OFC_L','pOFC_L','OFC_R','pOFC_R','IFG_L','IFG_R','OFC_L_merged','OFC_R_merged',...
    'AUDCORTEX_L','AUDCORTEX_R'};

normalization_option_tag = 'NORM'; % Or '3D'
scout_function_string = 'mean'; % 'mean','PCA','Max','Std'
dev_GAVR = 2; % 1 = Standard Deviation 2 = Standard Error
plot_baseline = -50;
plot_post = 300;
s_r = 1000;
time_samples=linspace(plot_baseline,plot_post,((((plot_baseline*(-1)) + plot_post)/1000)*s_r) +1);
% Colors
color_group = {[255 0 0],[0 0 0]}; % FE C
color_hem = {[0 0 255],[255 0 0]}; % L R
color_roi = {[255 0 255],[255 0 0],[150 0 0],...
    [255 128 0],[0 153 76],[0 0 153]}; % Depends on how many you put above (bot only one hem)
transparency = 0.2;
shaded_areas = 1; % 0 = NO shaded areas; 1 = YES
shaded_areas_scalp_Pitch_MMN_EEG = {[time_window_scalp_Pitch_MMN_EEG(1) time_window_scalp_Pitch_MMN_EEG(2) time_window_scalp_Pitch_MMN_EEG(2) time_window_scalp_Pitch_MMN_EEG(1)]}; 
shaded_areas_scalp_Pitch_MMN_MEG = {[time_window_scalp_Pitch_MMN_MEG(1) time_window_scalp_Pitch_MMN_MEG(2) time_window_scalp_Pitch_MMN_MEG(2) time_window_scalp_Pitch_MMN_MEG(1)]}; 
shaded_areas_scalp_Dur_MMN_EEG = {[time_window_scalp_Dur_MMN_EEG(1) time_window_scalp_Dur_MMN_EEG(2) time_window_scalp_Dur_MMN_EEG(2) time_window_scalp_Dur_MMN_EEG(1)]}; 
shaded_areas_scalp_Dur_MMN_MEG = {[time_window_scalp_Dur_MMN_MEG(1) time_window_scalp_Dur_MMN_MEG(2) time_window_scalp_Dur_MMN_MEG(2) time_window_scalp_Dur_MMN_MEG(1)]};
shaded_areas_source_left_Pitch_MMN = {[time_window_source_left_Pitch_MMN(1) time_window_source_left_Pitch_MMN(2) time_window_source_left_Pitch_MMN(2) time_window_source_left_Pitch_MMN(1)]}; % {[110 130 130 110]}; % [100 120 120 100] 
shaded_areas_source_left_Dur_MMN = {[time_window_source_left_Dur_MMN(1) time_window_source_left_Dur_MMN(2) time_window_source_left_Dur_MMN(2) time_window_source_left_Dur_MMN(1)]};
shaded_areas_source_right_Pitch_MMN = {[time_window_source_right_Pitch_MMN(1) time_window_source_right_Pitch_MMN(2) time_window_source_right_Pitch_MMN(2) time_window_source_right_Pitch_MMN(1)]}; % {[110 130 130 110]}; % [100 120 120 100]
shaded_areas_source_right_Dur_MMN = {[time_window_source_right_Dur_MMN(1) time_window_source_right_Dur_MMN(2) time_window_source_right_Dur_MMN(2) time_window_source_right_Dur_MMN(1)]};

save_figures_option = 0; % 0 = NO 1 = YES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot butterfly scalp (EEG/MEG(2Fig); C/FE; 5cond; 2x5)

% Change modality_data list depending on the choice of sources
modality_data = {'EEG','MEG'};

% Predefine which channels/clusters are going to be used from matrix
load([root_dir '/Scripts/Areas_channels.mat']);
% EEG channels
if size(choice_channel_EEG,2) > 1 % cluster was selected
    channel_pos_EEG = [];
    for i = 1:size(choice_channel_EEG,2) % Find all the indices
        channel_pos_EEG(i) = find(strcmp({Areas.Name}, choice_channel_EEG{i})==1);
    end
else % single channel
   channel_pos_EEG = find(strcmp({Areas.Name}, choice_channel_EEG)==1); 
end
channel_pos_EEG = channel_pos_EEG - 314; % to adjust for 62 channels
% MEG sensors
if size(choice_channel_MEG,2) > 1 % cluster was selected
    channel_pos_MEG = [];
    for i = 1:size(choice_channel_MEG,2) % Find all the indices
        channel_pos_MEG(i) = find(strcmp({Areas.Name}, choice_channel_MEG{i})==1);
    end
else % single sensor
   channel_pos_MEG = find(strcmp({Areas.Name}, choice_channel_MEG)==1); 
end

for mode = 1:length(modality_data)
    % New plot 
    figure('units','normalized','outerposition',[0 0 1 1]);
    num_subplots = 10; % For later
    h(1) = subplot(2,5,1); h(2) = subplot(2,5,2); h(3) = subplot(2,5,3); 
    h(4) = subplot(2,5,4);h(5) = subplot(2,5,5); h(6) = subplot(2,5,6);
    h(7) = subplot(2,5,7);h(8) = subplot(2,5,8); h(9) = subplot(2,5,9);
    h(10) = subplot(2,5,10);
    current_max = []; current_min = []; pos_ylim = 1; % To adjust ylim later
    % Define channel number
    eval(['chan_num = channel_pos_' modality_data{mode} ';'])
for c = 1:length(condition_mismatch_names)
for pg = 1:length(participant_group)
    % Load data matrix
    load([root_dir '/Scalp/' gavr_name '/Matrix_' modality_data{mode} '_' condition_mismatch_names{c} '_' participant_group{pg} '.mat']);
    % Ensure number of channels is correct
    if size(Matrix,1) ~= 306 && strcmp(modality_data{mode},'MEG')
        error(['MEG DATA MATRIX HAS MORE THAN 306 CHANNELS in ' condition_mismatch_names{c} ' ' participant_group{pg}]);
    elseif size(Matrix,1) ~= 62 && strcmp(modality_data{mode},'EEG')
        error(['EEG DATA MATRIX HAS MORE THAN 62 CHANNELS in ' condition_mismatch_names{c} ' ' participant_group{pg}]);
    end
    % Select participant labels and ensure they match with third dimension of Matrix
    good_participants_list = {};
    for p = 1:length(participant)
        % Identify patients/controls
        pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
        if ~strcmp(subject_array{pos_subj,2},participant_group{pg})
            continue; % so only include participants that correspond to the group
        end
        
        % If subject does not have any of these markers, don't include it in
        % the sweep count (because it may be a bad subject, etc)
        if ~strcmp(subject_array{pos_subj,3},'needs_forward') && ~strcmp(subject_array{pos_subj,3},'needs_sources') && ~strcmp(subject_array{pos_subj,3},'sources_finished')
            continue;
        end
        
        subject_row = {subject_array{pos_subj,:}}; 
        subject_row = subject_row(~cellfun('isempty', subject_row')); % to avoid empty cells in contains
        bad_signal = find(contains(subject_row,['bad_' modality_data{mode}]));
        if strcmp(modality_data{mode},'BIMODAL') % In case of bimodal, check if any EEG or MEG is bad
            bad_signal_EEG = find(contains(subject_row,['bad_EEG']));
            bad_signal_MEG = find(contains(subject_row,['bad_MEG']));
            % If either one of these is bad, mark bimodal as bad
            if ~isempty(bad_signal_EEG) || ~isempty(bad_signal_MEG)
                bad_signal = 1;
            end
        end
        if ~isempty(bad_signal)
            continue; % to next subject
        end     
        good_participants_list{p} = subject_array{pos_subj,1};
    end
    good_participants_list = good_participants_list(~cellfun('isempty', good_participants_list'));
    if length(good_participants_list) ~= size(Matrix,3)
        error('good_participants number does not match Matrix dimension')
    end
    % Determine position in figure based on iteration
    % (EEG/MEG(2Fig); C/FE; 5cond; 2x5)
    if strcmp(condition_mismatch_names{c},'Standard') && strcmp(participant_group{pg},'C') 
        plot_pos = 1;
    elseif strcmp(condition_mismatch_names{c},'DeviantPitch') && strcmp(participant_group{pg},'C') 
        plot_pos = 2;
    elseif strcmp(condition_mismatch_names{c},'DeviantDuration') && strcmp(participant_group{pg},'C') 
        plot_pos = 3;
    elseif strcmp(condition_mismatch_names{c},'DeviantPitch-Standard') && strcmp(participant_group{pg},'C') 
        plot_pos = 4;
    elseif strcmp(condition_mismatch_names{c},'DeviantDuration-Standard') && strcmp(participant_group{pg},'C') 
        plot_pos = 5;
    elseif strcmp(condition_mismatch_names{c},'Standard') && strcmp(participant_group{pg},'FE') 
        plot_pos = 6;
    elseif strcmp(condition_mismatch_names{c},'DeviantPitch') && strcmp(participant_group{pg},'FE') 
        plot_pos = 7;
    elseif strcmp(condition_mismatch_names{c},'DeviantDuration') && strcmp(participant_group{pg},'FE') 
        plot_pos = 8;
    elseif strcmp(condition_mismatch_names{c},'DeviantPitch-Standard') && strcmp(participant_group{pg},'FE') 
        plot_pos = 9;
    elseif strcmp(condition_mismatch_names{c},'DeviantDuration-Standard') && strcmp(participant_group{pg},'FE') 
        plot_pos = 10;
    
    end
    
    % Now plot
    for gp = 1:length(good_participants_list)
        % If cluster was selected, average amplitudes before plotting
        Amplitudes = Matrix(chan_num,:,gp);
        if size(Amplitudes,1) ~= 1
           Amplitudes = squeeze(mean(Amplitudes,1));
        else
           Amplitudes = squeeze(Amplitudes);
        end
        hold (h(plot_pos),'on')
        plot(h(plot_pos),time_samples, Amplitudes);
        current_max(pos_ylim) = max(Amplitudes);
        current_min(pos_ylim) = min(Amplitudes);
        pos_ylim = pos_ylim + 1; % iterations of participants and conditions
    end
    % Add title (unique of every subplot) and legend (unique of group)
    hold (h(plot_pos),'on')
    title (h(plot_pos),[abreviated_cond_titles{c} ' ' participant_group{pg}]);
    legend(h(plot_pos),good_participants_list);
end

end 
    % Once data is put in figure, adjust format
    y_lim_min = min(current_min);
    y_lim_max = max(current_max);
    for i = 1:num_subplots
        hold (h(i),'on')
        xlim(h(i),[plot_baseline,plot_post]);
        ylim(h(i),[y_lim_min,y_lim_max]);
        plot(h(i),xlim,[0 0], '-k')
        line(h(i),[0 0], [y_lim_min y_lim_max],'Color','black')
    end
    hold (h(6),'on')
    if strcmp(modality_data{mode},'EEG')
        ylabel(h(6),'Amplitude (µV)'); 
    elseif strcmp(modality_data{mode},'MEG')
        ylabel(h(6),'Magnetic field (fT)'); 
    end
    xlabel(h(6),'Time (ms)');
    current_title = ['All subjects_' modality_data{mode}];
    current_title = strrep(current_title,'_',' ');
    suptitle(current_title)
    pom=findobj('type','legend');
    delete(pom);
%     if strcmp(modality_data{mode},'EEG')
%         if size(choice_channel_EEG,2) > 1 % cluster was selected
%             legend(h(1),'Cluster');
%         else
%             legend(h(1),choice_channel_EEG{:});
%         end
%     elseif strcmp(modality_data{mode},'MEG')
%         if size(choice_channel_MEG,2) > 1 % cluster was selected
%             legend(h(1),'Cluster');
%         else
%             legend(h(1),choice_channel_MEG{:});
%         end
%     end
% Save figures if stated
% if save_figures_option == 1
%     % Ensure directory for this GAVR folder exists
%     if ~exist([root_dir '/Figures/' gavr_name],'dir')
%         mkdir([root_dir '/Figures/'], gavr_name);
%     end  
% end
% Patch time windows in gray (shaded areas)
    if shaded_areas == 1
        if strcmp(modality_data{mode},'EEG')
            shaded_areas_scalp_Pitch_MMN = shaded_areas_scalp_Pitch_MMN_EEG;
            shaded_areas_scalp_Dur_MMN = shaded_areas_scalp_Dur_MMN_EEG;
        elseif strcmp(modality_data{mode},'MEG')
            shaded_areas_scalp_Pitch_MMN = shaded_areas_scalp_Pitch_MMN_MEG;
            shaded_areas_scalp_Dur_MMN = shaded_areas_scalp_Dur_MMN_MEG; 
        end
    % Patch Pitch MMN plots
    for i = [4 9] 
        hold (h(i),'on')
        for shad = 1:length(shaded_areas_scalp_Pitch_MMN)
            gray = [0 0 0];
            patch(h(i),shaded_areas_scalp_Pitch_MMN{shad},[min(ylim)*[1 1] max(ylim)*[1 1]],gray,'FaceAlpha', transparency,'EdgeAlpha',0.1,'HandleVisibility','off')
        end
    end
    % Patch Duration MMN plots
    for i = [5 10] 
        hold (h(i),'on')
        for shad = 1:length(shaded_areas_scalp_Dur_MMN)
            gray = [0 0 0];
            patch(h(i),shaded_areas_scalp_Dur_MMN{shad},[min(ylim)*[1 1] max(ylim)*[1 1]],gray,'FaceAlpha', transparency,'EdgeAlpha',0.1,'HandleVisibility','off')
        end
    end
    end
end

% Reset modality_data variable to its original self
modality_data = {'EEG','MEG','BIMODAL'};

%% Plot GAVR scalp (overlay C/FE; EEG/MEG; 5cond; 1Fig 2x5)

% Change modality_data list depending on the choice of sources
modality_data = {'EEG','MEG'};

% Predefine which channels/clusters are going to be used from matrix
load([root_dir '/Scripts/Areas_channels.mat']);
% EEG channels
if size(choice_channel_EEG,2) > 1 % cluster was selected
    channel_pos_EEG = [];
    for i = 1:size(choice_channel_EEG,2) % Find all the indices
        channel_pos_EEG(i) = find(strcmp({Areas.Name}, choice_channel_EEG{i})==1);
    end
else % single channel
   channel_pos_EEG = find(strcmp({Areas.Name}, choice_channel_EEG)==1); 
end
channel_pos_EEG = channel_pos_EEG - 314; % to adjust for 62 channels
% MEG sensors
if size(choice_channel_MEG,2) > 1 % cluster was selected
    channel_pos_MEG = [];
    for i = 1:size(choice_channel_MEG,2) % Find all the indices
        channel_pos_MEG(i) = find(strcmp({Areas.Name}, choice_channel_MEG{i})==1);
    end
else % single sensor
   channel_pos_MEG = find(strcmp({Areas.Name}, choice_channel_MEG)==1); 
end

% New plot 
figure('units','normalized','outerposition',[0 0 1 1]);
num_subplots = 10; % For later
h(1) = subplot(2,5,1); h(2) = subplot(2,5,2); h(3) = subplot(2,5,3); 
h(4) = subplot(2,5,4);h(5) = subplot(2,5,5); h(6) = subplot(2,5,6);
h(7) = subplot(2,5,7);h(8) = subplot(2,5,8); h(9) = subplot(2,5,9);
h(10) = subplot(2,5,10);
for mode = 1:length(modality_data)
    % Define channel number
    eval(['chan_num = channel_pos_' modality_data{mode} ';'])
    current_max = []; current_min = []; pos_ylim = 1; % To adjust ylim later
for c = 1:length(condition_mismatch_names)
    % Determine position in figure based on iteration
    % (overlay C/FE; EEG/MEG; 5cond; 1Fig 2x5)
    if strcmp(condition_mismatch_names{c},'Standard') && strcmp(modality_data{mode},'EEG') 
        plot_pos = 1;
    elseif strcmp(condition_mismatch_names{c},'DeviantPitch') && strcmp(modality_data{mode},'EEG') 
        plot_pos = 2;
    elseif strcmp(condition_mismatch_names{c},'DeviantDuration') && strcmp(modality_data{mode},'EEG') 
        plot_pos = 3;
    elseif strcmp(condition_mismatch_names{c},'DeviantPitch-Standard') && strcmp(modality_data{mode},'EEG')  
        plot_pos = 4;
    elseif strcmp(condition_mismatch_names{c},'DeviantDuration-Standard') && strcmp(modality_data{mode},'EEG') 
        plot_pos = 5;
    elseif strcmp(condition_mismatch_names{c},'Standard') && strcmp(modality_data{mode},'MEG') 
        plot_pos = 6;
    elseif strcmp(condition_mismatch_names{c},'DeviantPitch') && strcmp(modality_data{mode},'MEG')
        plot_pos = 7;
    elseif strcmp(condition_mismatch_names{c},'DeviantDuration') && strcmp(modality_data{mode},'MEG')
        plot_pos = 8;
    elseif strcmp(condition_mismatch_names{c},'DeviantPitch-Standard') && strcmp(modality_data{mode},'MEG')
        plot_pos = 9;
    elseif strcmp(condition_mismatch_names{c},'DeviantDuration-Standard') && strcmp(modality_data{mode},'MEG') 
        plot_pos = 10;
    end
for pg = 1:length(participant_group)
    % Load data matrix
    load([root_dir '/Scalp/' gavr_name '/gavr/' modality_data{mode} '_' condition_mismatch_names{c} '_' participant_group{pg} '.mat']);
    % Ensure number of channels is correct
    if size(F,1) ~= 306 && strcmp(modality_data{mode},'MEG')
        error(['MEG DATA MATRIX HAS MORE THAN 306 CHANNELS in ' condition_mismatch_names{c} ' ' participant_group{pg}]);
    elseif size(F,1) ~= 62 && strcmp(modality_data{mode},'EEG')
        error(['EEG DATA MATRIX HAS MORE THAN 62 CHANNELS in ' condition_mismatch_names{c} ' ' participant_group{pg}]);
    end

    % Load gavr variable
    average = F(chan_num,:);
    % If cluster was selected, average amplitudes before plotting
    if size(average,1) ~= 1
       average = squeeze(mean(average,1));
    else
       average = squeeze(average);
    end
    % Load stdev or stderr
    if dev_GAVR == 1 % standard deviation
        load([root_dir '/Scalp/' gavr_name '/std_dev/' modality_data{mode} '_' condition_mismatch_names{c} '_' participant_group{pg} '.mat']);
        dev = F(chan_num,:);
    elseif dev_GAVR == 2 % standard error
        load([root_dir '/Scalp/' gavr_name '/std_err/' modality_data{mode} '_' condition_mismatch_names{c} '_' participant_group{pg} '.mat']);
        dev = F(chan_num,:);
    end  
    % Set data ready for plot
    curve1 = average + dev;
    curve2 = average - dev;
    time_samples_2 = [time_samples, fliplr(time_samples)];
    inBetween = [curve1, fliplr(curve2)];
    % Also, grab values for ylim later
    current_max(pos_ylim) = max(curve1);
    current_min(pos_ylim) = min(curve2);
    pos_ylim = pos_ylim + 1; % iterations of participants and conditions
    
    % Now plot
    hold (h(plot_pos),'on')
    fill(h(plot_pos),time_samples_2, inBetween, (color_group{pg}/256), 'FaceAlpha', transparency, 'LineStyle', 'none','HandleVisibility','off');
    plot(h(plot_pos),time_samples, average, 'color', (color_group{pg}/256), 'LineWidth', 1.5);
    
    % Add title (unique of every subplot)
    hold (h(plot_pos),'on')
    title (h(plot_pos),[abreviated_cond_titles{c} ' ' modality_data{mode}]);
end
end 
% We set ylim here separataley for EEG and MEG
y_lim_min = min(current_min);
y_lim_max = max(current_max);
if strcmp(modality_data{mode},'EEG')
    plot_iterations = 1:5; % EEG ones
elseif strcmp(modality_data{mode},'MEG')
    plot_iterations = 6:10; % MEG ones
end
for i = plot_iterations
    hold (h(i),'on')
    ylim(h(i),[y_lim_min,y_lim_max]);
    line(h(i),[0 0], [y_lim_min y_lim_max],'Color','black')
end
end
% Once data is put in figure, adjust format
for i = 1:num_subplots
    hold (h(i),'on')
    xlim(h(i),[plot_baseline,plot_post]);
    % legend(h(i),participant_group);
    plot(h(i),xlim,[0 0], '-k')
end
hold (h(1),'on')
ylabel(h(1),'Amplitude (µV)'); 
legend(h(1),participant_group); % If we wanted it in only one
hold (h(6),'on')
ylabel(h(6),'Magnetic field (fT)'); 
xlabel(h(6),'Time (ms)');
current_title = 'Grand average sensor level';
% current_title = strrep(current_title,'_',' ');
suptitle(current_title)

% Patch MMN areas
if shaded_areas == 1
% Patch Pitch MMN plots EEG
for i = [4] 
    hold (h(i),'on')
    for shad = 1:length(shaded_areas_scalp_Pitch_MMN_EEG)
        gray = [0 0 0];
        patch(h(i),shaded_areas_scalp_Pitch_MMN_EEG{shad},[min(ylim)*[1 1] max(ylim)*[1 1]],gray,'FaceAlpha', transparency,'EdgeAlpha',0.1,'HandleVisibility','off')
    end
end
% Patch Duration MMN plots EEG
for i = [5] 
    hold (h(i),'on')
    for shad = 1:length(shaded_areas_scalp_Dur_MMN_EEG)
        gray = [0 0 0];
        patch(h(i),shaded_areas_scalp_Dur_MMN_EEG{shad},[min(ylim)*[1 1] max(ylim)*[1 1]],gray,'FaceAlpha', transparency,'EdgeAlpha',0.1,'HandleVisibility','off')
    end
end
% Patch Pitch MMN plots MEG
for i = [9] 
    hold (h(i),'on')
    for shad = 1:length(shaded_areas_scalp_Pitch_MMN_MEG)
        gray = [0 0 0];
        patch(h(i),shaded_areas_scalp_Pitch_MMN_MEG{shad},[min(ylim)*[1 1] max(ylim)*[1 1]],gray,'FaceAlpha', transparency,'EdgeAlpha',0.1,'HandleVisibility','off')
    end
end
% Patch Duration MMN plots MEG
for i = [10] 
    hold (h(i),'on')
    for shad = 1:length(shaded_areas_scalp_Dur_MMN_MEG)
        gray = [0 0 0];
        patch(h(i),shaded_areas_scalp_Dur_MMN_MEG{shad},[min(ylim)*[1 1] max(ylim)*[1 1]],gray,'FaceAlpha', transparency,'EdgeAlpha',0.1,'HandleVisibility','off')
    end
end
end


% Reset modality_data variable to its original self
modality_data = {'EEG','MEG','BIMODAL'};

%% Plot butterfly source waveforms (MEG only; C/FE; L/R; 2cond; 6 ROIs)
% A1, LBelt, MBelt, PBelt, STSdp, STSvp
% 1FigxROI: 2x4 figure (2condvs2groupvs2hem)

% Change modality_data list depending on the choice of sources
if source_analysis == 1
    modality_data = {'EEG'};
elseif source_analysis == 2
    modality_data = {'MEG'};
elseif source_analysis == 3
    modality_data = {'BIMODAL'};
elseif source_analysis == 4
    modality_data = {'EEG','MEG','BIMODAL'};
end

for mode = 1:length(modality_data)
for roi = 1:length(ROIs)
    % New plot 
    figure('units','normalized','outerposition',[0 0 1 1]);
    num_subplots = 8; % For later
    h(1) = subplot(2,4,1); h(2) = subplot(2,4,2); h(3) = subplot(2,4,3); 
    h(4) = subplot(2,4,4);h(5) = subplot(2,4,5); h(6) = subplot(2,4,6);
    h(7) = subplot(2,4,7);h(8) = subplot(2,4,8);
    current_max = []; current_min = []; pos_ylim = 1; % To adjust ylim later
for pg = 1:length(participant_group)
for hem = 1:length(hemisphere)
    pos_source = find(strcmp(scout_positions,[ROIs{roi} '_' hemisphere{hem}]));
for c = 4:length(condition_mismatch_names) % ONLY MMN CONDITIONS HERE!!!!!!!!!!!!!!!
    load([root_dir '/Sources/Waveforms/' modality_data{mode} '/' gavr_name '_' normalization_option_tag '/Matrix_' scout_function_string '_' condition_mismatch_names{c} '_' participant_group{pg} '.mat'])
    % Determine position in matrix
    if length(scout_positions) ~= size(sources,1)
        error('Scout names described in this section do not match Matrix dimensions');
    end
    % Select participant labels and ensure they match with third dimension of Matrix
    good_participants_list = {};
    for p = 1:length(participant)
        % Identify patients/controls
        pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
        if ~strcmp(subject_array{pos_subj,2},participant_group{pg})
            continue; % so only include participants that correspond to the group
        end
        
        % If subject does not have any of these markers, don't include it in
        % the sweep count (because it may be a bad subject, etc)
        if ~strcmp(subject_array{pos_subj,3},'needs_forward') && ~strcmp(subject_array{pos_subj,3},'needs_sources') && ~strcmp(subject_array{pos_subj,3},'sources_finished')
            continue;
        end
        
        subject_row = {subject_array{pos_subj,:}}; 
        subject_row = subject_row(~cellfun('isempty', subject_row')); % to avoid empty cells in contains
        bad_signal = find(contains(subject_row,['bad_' modality_data{mode}]));
        if strcmp(modality_data{mode},'BIMODAL') % In case of bimodal, check if any EEG or MEG is bad
            bad_signal_EEG = find(contains(subject_row,['bad_EEG']));
            bad_signal_MEG = find(contains(subject_row,['bad_MEG']));
            % If either one of these is bad, mark bimodal as bad
            if ~isempty(bad_signal_EEG) || ~isempty(bad_signal_MEG)
                bad_signal = 1;
            end
        end
        if ~isempty(bad_signal)
            continue; % to next subject
        end     
        % Cov problems
        no_cov_sweeps = find(contains(subject_row,['insuf_cov_' modality_data{mode}]));
        if ~isempty(no_cov_sweeps)
            continue; % to next modality
        end
        % Sources-or-anatomy-specific problems
        if ~strcmp(subject_array{pos_subj,10},'corregistration_done') || ~strcmp(subject_array{pos_subj,11},'FEM_anatomy_ready')...
          || strcmp(subject_array{pos_subj,12},'No_HCP') || ~strcmp(subject_array{pos_subj,13},'EEG_projected')
            continue; % on to next subject
        end   
        good_participants_list{p} = subject_array{pos_subj,1};
    end
    good_participants_list = good_participants_list(~cellfun('isempty', good_participants_list'));
    if length(good_participants_list) ~= size(sources,3)
        error('good_participants number does not match Matrix dimension')
    end
    % Determine position in figure based on iteration
    % 4x4 figure (2condvs2groupvs2hem)
    if strcmp(condition_mismatch_names{c},'DeviantPitch-Standard') && strcmp(hemisphere{hem},'L') && strcmp(participant_group{pg},'C') 
        plot_pos = 1;
    elseif strcmp(condition_mismatch_names{c},'DeviantPitch-Standard') && strcmp(hemisphere{hem},'R') && strcmp(participant_group{pg},'C') 
        plot_pos = 2;
    elseif strcmp(condition_mismatch_names{c},'DeviantDuration-Standard') && strcmp(hemisphere{hem},'L') && strcmp(participant_group{pg},'C') 
        plot_pos = 3;
    elseif strcmp(condition_mismatch_names{c},'DeviantDuration-Standard') && strcmp(hemisphere{hem},'R') && strcmp(participant_group{pg},'C') 
        plot_pos = 4;
    elseif strcmp(condition_mismatch_names{c},'DeviantPitch-Standard') && strcmp(hemisphere{hem},'L') && strcmp(participant_group{pg},'FE') 
        plot_pos = 5;
    elseif strcmp(condition_mismatch_names{c},'DeviantPitch-Standard') && strcmp(hemisphere{hem},'R') && strcmp(participant_group{pg},'FE')
        plot_pos = 6;
    elseif strcmp(condition_mismatch_names{c},'DeviantDuration-Standard') && strcmp(hemisphere{hem},'L') && strcmp(participant_group{pg},'FE') 
        plot_pos = 7;
    elseif strcmp(condition_mismatch_names{c},'DeviantDuration-Standard') && strcmp(hemisphere{hem},'R') && strcmp(participant_group{pg},'FE')
        plot_pos = 8;
    end
    % Now plot
    for gp = 1:length(good_participants_list)
        hold (h(plot_pos),'on')
        plot(h(plot_pos),time_samples, sources(pos_source,:,gp));
        current_max(pos_ylim) = max(sources(pos_source,:,gp));
        current_min(pos_ylim) = min(sources(pos_source,:,gp));
        pos_ylim = pos_ylim + 1; % iterations of participants and conditions
    end
    % Add title (unique of every subplot)
    hold (h(plot_pos),'on')
    title (h(plot_pos),[abreviated_cond_titles{c} ' ' participant_group{pg} ' ' hemisphere{hem}]);
    legend(h(plot_pos),good_participants_list);
end
end
end 
    % Once data is put in figure, adjust format
    y_lim_min = min(current_min);
    y_lim_max = max(current_max);
    for i = 1:num_subplots
        hold (h(i),'on')
        xlim(h(i),[plot_baseline,plot_post]);
        ylim(h(i),[y_lim_min,y_lim_max]);
        plot(h(i),xlim,[0 0], '-k')
        line(h(i),[0 0], [0 y_lim_max],'Color','black')
    end
    hold (h(5),'on')
    ylabel(h(5),'Relative units'); 
    xlabel(h(5),'Time (ms)');
    current_title = ['All subjects_' ROIs{roi} '_' modality_data{mode}];
    current_title = strrep(current_title,'_',' ');
    suptitle(current_title)
    pom=findobj('type','legend');
    delete(pom);
    
    if shaded_areas == 1
    % Patch Pitch MMN plots
    for i = [1 5] 
        hold (h(i),'on')
        for shad = 1:length(shaded_areas_source_left_Pitch_MMN)
            gray = [0 0 0];
            patch(h(i),shaded_areas_source_left_Pitch_MMN{shad},[min(ylim)*[1 1] max(ylim)*[1 1]],gray,'FaceAlpha', transparency,'EdgeAlpha',0.1,'HandleVisibility','off')
        end
    end
    for i = [2 6] 
        hold (h(i),'on')
        for shad = 1:length(shaded_areas_source_right_Pitch_MMN)
            gray = [0 0 0];
            patch(h(i),shaded_areas_source_right_Pitch_MMN{shad},[min(ylim)*[1 1] max(ylim)*[1 1]],gray,'FaceAlpha', transparency,'EdgeAlpha',0.1,'HandleVisibility','off')
        end
    end
    % Patch Duration MMN plots
    for i = [3 7] 
        hold (h(i),'on')
        for shad = 1:length(shaded_areas_source_left_Dur_MMN)
            gray = [0 0 0];
            patch(h(i),shaded_areas_source_left_Dur_MMN{shad},[min(ylim)*[1 1] max(ylim)*[1 1]],gray,'FaceAlpha', transparency,'EdgeAlpha',0.1,'HandleVisibility','off')
        end
    end
    for i = [4 8] 
        hold (h(i),'on')
        for shad = 1:length(shaded_areas_source_left_Dur_MMN)
            gray = [0 0 0];
            patch(h(i),shaded_areas_source_left_Dur_MMN{shad},[min(ylim)*[1 1] max(ylim)*[1 1]],gray,'FaceAlpha', transparency,'EdgeAlpha',0.1,'HandleVisibility','off')
        end
    end
    end
end
end

% Reset modality_data variable to its original self
modality_data = {'EEG','MEG','BIMODAL'};

%% Plot GAVR source waveforms (MEG only; C/FE; L/R; 2cond; 6 ROIs)
% A1, LBelt, MBelt, PBelt, STSdp, STSvp
% 1) overlay C/FE, 2 condvs2hem; 1FigxROI: 6 2x2 figures
% 2) overlay L/R, 2condvs2group; 1FigxROI: 6 2x2 figures
% 3) overlay 6ROIs, 2condvs2groupvs2hem; 1 4x4 figure

% Change modality_data list depending on the choice of sources
if source_analysis == 1
    modality_data = {'EEG'};
elseif source_analysis == 2
    modality_data = {'MEG'};
elseif source_analysis == 3
    modality_data = {'BIMODAL'};
elseif source_analysis == 4
    modality_data = {'EEG','MEG','BIMODAL'};
end

% 1) overlay C/FE, 2 condvs2hem; 1FigxROI: 6 2x2 figures
for mode = 1:length(modality_data)
for roi = 1:length(ROIs)
    % New plot 
    figure; % ('units','normalized','outerposition',[0 0 1 1]);
    num_subplots = 4; % For later
    h(1) = subplot(2,2,1); h(2) = subplot(2,2,2); 
    h(3) = subplot(2,2,3); h(4) = subplot(2,2,4);
    current_max = []; current_min = []; pos_ylim = 1; % To adjust ylim later
for hem = 1:length(hemisphere)
    pos_source = find(strcmp(scout_positions,[ROIs{roi} '_' hemisphere{hem}]));
for c = 4:length(condition_mismatch_names) % ONLY MMN CONDITIONS HERE!!!!!!!!!!!!!!!
    % Determine position in figure based on iteration
    % overlay C/FE, 2 condvs2hem; 1FigxROI: 6 2x2 figures
    if strcmp(condition_mismatch_names{c},'DeviantPitch-Standard') && strcmp(hemisphere{hem},'L')
        plot_pos = 1;
    elseif strcmp(condition_mismatch_names{c},'DeviantPitch-Standard') && strcmp(hemisphere{hem},'R')
        plot_pos = 2;
    elseif strcmp(condition_mismatch_names{c},'DeviantDuration-Standard') && strcmp(hemisphere{hem},'L')
        plot_pos = 3;
    elseif strcmp(condition_mismatch_names{c},'DeviantDuration-Standard') && strcmp(hemisphere{hem},'R')
        plot_pos = 4;
    end  
for pg = 1:length(participant_group)
    load([root_dir '/Sources/Waveforms/' modality_data{mode} '/' gavr_name '_' normalization_option_tag '/gavr/' scout_function_string '_' condition_mismatch_names{c} '_' participant_group{pg} '.mat'])
    % Determine if scout position is going to be accurate
    if length(scout_positions) ~= size(F,1)
        error('Scout names described in this section do not match variable dimensions');
    end
    % Load gavr variable
    average = F(pos_source,:);
    % Load stdev or stderr
    if dev_GAVR == 1 % standard deviation
        load([root_dir '/Sources/Waveforms/' modality_data{mode} '/' gavr_name '_' normalization_option_tag '/std_dev/' scout_function_string '_' condition_mismatch_names{c} '_' participant_group{pg} '.mat']);
        dev = F(pos_source,:);
    elseif dev_GAVR == 2 % standard error
        load([root_dir '/Sources/Waveforms/' modality_data{mode} '/' gavr_name '_' normalization_option_tag '/std_err/' scout_function_string '_' condition_mismatch_names{c} '_' participant_group{pg} '.mat']);
        dev = F(pos_source,:);
    end  
    % Set data ready for plot
    curve1 = average + dev;
    curve2 = average - dev;
    time_samples_2 = [time_samples, fliplr(time_samples)];
    inBetween = [curve1, fliplr(curve2)];
    % Also, grab values for ylim later
    current_max(pos_ylim) = max(curve1);
    current_min(pos_ylim) = min(curve2);
    pos_ylim = pos_ylim + 1; % iterations of participants and conditions
    
    % Now plot
    hold (h(plot_pos),'on')
    fill(h(plot_pos),time_samples_2, inBetween, (color_group{pg}/256), 'FaceAlpha', transparency, 'LineStyle', 'none','HandleVisibility','off');
    plot(h(plot_pos),time_samples, average, 'color', (color_group{pg}/256), 'LineWidth', 1.5);
end
    % Add title (unique of every subplot)
    hold (h(plot_pos),'on')
    title (h(plot_pos),[abreviated_cond_titles{c} ' ' hemisphere{hem}])  
end
end 
    % Once data is put in figure, adjust format
    y_lim_min = min(current_min);
    y_lim_max = max(current_max);
    for i = 1:num_subplots
        hold (h(i),'on')
        xlim(h(i),[plot_baseline,plot_post]);
        ylim(h(i),[y_lim_min,y_lim_max]);
        % legend(h(i),participant_group);
        plot(h(i),xlim,[0 0], '-k','HandleVisibility','off')
        line(h(i),[0 0], [0 y_lim_max],'Color','black','HandleVisibility','off')
    end
    hold (h(3),'on')
    ylabel(h(3),'Relative units'); 
    xlabel(h(3),'Time (ms)');
    hold (h(1),'on')
    legend(h(1),participant_group); % If we wanted it in only one
    current_title = [ROIs{roi} '_' modality_data{mode}];
    current_title = strrep(current_title,'_',' ');
    suptitle(current_title)
    
    if shaded_areas == 1
    % Patch Pitch MMN plots
    for i = 1 
        hold (h(i),'on')
        for shad = 1:length(shaded_areas_source_left_Pitch_MMN)
            gray = [0 0 0];
            patch(h(i),shaded_areas_source_left_Pitch_MMN{shad},[min(ylim)*[1 1] max(ylim)*[1 1]],gray,'FaceAlpha', transparency,'EdgeAlpha',0.1,'HandleVisibility','off')
        end
    end
    for i = 2 
        hold (h(i),'on')
        for shad = 1:length(shaded_areas_source_right_Pitch_MMN)
            gray = [0 0 0];
            patch(h(i),shaded_areas_source_right_Pitch_MMN{shad},[min(ylim)*[1 1] max(ylim)*[1 1]],gray,'FaceAlpha', transparency,'EdgeAlpha',0.1,'HandleVisibility','off')
        end
    end
    % Patch Duration MMN plots
    for i = 3 
        hold (h(i),'on')
        for shad = 1:length(shaded_areas_source_left_Dur_MMN)
            gray = [0 0 0];
            patch(h(i),shaded_areas_source_left_Dur_MMN{shad},[min(ylim)*[1 1] max(ylim)*[1 1]],gray,'FaceAlpha', transparency,'EdgeAlpha',0.1,'HandleVisibility','off')
        end
    end
    for i = 4 
        hold (h(i),'on')
        for shad = 1:length(shaded_areas_source_right_Dur_MMN)
            gray = [0 0 0];
            patch(h(i),shaded_areas_source_right_Dur_MMN{shad},[min(ylim)*[1 1] max(ylim)*[1 1]],gray,'FaceAlpha', transparency,'EdgeAlpha',0.1,'HandleVisibility','off')
        end
    end
    end
end
end

% 2) overlay L/R, 2condvs2group; 1FigxROI: 6 2x2 figures
for mode = 1:length(modality_data)
for roi = 1:length(ROIs)
    % New plot 
    figure; % ('units','normalized','outerposition',[0 0 1 1]);
    num_subplots = 4; % For later
    h(1) = subplot(2,2,1); h(2) = subplot(2,2,2); 
    h(3) = subplot(2,2,3); h(4) = subplot(2,2,4);
    current_max = []; current_min = []; pos_ylim = 1; % To adjust ylim later
for c = 4:length(condition_mismatch_names) % ONLY MMN CONDITIONS HERE!!!!!!!!!!!!!!!
for pg = 1:length(participant_group)
    % Determine position in figure based on iteration
    % overlay L/R, 2condvs2group; 1FigxROI: 6 2x2 figures
    if strcmp(condition_mismatch_names{c},'DeviantPitch-Standard') && strcmp(participant_group{pg},'C')
        plot_pos = 1;
    elseif strcmp(condition_mismatch_names{c},'DeviantDuration-Standard') && strcmp(participant_group{pg},'C')
        plot_pos = 2;
    elseif strcmp(condition_mismatch_names{c},'DeviantPitch-Standard') && strcmp(participant_group{pg},'FE')
        plot_pos = 3;
    elseif strcmp(condition_mismatch_names{c},'DeviantDuration-Standard') && strcmp(participant_group{pg},'FE')
        plot_pos = 4;
    end    
for hem = 1:length(hemisphere)
    pos_source = find(strcmp(scout_positions,[ROIs{roi} '_' hemisphere{hem}]));
    load([root_dir '/Sources/Waveforms/' modality_data{mode} '/' gavr_name '_' normalization_option_tag '/gavr/' scout_function_string '_' condition_mismatch_names{c} '_' participant_group{pg} '.mat'])
    % Determine if scout position is going to be accurate
    if length(scout_positions) ~= size(F,1)
        error('Scout names described in this section do not match variable dimensions');
    end
    % Load gavr variable
    average = F(pos_source,:);
    % Load stdev or stderr
    if dev_GAVR == 1 % standard deviation
        load([root_dir '/Sources/Waveforms/' modality_data{mode} '/' gavr_name '_' normalization_option_tag '/std_dev/' scout_function_string '_' condition_mismatch_names{c} '_' participant_group{pg} '.mat']);
        dev = F(pos_source,:);
    elseif dev_GAVR == 2 % standard error
        load([root_dir '/Sources/Waveforms/' modality_data{mode} '/' gavr_name '_' normalization_option_tag '/std_err/' scout_function_string '_' condition_mismatch_names{c} '_' participant_group{pg} '.mat']);
        dev = F(pos_source,:);
    end  
    % Set data ready for plot
    curve1 = average + dev;
    curve2 = average - dev;
    time_samples_2 = [time_samples, fliplr(time_samples)];
    inBetween = [curve1, fliplr(curve2)];
    % Also, grab values for ylim later
    current_max(pos_ylim) = max(curve1);
    current_min(pos_ylim) = min(curve2);
    pos_ylim = pos_ylim + 1; % iterations of participants and conditions
    
    % Now plot
    hold (h(plot_pos),'on')
    fill(h(plot_pos),time_samples_2, inBetween, (color_hem{hem}/256), 'FaceAlpha', transparency, 'LineStyle', 'none','HandleVisibility','off');
    plot(h(plot_pos),time_samples, average, 'color', (color_hem{hem}/256), 'LineWidth', 1.5);
end
    % Add title (unique of every subplot)
    hold (h(plot_pos),'on')
    title (h(plot_pos),[abreviated_cond_titles{c} ' ' participant_group{pg}])
end
end 
    % Once data is put in figure, adjust format
    y_lim_min = min(current_min);
    y_lim_max = max(current_max);
    for i = 1:num_subplots
        hold (h(i),'on')
        xlim(h(i),[plot_baseline,plot_post]);
        ylim(h(i),[y_lim_min,y_lim_max]);
        % legend(h(i),participant_group);
        plot(h(i),xlim,[0 0], '-k','HandleVisibility','off')
        line(h(i),[0 0], [0 y_lim_max],'Color','black','HandleVisibility','off')
    end
    hold (h(3),'on')
    ylabel(h(3),'Relative units'); 
    xlabel(h(3),'Time (ms)');
    hold (h(1),'on')
    legend(h(1),hemisphere); % If we wanted it in only one
    current_title = [ROIs{roi} '_' modality_data{mode}];
    current_title = strrep(current_title,'_',' ');
    suptitle(current_title)
end
end

% 3) overlay 6ROIs, 2condvs2groupvs2hem; 1 4x4 figure
for mode = 1:length(modality_data)
    % New plot 
    figure; % ('units','normalized','outerposition',[0 0 1 1]);
    num_subplots = 8; % For later
    h(1) = subplot(2,4,1); h(2) = subplot(2,4,2); h(3) = subplot(2,4,3); 
    h(4) = subplot(2,4,4);h(5) = subplot(2,4,5); h(6) = subplot(2,4,6);
    h(7) = subplot(2,4,7);h(8) = subplot(2,4,8);
    current_max = []; current_min = []; pos_ylim = 1; % To adjust ylim later
for c = 4:length(condition_mismatch_names) % ONLY MMN CONDITIONS HERE!!!!!!!!!!!!!!!
for pg = 1:length(participant_group)
for hem = 1:length(hemisphere)
    % Determine position in figure based on iteration
    % 4x4 figure (2condvs2groupvs2hem)
    if strcmp(condition_mismatch_names{c},'DeviantPitch-Standard') && strcmp(hemisphere{hem},'L') && strcmp(participant_group{pg},'C') 
        plot_pos = 1;
    elseif strcmp(condition_mismatch_names{c},'DeviantPitch-Standard') && strcmp(hemisphere{hem},'R') && strcmp(participant_group{pg},'C') 
        plot_pos = 2;
    elseif strcmp(condition_mismatch_names{c},'DeviantDuration-Standard') && strcmp(hemisphere{hem},'L') && strcmp(participant_group{pg},'C') 
        plot_pos = 3;
    elseif strcmp(condition_mismatch_names{c},'DeviantDuration-Standard') && strcmp(hemisphere{hem},'R') && strcmp(participant_group{pg},'C') 
        plot_pos = 4;
    elseif strcmp(condition_mismatch_names{c},'DeviantPitch-Standard') && strcmp(hemisphere{hem},'L') && strcmp(participant_group{pg},'FE') 
        plot_pos = 5;
    elseif strcmp(condition_mismatch_names{c},'DeviantPitch-Standard') && strcmp(hemisphere{hem},'R') && strcmp(participant_group{pg},'FE')
        plot_pos = 6;
    elseif strcmp(condition_mismatch_names{c},'DeviantDuration-Standard') && strcmp(hemisphere{hem},'L') && strcmp(participant_group{pg},'FE') 
        plot_pos = 7;
    elseif strcmp(condition_mismatch_names{c},'DeviantDuration-Standard') && strcmp(hemisphere{hem},'R') && strcmp(participant_group{pg},'FE')
        plot_pos = 8;
    end
for roi = 1:length(ROIs)
    load([root_dir '/Sources/Waveforms/' modality_data{mode} '/' gavr_name '_' normalization_option_tag '/gavr/' scout_function_string '_' condition_mismatch_names{c} '_' participant_group{pg} '.mat'])
    % Determine if scout position is going to be accurate
    if length(scout_positions) ~= size(F,1)
        error('Scout names described in this section do not match variable dimensions');
    end
    pos_source = find(strcmp(scout_positions,[ROIs{roi} '_' hemisphere{hem}])); 
    % Load gavr variable
    average = F(pos_source,:);
    % Load stdev or stderr
    if dev_GAVR == 1 % standard deviation
        load([root_dir '/Sources/Waveforms/' modality_data{mode} '/' gavr_name '_' normalization_option_tag '/std_dev/' scout_function_string '_' condition_mismatch_names{c} '_' participant_group{pg} '.mat']);
        dev = F(pos_source,:);
    elseif dev_GAVR == 2 % standard error
        load([root_dir '/Sources/Waveforms/' modality_data{mode} '/' gavr_name '_' normalization_option_tag '/std_err/' scout_function_string '_' condition_mismatch_names{c} '_' participant_group{pg} '.mat']);
        dev = F(pos_source,:);
    end 

    % Set data ready for plot
    curve1 = average + dev;
    curve2 = average - dev;
    time_samples_2 = [time_samples, fliplr(time_samples)];
    inBetween = [curve1, fliplr(curve2)];
    % Also, grab values for ylim later
    current_max(pos_ylim) = max(curve1);
    current_min(pos_ylim) = min(curve2);
    pos_ylim = pos_ylim + 1; % iterations of participants and conditions
    
    % Now plot
    hold (h(plot_pos),'on')
    fill(h(plot_pos),time_samples_2, inBetween, (color_roi{roi}/256), 'FaceAlpha', transparency, 'LineStyle', 'none','HandleVisibility','off');
    plot(h(plot_pos),time_samples, average, 'color', (color_roi{roi}/256), 'LineWidth', 1.5);
end
    % Add title (unique of every subplot)
    hold (h(plot_pos),'on')
    title (h(plot_pos),[abreviated_cond_titles{c} ' ' participant_group{pg} ' ' hemisphere{hem}])
end
end 
    % Once data is put in figure, adjust format
    y_lim_min = min(current_min);
    y_lim_max = max(current_max);
    for i = 1:num_subplots
        hold (h(i),'on')
        xlim(h(i),[plot_baseline,plot_post]);
        ylim(h(i),[y_lim_min,y_lim_max]);
        % legend(h(i),participant_group);
        plot(h(i),xlim,[0 0], '-k','HandleVisibility','off')
        line(h(i),[0 0], [0 y_lim_max],'Color','black','HandleVisibility','off')
    end
    hold (h(5),'on')
    ylabel(h(5),'Relative units'); 
    xlabel(h(5),'Time (ms)');
    hold (h(1),'on')
    legend(h(1),ROIs); % If we wanted it in only one
    current_title = ['Averaged sources_' modality_data{mode}];
    current_title = strrep(current_title,'_',' ');
    suptitle(current_title)
end

if shaded_areas == 1
% Patch Pitch MMN plots
for i = [1 5] 
    hold (h(i),'on')
    for shad = 1:length(shaded_areas_source_left_Pitch_MMN)
        gray = [0 0 0];
        patch(h(i),shaded_areas_source_left_Pitch_MMN{shad},[min(ylim)*[1 1] max(ylim)*[1 1]],gray,'FaceAlpha', transparency,'EdgeAlpha',0.1,'HandleVisibility','off')
    end
end
for i = [2 6] 
    hold (h(i),'on')
    for shad = 1:length(shaded_areas_source_right_Pitch_MMN)
        gray = [0 0 0];
        patch(h(i),shaded_areas_source_right_Pitch_MMN{shad},[min(ylim)*[1 1] max(ylim)*[1 1]],gray,'FaceAlpha', transparency,'EdgeAlpha',0.1,'HandleVisibility','off')
    end
end
% Patch Duration MMN plots
for i = [3 7] 
    hold (h(i),'on')
    for shad = 1:length(shaded_areas_source_left_Dur_MMN)
        gray = [0 0 0];
        patch(h(i),shaded_areas_source_left_Dur_MMN{shad},[min(ylim)*[1 1] max(ylim)*[1 1]],gray,'FaceAlpha', transparency,'EdgeAlpha',0.1,'HandleVisibility','off')
    end
end
for i = [4 8] 
    hold (h(i),'on')
    for shad = 1:length(shaded_areas_source_right_Dur_MMN)
        gray = [0 0 0];
        patch(h(i),shaded_areas_source_right_Dur_MMN{shad},[min(ylim)*[1 1] max(ylim)*[1 1]],gray,'FaceAlpha', transparency,'EdgeAlpha',0.1,'HandleVisibility','off')
    end
end
end

end

% Reset modality_data variable to its original self
modality_data = {'EEG','MEG','BIMODAL'};

