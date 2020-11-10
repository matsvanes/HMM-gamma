% Preprocessing script for TMS-MEG data from CZ and CS. This scripts is
% adapted from CZ's MEG_Gamma_HandOver.m. Sections 1-4 have been executed
% by CZ and the resulting data are the starting point for MvE
%% DEFINE JOBS
if ~exist('run','var') || isempty(run)
  run.preproc.filt        = 0;
  run.preproc.ica_run     = 0;
  run.preproc.ica_id      = 0;
  run.preproc.trigger     = 0;
  
  run.preproc.copy        = 1;
  run.preproc.normalise   = 1;
  run.preproc.whiten      = 1;
  run.preproc.bpfilt      = 0;
  run.preproc.coreg       = 1;
  run.preproc.check_coreg = 0;
  run.preproc.parcel      = 0;
  run.preproc.beamform    = 0;
  run.preproc.epoch       = 0;
  
  
  run.keep = 0;
  run.preproc
  disp(['%%% Check if correct jobs are defined! %%%'])
end

%% DEFINE PATHs
if isdir('/ohba/pi/mwoolrich/')
  PATH_BASE = '/ohba/pi/mwoolrich/';
  PATH_ORIG = [PATH_BASE, 'datasets/CZ_Gamma/MEG_Gamma_HandOver/'];
  PATH_BASE = [PATH_BASE, 'mvanes/'];
  islocal = 0;
else
  PATH_BASE = '/Volumes/T5_OHBA/'; % if on a local desktop
  warning('RAWPATH not accessible')
  PATH_ORIG = [];
  islocal = 1;
end
PATH_ANALYSIS = [PATH_BASE, 'analysis/HMM-gamma/'];

PATH_ORIGDATA = [PATH_ORIG, 'data/'];
PATH_RT = [PATH_ORIG 'data/rt'];
PATH_TF = [PATH_ANALYSIS, 'TF/'];
PATH_HMM = [PATH_ANALYSIS 'HMM/'];
PATH_DATA =  [PATH_ANALYSIS 'data/'];
PATH_SCRIPT = [PATH_BASE, 'scripts/HMM-gamma/'];

if islocal
  files=dir([PATH_DATA 'fd_*.mat']);
else
  files=dir([PATH_ORIGDATA 'fd_*.mat']);
end

subinfo;
prefix = 'fd';
%% PARAMS
% subjects
if ~exist('subs', 'var'), subs = 1:length(files); end

% Timing
time_epoch  = [-2000 2000];
time_bl     = [-1 -0.5];
time_peri   = [0 0.5];
time_pre    = [-0.5 0];
time_post   = [1 1.5];

%% 1. DOWNSAMPLE & FILTER
if run.preproc.filt
  files=dir([PATH_ORIGDATA, 'case*_tsss.mat']);
  for s = subs
    meg_tms_preproc01_df(PATH_ORIGDATA,files(s).name);
  end
end

%% 2. ICA
if run.preproc.ica_run || run.preproc.ica_id
  files=dir([PATH_ORIGDATA 'fd_*.mat']);
  for s = subs
    D = spm_eeg_load([PATH_ORIGDATA files(s).name]);
    if run.preproc.ica_run
      D = osl_detect_artefacts(D,'modalities',{'MEGMAG','MEGPLANAR'},'badchannels',true,'badtime',true,'dummy_epoch_tsize',1);
      D = osl_africa(D,'artefact_channels',{'EOG','ECG'},'used_maxfilter',true,'do_ica',true,'do_ident','auto','do_remove',true,'precompute_topos',true);
    elseif run.preproc.ica_id
      D.montage('switch',1);
      D = osl_africa(D,'artefact_channels',{'EOG','ECG'},'do_ica',false','do_ident','manual','do_remove',true');
    end
    D.save
  end
  % Montage 1 = automatic IC identification
  % Montage 2 = manual IC identification
end

%% 3. TRIGGER
if run.preproc.trigger
  files=dir([PATH_ORIGDATA 'fd_*.mat']);
  
  for s = subs
    load([PATH_RT 'emgf' files(s).name '/Timings.mat']);
    
    D = spm_eeg_load([PATH_ORIGDATA files(s).name]);
    Events = D.events;
    % remove old events (from previous iterations of this script)
    to_be_deleted=[];
    for e=1:length(Events)
      if isnumeric(Events(e).value) && (Events(e).value==66 || Events(e).value==99)
        to_be_deleted(end+1)=e;
      end
    end
    Events(to_be_deleted)=[];
    
    % add new move_onset events
    Events_new_start = struct([]);
    for j = 1:length(move_on_eve)
      Events_new_start(end+1).type   = 'Move_on';
      Events_new_start(end).value    = 66;
      Events_new_start(end).time     = move_on_eve(j);
      Events_new_start(end).duration = [];
      Events_new_start(end).offset   = 0;
    end
    if ~isempty(Events_new_start)
      Events = [Events(:); Events_new_start(:)];
    end
    
    % Merge new events with previous
    D = D.events(1,Events);
    D.save
  end %subj
end %if

%% 4. COPY
if run.preproc.copy
  %{
  % copy data
    if isdir(PATH_ORIGDATA)
        rmdir(PATH_ORIGDATA,'s');
    end
    mkdir(PATH_ORIGDATA);
    
    files=dir([PATH_ORIGDATA 'fd_*.mat']);
    for s = subs
        source = fullfile(PATH_ORIGDATA,files(s).name);
        copyfile(source,PATH_ORIGDATA);
    end
end
  %}
  for     s = subs
    files=dir([PATH_ORIGDATA sprintf('%s_*%s*.mat', prefix, sub(s).id)]);
    S           = [];
    S.D         = [PATH_ORIGDATA files.name];
    D           = spm_eeg_load(S.D);
    D = data_to_spm_object(D, D.montage('switch',2), fullfile(PATH_DATA,D.fname));
    D.save
  end
end

%% OPTIONAL - FILTER
if run.preproc.bpfilt
  for s=subs
    files=dir([PATH_DATA sprintf('%s_*%s*.mat', prefix, sub(s).id)]);
    S = [];
    S.D = [PATH_DATA files.name];
    S.band = 'bandpass';
    S.freq = [60 90];
    S.prefix = 'f';
    D = spm_eeg_filter(S);
    D.save;
    if ~keep, delete(S.D);  end
  end
  prefix = ['f' prefix];
end

%% 5. NORMALISE
if run.preproc.normalise
  
  for s = subs
    files=dir([PATH_DATA sprintf('%s_*%s*.mat', prefix, sub(s).id)]);
    
    D = spm_eeg_load([PATH_DATA files.name]);
    
    S               = [];
    S.D             = D;
    S.modalities    = {'MEGMAG','MEGPLANAR'};
    S.do_plots      = 0;
    S.samples2use   = good_samples(D,D.indchantype(S.modalities,'GOOD'));
    S.trials        = 1;
    S.pca_dim       = 99;
    S.force_pca_dim = 0;
    S.normalise_method = 'min_eig';
    D = normalise_sensor_data(S);
    tmp1 = D.chantype;
    tmp2 = D.chanlabels;
    tmp3 = D.montage('getmontage',1);
    
    
    D = data_to_spm_object(D.montage('switch',0), D.montage('switch',1), fullfile(PATH_DATA,D.fname));
    chanunits = {tmp3.channels(:).units};
    sensors = D.sensors('MEG');
    sensors.chanunit(:) = {'T/mm'};% = chanunits; %
    D = D.sensors('MEG',sensors);
    D=units(D,1:306, 'T/mm');
    D = chantype(D,1:D.nchannels,tmp1);
    D = chanlabels(D,1:D.nchannels,tmp2);
    D.save
  end
end

%% 6. SPECTRAL WHITENING (Demanuele et al)
% replace orig data with whitened data in D obj (problem: osl_change_spm_eeg_data does not work properly on data
% with online montag (likely only used after beamforming)
if run.preproc.whiten
  for s = subs
    files=dir([PATH_DATA sprintf('%s_*%s*.mat', prefix, sub(s).id)]);
    D = spm_eeg_load([PATH_DATA, files.name]);
    
    % Copy data file from montag 0 containing original sensor info
    S_new           = [];
    S_new.outfile   = [PATH_DATA, 'w',files.name];
    S_new.D         = D;
    D_new           = spm_eeg_copy(S_new);
    
    % do whitening
    ind_chan = D_new.indchantype('MEGANY');
    data = D_new(ind_chan,:,:);
    data = [diff(data,1,2),zeros(length(ind_chan),1)]./(1/D.fsample);
    
    % put new data in D obj.
    D_new(ind_chan,:,:) = data;
    
    % Write chages to disk
    D_new.save;
  end
  prefix = ['w', prefix];
end


%% 7. COREGISTRATION
if run.preproc.coreg
  for s = subs
    files=dir([PATH_DATA sprintf('%s_*%s*.mat', prefix, sub(s).id)]);
    D = spm_eeg_load([PATH_DATA files.name]);
    
    % Coregistration (see comment in Quinn Frontiers script regarding use headshape, but also consider that here a template is used)
    coreg_set = struct;
    coreg_set.D = D;
    coreg_set.mri = []; %'/home/mvanes/software/osl/std_masks/MNI152_T1_8mm_brain.nii.gz'];
    coreg_set.useheadshape = false;
    coreg_set.forward_meg = 'Single Shell';
    coreg_set.use_rhino = true; % if set to false it uses spm coreg
    coreg_set.fid.label.nasion = 'Nasion';
    coreg_set.fid.label.lpa = 'LPA';
    coreg_set.fid.label.rpa = 'RPA';
    D = osl_headmodel(coreg_set);
    if run.preproc.check_coreg
      h=report.coreg(D);
    end
    D.save
  end
end


%% 8. BEAMFORM
if run.preproc.parcel
  for s = subs
    files=dir([PATH_DATA sprintf('%s_*%s*.mat', prefix, sub(s).id)]);
    S = [];
    S.D = [PATH_DATA files.name];
    p = parcellation('dk_full'); % alternative ('fmri_d100_parcellation_with_PCC_tighterMay15_v2_8mm')
    D = osl_inverse_model(S.D,p.template_coordinates, 'pca_order',55);
  end
end

%% 9. PARCELATION
% NOTE: now it takes the PC from each parcel. Instead, potentially use the
% vertex with maximum gamma increase (subject specific).
if run.preproc.beamform
  for s = subs
    files=dir([PATH_DATA sprintf('%s_*%s*.mat', prefix, sub(s).id)]);
    D = spm_eeg_load([PATH_DATA files.name]);
    D = D.montage('switch',2);
    D = ROInets.get_node_tcs(D,p.parcelflag,'spatialBasis','Giles');
    % 21 - Left postcentral S1, 23 - Left precentral M1
    D.save
  end
end

%% 10. CHECK COREGISTRATION
if run.preproc.check_coreg
  for s = subs
    D = spm_eeg_load([PATH_DATA files(s).name]);
    h = report.coreg(D);
    %         keyboard;
    %         close (h);
  end
end

%% 11. EPOCH
if run.preproc.epoch
  for s = subs
    files=dir([PATH_DATA sprintf('%s_*%s*.mat', prefix, sub(s).id)]);
    D_cont = spm_eeg_load([PATH_DATA, files.name]);
    S = [];
    S.timewin = [time_epoch(1) time_epoch(2)];
    S.trialdef(1).conditionlabel = 'Move_on';
    S.trialdef(1).eventtype = 'Move_on';
    S.trialdef(1).eventvalue = 66;
    S.reviewtrials = 0;
    S.save = 0;
    S.epochinfo.padding = 0;
    S.event = D_cont.events;
    S.fsample = D_cont.fsample;
    S.timeonset = D_cont.timeonset;
    S.D = D_cont;
    [epochinfo.trl, epochinfo.conditionlabels, S2] = spm_eeg_definetrial(S);
    
    % do epoching
    S2 = [];
    S2 = epochinfo;
    S2.D = D_cont;
    S2.prefix = 'e';
    D = osl_epoch(S2);
    if ~keep, delete(S2.D);  end
  end
end

