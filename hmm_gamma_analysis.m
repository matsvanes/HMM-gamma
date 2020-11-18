%% JOBS
if ~exist('run', 'var') || isempty(run)
  run.preproc.bpfilt  = 0;
  run.preproc.whiten  = 1;
  run.remove_parc     = 0;
  run.TF.run          = 1;
  run.TF.eval         = 0;
  run.HMM.prep        = 0;
  run.HMM.run         = 0;
  run.HMM.eval        = 0;
  run.HMM.ST          = 0;
  if ~isfield(run.TF, 'ROI') || isempty(run.ROI)
    run.ROI          = input('which ROIs do you want to analyze?{sensor, parc, M1'); % can be 'sensor', 'parc', 'M1'
  end
  run.TF
  run.HMM
  disp(['%%% Check if correct jobs are defined! %%%']);keyboard
else
  run.TF
  run.HMM
  disp(['%%% The following jobs are defined! %%%']);
end


%% DEFINE
if exist('/ohba/pi/mwoolrich/', 'dir'), islocal=0; else, islocal=1; end
if islocal
  PATH.BASE = '/Volumes/T5_OHBA/'; % if on a local desktop
  PATH.ORIG = [];
  islocal   = 1;
  warning('PATH.ORIG not accessible')
else
  PATH.BASE = '/ohba/pi/mwoolrich/';
  PATH.ORIG = [PATH.BASE, 'datasets/CZ_Gamma/MEG_Gamma_HandOver/'];
  PATH.BASE = [PATH.BASE, 'mvanes/'];
  islocal   = 0;
end
PATH.ANALYSIS = [PATH.BASE, 'analysis/HMM-gamma/'];
PATH.DATA     = [PATH.ANALYSIS 'data/'];
PATH.ORIGDATA = [PATH.ORIG, 'data/'];
PATH.TF       = [PATH.ANALYSIS, 'TF/'];
PATH.HMM      = [PATH.ANALYSIS 'HMM/'];

PATH.SCRIPT = [PATH.BASE, 'scripts/HMM-gamma/'];
 if ~exist('userdir', 'var'), userdir = []; else, if ~strcmp(userdir(end),'/'), userdir = [userdir, '/']; end, end

files=dir([PATH.DATA 'efd_*.mat']);
subinfo;
prefix = 'fd';
if run.preproc.bpfilt,  prefix = ['f', prefix]; end
if run.preproc.whiten,  prefix = ['w', prefix]; end
prefix = ['e', prefix];
p = parcellation('dk_full');
%% PARAMS

% TF
keep = 0;
megtype = 'neuromag_sss';
if ~exist('time_bl', 'var'), time_bl = [-1 -0.5]; end
if ~exist('freq',    'var'), freq    = 5:1:145;   end
if ~exist('freqres', 'var'), freqres = 2.5;       end
if ~exist('timres',  'var'), timres  = 400;       end

sensor = 0;

% HMM
round_factor = 1000;
order        = 5;
N_states     = [3:1:10];
realization  = [1:10];

if ~exist('subs', 'var'), subs = 1:length(files); end

%% RUN TF
if run.TF.run
  
  
  for rois = 1:numel(run.ROI)
    
    PATH.TARGET = [PATH.TF, run.ROI{rois}, '/', userdir];
    PATH.DIPOLE   = [PATH.TARGET, 'peaks_60_90Hz/'];
    dipole_group = load([PATH.TARGET, prefix, '_dip_index_sorted.mat']);
    for s = subs
      s
      files=dir([PATH.DATA sprintf('%s_*%s*.mat', prefix, sub(s).id)]);
      dipole_subject = load([PATH.DIPOLE, sprintf('%s_case_%s', prefix, sub(s).id), '.mat']);
      D = hmm_gamma_preparedata(PATH, files.name, run.ROI{rois}, run.remove_parc, p, dipole_group, dipole_subject);
      
      % TF
      S = [];
      S.D = D;
      S.frequencies = freq;
      S.timewin = [-inf inf];
      S.phase = 0;
      S.method = 'mtmconvol';
      S.settings.taper = 'dpss';
      S.settings.timeres = timres;
      S.settings.timestep = 50;
      S.settings.freqres = freqres;
      D = spm_eeg_tf(S);
      if ~keep, delete(S.D);  end
      
      % Crop
      S = [];
      S.D = D;
      S.timewin = [-1800 1800];
      D = spm_eeg_crop(S);
      if ~keep, delete(S.D);  end
      
      % Average
      if ~strcmp(run.ROI{rois}, 'M1')
        S = [];
        S.D = D;
        S.robust.ks = 5;
        S.robust.bycondition = false;
        S.robust.savew = false;
        S.robust.removebad = 0;
        S.circularise = false;
        D = spm_eeg_average(S);
        if ~keep, delete(S.D);  end
        
        % Log
        S = [];
        S.D = D;
        S.method = 'LogR';
        S.timewin = [time_bl(1)*1000 time_bl(2)*1000];
        S.pooledbaseline = 1;
        D = spm_eeg_tf_rescale(S);
      end
    end
  end
end

%% EVAL TF
if  run.TF.eval
  plot_M1_TF
end

%% PREP HMM
if run.HMM.prep
  
  for rois = 1:numel(run.ROI)
    PATH.TARGET = [PATH.HMM, run.ROI{rois}, '/', userdir];
      dipole_group = load([PATH.TF, run.ROI{rois}, '/','optimised/', 'effd_dip_index_sorted.mat']);
    for s = subs
      s
      files=dir([PATH.DATA sprintf('%s_*%s*.mat', prefix, sub(s).id)]);
      [D, dat] = hmm_gamma_preparedata(PATH, files.name, run.ROI{rois}, run.remove_parc, p, dipole_group);
      
      % select parcel/voxel of interest.
      if strcmp(run.ROI{rois}, 'parc')
        if run.remove_parc, POI = 14; else, POI = 23; end
      elseif strcmp(run.ROI{rois},'M1')
        dipole_subject = load([PATH.TF, run.ROI{rois}, '/','optimised/peaks_60_90Hz/', sprintf('effd_case_%s_sel.mat', sub(s).id)]);
        POI = dipole_subject.voxel;
      end
      
      PAC{s} = squeeze(dat(POI,:,:));
      t{s} = D.time(1):1/D.fsample:D.time(end);
      t{s} = round(t{s}*round_factor)/round_factor; % ensure that it's rounded properly for 'find' function
      
      % PC_1st: 1xlength(files) cell, whereby each cell contains a samples(currently 301)xtrials matrix
      nsamples{s} = size(PAC{s},1);
      ntrials{s} = size(PAC{s},2);
      
      %format for HMM
      X{s}(:,1) = PAC{s}(:); % change format. Within each subject/file cell from matrix (smaplesx trials) to vector (concatenate all trials)
      T{s} = repmat(nsamples{s},1,ntrials{s})';
    end %subj
    if run.remove_parc
      fname = 'PREP_HMM_sel.mat';
    else
      fname = 'PREP_HMM.mat';
    end
    save([PATH.HMM fname],'X','T','ntrials','t','round_factor','order','D');
    
    % MAR check
%     X_all=cell2mat(X');
%     mkdir([PATH.HMM 'order_check/'])
%     cd([PATH.SCRIPT 'spectra_check/']);
%     for oo=1:20
%       hmmmar_spectra_check(X_all, oo, D.fsample);
%       print('-dpng',[PATH.HMM,'order_check/order_',sprintf('%03d.png',oo)]);
%       close(gcf)
%     end
%     hmmmar_spectra_check(X_all, order, D.fsample);
  end %if prep HMM
end

%% RUN HMMs
if run.HMM.run
  for n=1:length(N_states)
    for r=1:length(realization)
      PATH.HMM_PREC = [PATH.HMM 'order_' num2str(order) '/' num2str(N_states(n)) '_states/real_' num2str(realization(r)) '/'];
      mkdir(PATH.HMM_PREC);
      disp(['%%% run HMM: order ',num2str(order),', ',num2str(N_states(n)),' states, real ',num2str(realization(r)),' %%%']);
      
      load([PATH.HMM 'PREP_HMM.mat']);
      
      % define HMM params
      load([PATH.SCRIPT 'meg_TMS_hmm_options.mat']);
      options.K = N_states(n);
      options.Fs = D.fsample;
      options.order = order;
      
      % run HMM
      clear hmm Gamma
      [hmm{1},Gamma] = hmmmar(X,T,options);
      
      % get time vector
      for s=1:length(T)
        t{s}=D.time(1):1/D.fsample:D.time(end);
        for o=1:ntrials{s}
          t_Gamma{s}{o}=t{s}(hmm{1}.train.order+1):1/D.fsample:t{s}(end);
        end
      end
      
      % estimate gammas & sepctrum
      [Gamma_mat_Gavg_inf,Gamma_mat_avg_inf,Gamma_mat_inf,Gamma_t] = hmm_get_gammas(X,T,t_Gamma,hmm,options,[],ntrials);
      [spectra_t,options_mar,options_mt] = hmm_get_spectra(X,T,D.fsample,Gamma_t,hmm,1,options,256);
      
      % save Vars
      save([PATH.HMM_PREC 'POST_HMM'],'hmm','X','T','ntrials','t','round_factor','order','D','t_Gamma','Gamma*','spectra*','options*','-v7.3');
    end %real
  end %states
end

