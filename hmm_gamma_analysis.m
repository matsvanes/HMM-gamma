%% JOBS
if ~exist('run', 'var') || isempty(run)
  run.remove_parc     = 1;
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
  PATH_BASE = '/Volumes/T5_OHBA/'; % if on a local desktop
  warning('PATH_ORIG not accessible')
  PATH_ORIG = [];
  islocal = 1;
else
  PATH_BASE = '/ohba/pi/mwoolrich/';
  PATH_ORIG = [PATH_BASE, 'datasets/CZ_Gamma/MEG_Gamma_HandOver/'];
  PATH_BASE = [PATH_BASE, 'mvanes/'];
  islocal = 0;
end
PATH_ANALYSIS = [PATH_BASE, 'analysis/HMM-gamma/'];
PATH_DATA =  [PATH_ANALYSIS 'data/'];
PATH_ORIGDATA = [PATH_ORIG, 'data/'];
PATH_TF = [PATH_ANALYSIS, 'TF/'];
PATH_HMM = [PATH_ANALYSIS 'HMM/'];

PATH_SCRIPT = [PATH_BASE, 'scripts/HMM-gamma/'];

files=dir([PATH_DATA 'efd_*.mat']);

%% PARAMS

time_bl     = [-1 -0.5];


% TF
keep = 0;
megtype = 'neuromag_sss';
freq = 5:1:145;
res = 2.5;
sensor = 0;

% HMM
round_factor = 1000;
order = 5;
N_states = [3:1:10];
realization = [1:10];

if ~exist('subs', 'var'), subs = 1:length(files); end

%% RUN TF
if run.TF.run
  files=dir([PATH_DATA 'efd_*.mat']);
  
  for rois = 1:numel(run.ROI)
    tmpPATH_TF = [PATH_TF, run.ROI{rois}, '/'];
    
    TF_avg = [];
    for s = subs
      [D, POI] = hmm_gamma_preparedata(PATH_DATA, tmpPATH_TF, files(s).name, run.ROI{rois}, run.remove_parc);
            
      % TF
      S = [];
      S.D = D;
      S.frequencies = freq;
      S.timewin = [-Inf Inf];
      S.phase = 0;
      S.method = 'mtmconvol';
      S.settings.taper = 'dpss';
      S.settings.timeres = 400;
      S.settings.timestep = 50;
      S.settings.freqres = res;
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
      end
      
      % Log
      S = [];
      S.D = D;
      if strcmp(run.ROI{rois}, 'M1')
        S.method = 'Diff';
      else
        S.method = 'LogR';
      end
      S.timewin = [time_bl(1)*1000 time_bl(2)*1000];
      S.pooledbaseline = 1;
      D = spm_eeg_tf_rescale(S);
      if strcmp(run.ROI{rois}, 'parc')
        if run.remove_parc
          TF_avg(s,:,:)=squeeze(D(14,:,:,:)); %POI = 14; % Precentral after removal of some ROIs, otherwise POI=23;
        else
          TF_avg(s,:,:)=squeeze(D(23,:,:,:));
        end
      elseif strcmp(run.ROI{rois}, 'sensor')
        TF_avg(s,:,:,:)=D(:,:,:);
      end
    end
    
    if strcmp(run.ROI{rois}, 'parc') || strcmp(run.ROI{rois}, 'sensor')
      save([tmpPATH_TF 'TF_avg.mat'],'TF_avg','files');
    end
  end
end

%% EVAL TF
if  run.TF.eval
  load ([PATH_TF 'TF_avg.mat']);
  % visulaisation (GA,SS)
  meg_tms_tf_vis(TF_avg,D,freq,time_bl,files,PATH_ORIGDATA,PATH_TF)
  keyboard;
  % peak frequency
  [freq_peak,power]=meg_tms_tf_peak(TF_avg,D,freq,time_pre,time_peri,time_post,0,1,1);
  % eval and corr
  meg_tms_tf_eval(PATH_BASE,PATH_TF,files,freq_peak,power)
end

%% PREP HMM
if run.HMM.prep
files=dir([PATH_DATA 'efd_*.mat']);
for rois = 1:numel(run.ROI)
    tmpPATH_HMM = [PATH_HMM, run.ROI{rois}, '/'];
  for s = subs
    [D, POI, dat] = hmm_gamma_preparedata(PATH_DATA, tmpPATH_TF, files(s).name, run.ROI{rois}, run.remove_parc);

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
  save([PATH_HMM 'PREP_HMM.mat'],'X','T','ntrials','t','round_factor','order','D');
  
  % MAR check
  X_all=cell2mat(X');
  mkdir([PATH_HMM 'order_check/'])
  cd([PATH_SCRIPT 'spectra_check/']);
  for oo=1:20
    hmmmar_spectra_check(X_all, oo, D.fsample);
    print('-dpng',[PATH_HMM,'order_check/order_',sprintf('%03d.png',oo)]);
    close(gcf)
  end
  hmmmar_spectra_check(X_all, order, D.fsample);
end %if prep HMM
end

%% RUN HMMs
if run.HMM.run
  for n=1:length(N_states)
    for r=1:length(realization)
      PATH_HMM_PREC = [PATH_HMM 'order_' num2str(order) '/' num2str(N_states(n)) '_states/real_' num2str(realization(r)) '/'];
      mkdir(PATH_HMM_PREC);
      disp(['%%% run HMM: order ',num2str(order),', ',num2str(N_states(n)),' states, real ',num2str(realization(r)),' %%%']);
      
      load([PATH_HMM 'PREP_HMM.mat']);
      
      % define HMM params
      load([PATH_SCRIPT 'meg_TMS_hmm_options.mat']);
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
      save([PATH_HMM_PREC 'POST_HMM'],'hmm','X','T','ntrials','t','round_factor','order','D','t_Gamma','Gamma*','spectra*','options*','-v7.3');
    end %real
  end %states
end

