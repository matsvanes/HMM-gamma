%% JOBS
if ~exist('todo', 'var') || isempty(todo)
  todo.preproc.bpfilt  = 0;
  todo.preproc.whiten  = 1;
  todo.remove_parc     = 0;
  todo.TF.run          = 0;
  todo.TF.eval         = 0;
  todo.HMM.prep        = 1;
  todo.HMM.run         = 1;
  todo.HMM.model       = 'tde'; % can be 'mar', 'tde'
  if ~isfield(todo.TF, 'ROI') || isempty(todo.ROI)
    todo.ROI          = {'M1'};%input('which ROIs do you want to analyze?{sensor, parc, M1'); % can be 'sensor', 'parc', 'M1'
  end
  todo.TF
  todo.HMM
  disp(['%%% Check if correct jobs are defined! %%%']);keyboard
else
  todo.TF
  todo.HMM
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
if ~exist('userdir', 'var') || isempty(userdir), userdir = []; else, if ~strcmp(userdir(end),'/'), userdir = [userdir, '/']; end, end

files=dir([PATH.DATA 'efd_*.mat']);
subinfo;
prefix = 'fd';
if todo.preproc.bpfilt,  prefix = ['f', prefix]; end
if todo.preproc.whiten,  prefix = ['w', prefix]; end
prefix = ['e', prefix];
p = parcellation('dk_full');
%% PARAMS

% TF
keep = 0;
megtype = 'neuromag_sss';
if ~exist('time_bl', 'var'), time_bl = [-.6 -.2]; end %  [-1 -0.5];
if ~exist('freq',    'var'), freq    = 60:90;     end % 5:1:145
if ~exist('freqres', 'var'), freqres = 5;         end % 2.5
if ~exist('timres',  'var'), timres  = 200;       end % 400

sensor = 0;
if ~exist('timewin',  'var'), timewin  = [];       end

% HMM
if todo.HMM.run==1
  % general settings
  if ~exist('N_states', 'var'),    N_states      = [3:1:10]; end
  if ~exist('realization', 'var'), realization   = [1:10];   end
  if ~exist('round_factor', 'var'), round_factor = 1000; end
  if ~exist('rungroup', 'var'),     rungroup     = 1; end
  % model specific settings
  if strcmp(todo.HMM.model, 'mar')
    PATH.HMM = [PATH.HMM 'MAR/'];
    if ~exist('order', 'var'),        order        = 5;    end
  elseif strcmp(todo.HMM.model, 'tde')
    PATH.HMM = [PATH.HMM 'TDE/'];
    order = 0;
    covtype = 'full';
    zeromean = 1;
    if ~exist('embeddedlags', 'var'), embeddedlags = -3:3; end
    if ~exist('pca', 'var'),          pca = 0;             end
  end
end

if ~exist('subs', 'var'), subs = 1:length(sub); end

%% RUN TF
if todo.TF.run
  
  
  for rois = 1:numel(todo.ROI)
    
    PATH.TARGET = [PATH.TF, todo.ROI{rois}, '/', userdir];
    if strcmp(todo.ROI{rois}, 'M1') && todo.remove_parc
      PATH.DIPOLE   = [PATH.TARGET, 'peaks_60_90Hz/'];
      dipole_group = load([PATH.TARGET, prefix, '_dip_index_sorted.mat']);
    else
      dipole_group = [];
    end
    for s = subs
      s
      files=dir([PATH.DATA sprintf('%s_*%s*.mat', prefix, sub(s).id)]);
      if strcmp(todo.ROI{rois}, 'M1') && todo.remove_parc
        dipole_subject = load([PATH.DIPOLE, sprintf('%s_case_%s', prefix, sub(s).id), '.mat']);
      end
      D = hmm_gamma_preparedata(PATH, files.name, todo.ROI{rois}, todo.remove_parc, p, dipole_group);
      
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
      if ~strcmp(todo.ROI{rois}, 'M1')
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
if  todo.TF.eval
  plot_M1_TF
end

%% PREP HMM
if todo.HMM.prep
  
  for rois = 1:numel(todo.ROI)
    PATH.TARGET = [PATH.HMM, todo.ROI{rois}, '/', userdir];
    dipole_group = load([PATH.TF, todo.ROI{rois}, '/', 'effd_dip_index_sorted.mat']);
    for s = subs
      s
      files=dir([PATH.DATA sprintf('%s_*%s*.mat', prefix, sub(s).id)]);
      [D, dat] = hmm_gamma_preparedata(PATH, files.name, todo.ROI{rois}, todo.remove_parc, p, dipole_group, timewin);
      
      % select parcel/voxel of interest.
      if strcmp(todo.ROI{rois}, 'parc')
        if todo.remove_parc, POI = 14; else, POI = 23; end
      elseif strcmp(todo.ROI{rois},'M1')
        if isfield(todo, 'orig') && todo.orig==1
          dipole_subject = load([PATH.TF, todo.ROI{rois}, '/','peaks_60_90Hz_orig/', sprintf('case_%s.mat', sub(s).id)]);
        else
          if todo.remove_parc
            dipole_subject = load([PATH.TF, todo.ROI{rois}, '/','peaks_60_90Hz/', sprintf('effd_case_%s_sel.mat', sub(s).id)]);
          else
            dipole_subject = load([PATH.TF, todo.ROI{rois}, '/','peaks_60_90Hz/', sprintf('effd_case_%s.mat', sub(s).id)]);
          end
        end
        POI = dipole_subject.voxel;
      end
      
      PAC{s} = normalize(squeeze(dat(POI,:,:)),1);
      t{s} = D.time(1):1/D.fsample:D.time(end);
      t{s} = round(t{s}*round_factor)/round_factor; % ensure that it's rounded properly for 'find' function
      
      % PC_1st: 1xlength(files) cell, whereby each cell contains a samples(currently 301)xtrials matrix
      nsamples{s} = size(PAC{s},1);
      ntrials{s} = size(PAC{s},2);
      
      %format for HMM
      X{s}(:,1) = PAC{s}(:); % change format. Within each subject/file cell from matrix (smaplesx trials) to vector (concatenate all trials)
      T{s} = repmat(nsamples{s},1,ntrials{s})';
    end %subj
    fname = [prefix '_PREP_HMM'];
    if todo.remove_parc,   fname = [fname, '_sel'];    end
    if isfield(todo, 'orig') && todo.orig==1,    fname = [fname ,'_orig'];    end
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
if todo.HMM.run
  if rungroup==1, runsubs = 1; else, runsubs = subs; end
  for s=runsubs
    for n=1:length(N_states)
      %     for r=1:length(realization)
      %       PATH.HMM_PREC = [PATH.HMM 'order_' num2str(order) '/' num2str(N_states(n)) '_states/real_' num2str(realization(r)) '/'];
      if strcmp(todo.HMM.model, 'mar')
        PATH.HMM_PREC = [PATH.HMM 'order_' num2str(order) '/' num2str(N_states(n)) '_states/'];
      elseif strcmp(todo.HMM.model, 'tde')
        PATH.HMM_PREC = [PATH.HMM num2str(N_states(n)) '_states/'  'lag_' num2str(embeddedlags(end)) '/' ];
      end
      mkdir(PATH.HMM_PREC);
      disp(['%%% run HMM: order ',num2str(order),', ',num2str(N_states(n)),' states', ' %%%']);
      
      fname = [PATH.HMM, prefix, '_PREP_HMM'];
      if todo.remove_parc,  fname = [fname, '_sel']; end
      if isfield(todo, 'orig') && todo.orig==1,   fname = [fname, '_orig'];   end
      load(fname)
      
      if ~rungroup, X = X(s); T = T(s); ntrials = ntrials(s); t = t(s); end
      
      % define HMM params
      options = [];
      options.K = N_states(n);
      options.Fs = D.fsample;
      options.order = order;
      options.repetitions = numel(realization);
      options.useParallel = false;
      if strcmp(todo.HMM.model, 'tde')
        options.covtype = covtype;
        options.zeromean = zeromean;
        options.embeddedlags = embeddedlags;
        options.pca = pca;
      end
      
      % run HMM
      clear hmm Gamma
      [hmm,Gamma] = hmmmar(X,T,options);
      
      % get time vector
      for j=1:length(T)
        t{j}=D.time(1):1/D.fsample:D.time(end);
        for o=1:ntrials{j}
          t_Gamma{j}{o}=t{j}(hmm.train.order+1):1/D.fsample:t{j}(end);
        end
      end
      
      % estimate gammas & sepctrum
      %       [Gamma_mat_Gavg_inf,Gamma_mat_avg_inf,Gamma_mat_inf,Gamma_t] = hmm_get_gammas(X,T,t_Gamma,hmm,options,[],ntrials);
      %       [spectra_t,options_mar,options_mt] = hmm_get_spectra(X,T,D.fsample,Gamma_t,hmm,1,options,256);
      filename = [PATH.HMM_PREC, prefix, '_POST_HMM'];
      if ~rungroup, filename = [filename '_', sub(s).id]; end
      if todo.remove_parc,     filename = [filename, '_sel'];    end
      if isfield(todo, 'orig') && todo.orig==1,   filename = [fname, '_orig']; end
      
      save(filename,'hmm','X','T','ntrials','t','round_factor','order','D','t_Gamma','Gamma','options', '-v7.3');
      try
        [Gamma, MLGamma, dynamics, spectramt, tfmt, tfmt_avg] = hmm_gamma_hmm_post(X, Gamma, hmm, T, t, options, 1, filename);
        save(filename,'Gamma', 'MLGamma', 'dynamics', 'spectramt', 'tfmt', 'tfmt_avg', '-append');
      catch
        warning('hmm_gamma_hmm_post failed')
      end
    end %states
  end
end

