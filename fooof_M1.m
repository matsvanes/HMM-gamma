%% settings

PATH  = ['/Volumes/T5_OHBA/analysis/HMM-gamma/TF/M1/'];
cmap = inferno(64);
if ~exist('subs', 'var'),        subs        = 1:33;                         end
if ~exist('doplot', 'var'),      doplot      = 1;                            end
if ~exist('N_states', 'var'),    N_states    = 4;                            end
if ~exist('lag', 'var'),         lag         = 3;                            end
if ~exist('prefix', 'var'),      prefix      = input('what is the prefix?'); end
subinfo;
% FOOOF settings
f_range = [60 90];
settings = struct();
settings.max_n_peaks=2;

%% Prepare Traditional TF
tf_group{1} = nan(33, 31, 41);
for s=subs
  d = dir([PATH, sprintf('ptf_%s*%s*.mat', prefix, sub(s).id)]);
  s
  D = spm_eeg_load([PATH, d.name]);
  if s==1
    t{1} = D.time;
    t0 = nearest(t{1}, -1);
    t1 = nearest(t{1}, 1);
    t{1} = t{1}(t0:t1);
    f{1} = D.frequencies;
    f0 = nearest(f{1},60);
    f1 = nearest(f{1},90);
    f{1} = f{1}(f0:f1);
  end
  tf = D(:,f0:f1,t0:t1,:);
  chanlabels = D.chanlabels;
  
  % max voxel in spectrospatial space
  tmp = extractBetween(d.name, sprintf('%s_', prefix), '_go');
  filename = [PATH,'peaks_60_90Hz/', prefix,'_', tmp{1}];
  load([filename, '.mat'],'voxel');
  
  tf_group{1}(s,:,:) = squeeze(nanmean(tf(voxel, :,:,:),4));
end

%% Prepare HMM TF
load(sprintf('/Volumes/T5_OHBA/analysis/HMM-gamma/HMM/TDE/%d_states/lag_%d/%s_POST_HMM.mat', N_states, lag, prefix), 'spectramt', 'tfmt', 'tfmt_avg')

t{2} = -1:1/500:1;
f{2}= spectramt.state(1).f;
f0 = nearest(f{2},60);
f1 = nearest(f{2},90);
f{2} = f{2}(f0:f1);

% make sure the input is as expected
tf_group{2} = permute(cat(3, tfmt{:}), [3,2,1]);
tf_group{2} = tf_group{2}(:,f0:f1,:);
if sign(min(tf_group{2}(:)))==-1
  tf_group{2} = tf_group{2} - 2*min(tf_group{2}(:)); % make sure all values are >0 --> probably needed because frequencies are centered
end
tfmt_avg = tfmt_avg - 2*min(tfmt_avg(:));


%% FOOOF
% fit on individual subjects
for ii=1:numel(tf_group)
  [s1,s2,s3] = size(tf_group{ii});
  power_spectrum{ii} = zeros(s1,s2,s3);
  fooofed_spectrum{ii} = zeros(s1,s2,s3);
  ap_fit{ii} = zeros(s1,s2,s3);
  tf_min_apavg{ii} = zeros(s1,s2,s3);
  t1(ii) = nearest(t{ii},-0.5);
  t2(ii) = nearest(t{ii},0);
  t3(ii) = nearest(t{ii},1);
  for s=subs
    % fit on time-average over baseline and movement
    psd_off = mean(tf_group{ii}(s,:,t1(ii):t2(ii)),3);
    psd_on = mean(tf_group{ii}(s,:,t2(ii):t3(ii)),3);
    ff_on = fooof(f{ii}, psd_on, f_range, settings, true);
    ff_off = fooof(f{ii}, psd_off, f_range, settings, true);
    
    for k=1:numel(t{ii})
      % fit on every time point
      psd = tf_group{ii}(s,:,k);
      x = fooof(f{ii}, psd, f_range, settings, true);
      power_spectrum{ii}(s,:,k) = x.power_spectrum;
      fooofed_spectrum{ii}(s,:,k) = x.fooofed_spectrum;
      ap_fit{ii}(s,:,k) = x.ap_fit;
      
      % subtract temporal average ap-fit from original tf
      if k<t2
        tf_min_apavg{ii}(s,:,k) = log10(tf_group{ii}(s,:,k))-ff_off.ap_fit;
      else
        tf_min_apavg{ii}(s,:,k) = log10(tf_group{ii}(s,:,k))-ff_on.ap_fit;
      end
    end
  end
  tf_min_ap{ii} = power_spectrum{ii}-ap_fit{ii}; % subtract ap-fit from original tf
  ff_min_ap{ii} = fooofed_spectrum{ii}-ap_fit{ii}; % subtract ap-fit from modelled tf
  
  % Average over subjects
  tf_min_apavg_group{ii} = squeeze(mean(tf_min_apavg{ii},1));
  ff_min_ap_group{ii} = squeeze(mean(ff_min_ap{ii},1));
  tf_min_ap_group{ii} = squeeze(mean(tf_min_ap{ii},1));
end

%%
if doplot
  str = {'traditional TF','HMM TF'};
  figure(34);
  for ii=1:numel(tf_group)
    subplot(2,3,(ii-1)*3+1)
    % relative change from baseline
    X = 10.^ff_min_ap_group{ii};
    X=X./mean(X(:,t1(ii):t2(ii)),2)-1;
    plot_TF(t{ii},f{ii}, X, ...
      {sprintf('%s', str{ii}), ' individual fit per timepoint', 'fooofed spectrum minus aperiodic fit'},...
      [], cmap)
    subplot(2,3,(ii-1)*3+2)
    X = 10.^tf_min_ap_group{ii};
    X=X./mean(X(:,t1(ii):t2(ii)),2)-1;
    plot_TF(t{ii},f{ii}, X, ...
      {sprintf('%s', str{ii}), 'individual fit per timepoint', 'power spectrum minus aperiodic fit'},...
      [], cmap)
    subplot(2,3,(ii-1)*3+3)
    X = 10.^tf_min_apavg_group{ii};
    X=X./mean(X(:,t1(ii):t2(ii)),2)-1;
    plot_TF(t{ii},f{ii}, X, ...
      {sprintf('%s', str{ii}), 'fit on pre and post', 'power spectrum minus aperiodic fit'},...
      [], cmap)
    suptitle('FOOOF group')
  end
  
  for s=subs
    figure(s);
    for ii=1:numel(tf_group)
      subplot(2,3,(ii-1)*3+1)
      X = 10.^squeeze(ff_min_ap{ii}(s,:,:));
      X=X./mean(X(:,t1(ii):t2(ii)),2)-1;
      plot_TF(t{ii},f{ii}, X, ...
        {sprintf('%s', str{ii}), ' individual fit per timepoint', 'fooofed spectrum minus aperiodic fit'},...
        [], cmap)
      subplot(2,3,(ii-1)*3+2)
      X = 10.^squeeze(tf_min_ap{ii}(s,:,:));
      X=X./mean(X(:,t1(ii):t2(ii)),2)-1;
      plot_TF(t{ii},f{ii}, X, ...
        {sprintf('%s', str{ii}), 'individual fit per timepoint', 'power spectrum minus aperiodic fit'},...
        [], cmap)
      subplot(2,3,(ii-1)*3+3)
      X = 10.^squeeze(tf_min_apavg{ii}(s,:,:));
      X=X./mean(X(:,t1(ii):t2(ii)),2)-1;
      plot_TF(t{ii},f{ii}, X, ...
        {sprintf('%s', str{ii}), 'fit on pre and post', 'power spectrum minus aperiodic fit'},...
        [], cmap)
    end
    suptitle(sprintf('FOOOF subject %d',s), 'Interpreter', 'none')
  end
end







