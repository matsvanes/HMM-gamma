function [spectra, tfmt, tfmt_avg, time_tf] = hmmspectramt_ft(X, T, time, Gamma_pad, options, keeptrials)

if ~exist('keeptrials', 'var') || isempty(keeptrials), keeptrials = 0; end
for k=1:numel(T)
  ntrl(k) = numel(T{k});
end
X = cat(1,X{:});
T = cat(1,T{:});

% combine data and gamma's
X = X .* Gamma_pad;

% reshape into trials and time and put in FT struct.
dat = [];
for k=1:options.K
  dat.label{k} = sprintf('state %d',k);
end
tmptrl = permute(reshape(X, T(1), [], options.K), [2,3,1]);
for k=1:size(tmptrl,1)
  dat.trial{k,1} = squeeze(tmptrl(k,:,2:end));
  dat.time{k,1} = time(2:end);
end

cfg=[];
cfg.method = 'mtmfft';
cfg.taper = 'dpss';
cfg.foi = options.fpass;
cfg.tapsmofrq = options.smo;
f = ft_freqanalysis(cfg, dat);

cfg.method = 'mtmconvol';
cfg.t_ftimwin = 0.2*ones(1, numel(cfg.foi));
cfg.toi = -1:0.05:1;
cfg.keeptrials = 'yes';
tf = ft_freqanalysis(cfg, dat);

% shape output as spectramt
spectra=[];
for k=1:options.K
  spectra.state(k).f = f.freq;
  spectra.state(k).psd = f.powspctrm(k,:)';
end

% average over trials per subject, than over subjects for group average
time_tf = tf.time;
f_tf = tf.freq;
for k=1:numel(ntrl)
  tfmt{k} = permute(squeeze(mean(tf.powspctrm(sum(ntrl(1:k))-ntrl(k)+1:sum(ntrl(1:k)),:,:,:),2)), [1,3,2]);
  tfmt_avg{k} = squeeze(mean(tfmt{k}));
  if ~keeptrials
    tfmt{k} = tfmt_avg{k};
  end
end
tfmt_avg = mean(cat(3,tfmt_avg{:}),3);


