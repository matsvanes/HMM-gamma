function tfconvol = hmmtimefreqconvol(cfg, Gamma, T, time, F0)
% obtains a time-frequency representation of state time courses, using
% convolution (in contrast to hmmtimefreq, which uses the HMM spectra). 
% The argument Gamma contains the state time courses. 
%
% The output arguments are
%  psd_tf: (time points by no. of frequency bins by no.regions) 
%
% Author: Mats van Es, OHBA, University of Oxford (2021)


if ~isfield(cfg, 'baseline'), cfg.baseline = [-.6 -.2]; end
if ~isfield(cfg, 'foi'), cfg.foi    = 60:90;     end 
if ~isfield(cfg, 'tapsmofrq'), cfg.tapsmofrq = 5;         end 
if ~isfield(cfg, 't_ftimwin'), cfg.t_ftimwin  = 0.200*ones(1,length(cfg.foi));       end 
cfg.method = 'mtmconvol';
cfg.output = 'pow';
cfg.taper = 'dpss';
cfg.toi = time(1):0.05:time(end);
cfg.pad = 2^nextpow2(time(end)-time(1)+diff(time(1:2)));
cfg.keeptrials = 'yes';

remember = [osldir, '/osl-core/osl_startup'];
osl_shutdown()
pathinfo; addpath(fieldtrippath)
ft_defaults
nstates = size(Gamma,3);

dat=[];
dat.trial = permute(Gamma, [2,3,1]);
for k=1:nstates
  dat.label{k} = sprintf('chan%02d', k);
end
dat.time = time;

freq = ft_freqanalysis(cfg, dat);
restoredefaultpath;
run(startupdir)
run(remember);

for k=1:numel(T)
  ntrl(k) = numel(T{k});
end

q = 1;
while numel(downsample(F0(:,1),q)) ~= numel(freq.time)
  q=q+1;
end
F0ds = downsample(F0,q);

idx(1) = nearest(dat.time, cfg.baseline(1));
idx(2) = nearest(dat.time, cfg.baseline(2));
for k=1:numel(ntrl)
  psd_tf{k} = permute(squeeze(nanmean(freq.powspctrm(sum(ntrl(1:k))-ntrl(k)+1:sum(ntrl(1:k)),:,:,:))), [3,2,1]);
  [s1,s2,~] = size(psd_tf{k});
%   psd_tf{k} = psd_tf{k}-repmat(nanmean(psd_tf{k},1),[s1,1,1]);
  psd_tf_norm{k} = mean(psd_tf{k}.*permute(repmat(F0ds, [1,1,s2]), [1,3,2]),3);
end

tfconvol = [];
tfconvol.tf = psd_tf;
tfconvol.tf_norm = psd_tf_norm;
tfconvol.time = freq.time;
tfconvol.freq = freq.freq;
