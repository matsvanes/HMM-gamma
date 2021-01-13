% combine sensor level TF from all subjects.


PATH_BASE = '/Volumes/T5_OHBA/analysis/HMM-gamma/';
PATH  = [PATH_BASE 'TF/sensor/'];
subinfo;
subs=1:33;
prefix = 'efd'; % 'effd' or 'ewfd'

tf = zeros(33,306,31,72);
for s=subs
  files = dir([PATH, sprintf('rmptf_%s*%s*.mat', prefix, sub(s).id)]);
  D = spm_eeg_load([PATH, files.name]);
  chidx = find(contains(D.chantype, 'MEG'));
  tf(s,:,:,:) = D(chidx,:,1:72);
  if s==1
    time = D.time;
    label = D.chanlabels;
    freq = D.frequencies;
  end 
end

save([files.folder, sprintf('/%s_TF_sensor_groupAvg.mat', prefix)], 'tf', 'time', 'label', 'freq')