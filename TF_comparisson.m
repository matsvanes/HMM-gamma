% Compare the TF from different analyses: sensor level, M1 dipole, and HMM.
% (can be extended with parcel level)
subinfo

% general parameters
blwindow = [-.6 -.2];
tsel=[-0.5 0.5];
prefix = 'efd'; % 'effd';
cmap = inferno(64);
%% sensor
fname = sprintf('/Volumes/T5_OHBA/analysis/HMM-gamma/TF/sensor/%s_TF_sensor_groupAvg.mat', prefix);
load(fname)

chans = {'MEG1811', 'MEG0441', 'MEG0411', 'MEG1821', 'MEG0431', 'MEG0421', 'MEG0631', 'MEG0741', 'MEG0711', 'MEG1812', 'MEG1822', 'MEG0432', 'MEG0742', 'MEG0712', 'MEG1813', 'MEG1823', 'MEG0433', 'MEG0423', 'MEG0633', 'MEG0743', 'MEG0713'};
idx = find(contains(label, chans));

t1 = nearest(time, tsel(1));
t2 = nearest(time, tsel(2));
t_short = time(t1:t2);
tf_short = squeeze(mean(tf(:,idx,:,t1:t2),2));

% subjects
for s=1:numel(sub)
  figure(s)
  subplot(1,3,1)%subplot(2,2,1)
  plot_TF(t_short, freq, squeeze(tf_short(s,:,:)), 'Fourier sensor', [], cmap)
end

% Group
figure(34);
subplot(1,3,1)%subplot(2,2,1)
plot_TF(t_short, freq, squeeze(mean(tf_short,1)),'Fourier sensor', [], cmap)

%% M1
d = '/Volumes/T5_OHBA/analysis/HMM-gamma/TF/M1/peaks_60_90Hz/';
tfm1 = zeros(33,31,21);
t_short = -.5:0.05:.5;
freq=60:90;
for s=1:numel(sub)
  tmp = load([d, sprintf('%s_case_%s.mat', prefix, sub(s).id)]);
  tfm1(s,:,:) = tmp.tf_2d(:,t1:t2);
end

% subjects
for s=1:numel(sub)
  figure(s)
  subplot(1,3,2)%subplot(2,2,2)
  plot_TF(t_short, freq, squeeze(tfm1(s,:,:)), 'Fourier M1 dipole', [], cmap)
end

% Group
figure(34);
subplot(1,3,2)% subplot(2,2,2)
plot_TF(t_short, freq, squeeze(mean(tfm1,1)),'Fourier M1 dipole', [], cmap)


%% HMM
d = load(sprintf('/Volumes/T5_OHBA/analysis/HMM-gamma/HMM/TDE/4_states/lag_3/%s_POST_HMM.mat', prefix), 't', 'tf','tfmt','tfmt_avg', 'spectramt');

f=d.spectramt.state(1).f;
t=d.t{1};
tf=d.tfmt;

f1=nearest(f,60);
f2=nearest(f,90);

t1 = nearest(t, tsel(1));
t2 = nearest(t, tsel(2));

tfhmm = permute(cat(3,tf{:}),[3 2 1]);


% subjects
for s=1:numel(sub)
  figure(s)
  subplot(1,3,3)%   subplot(2,2,3+k-1)
  plot_TF(t(t1:t2), f(f1:f2), squeeze(tfhmm(s,f1:f2,:)), 'HMM M1 dipole', [], cmap)
end

% Group
figure(34);
subplot(1,3,3) %subplot(2,2,3+k-1)
plot_TF(t(t1:t2), f(f1:f2), squeeze(mean(tfhmm(:,f1:f2,:),1)),'HMM M1 dipole', [], cmap)


%% Save
savedir = '/Volumes/T5_OHBA/analysis/HMM-gamma/TF/';
for s=1:numel(sub)
  figure(s)
  suptitle(sprintf('TF - %s', sub(s).id))
  print('-dpng',[savedir, sprintf('%s_TF_comparisson_%s', prefix, sub(s).id)]);pause(1);
end

figure(34);
suptitle('TF - Group average')
print('-dpng',[savedir, sprintf('%s_TF_comparisson_groupAvg', prefix)]);pause(1);


%% Sensor level - no preprocessing, BP filter, prewhitening
prefixes = {'efd', 'effd', 'ewfd'};
for i = 1:numel(prefixes)
  
  fname = sprintf('/Volumes/T5_OHBA/analysis/HMM-gamma/TF/sensor/%s_TF_sensor_groupAvg.mat', prefixes{i});
  load(fname)
  
  chans = {'MEG1811', 'MEG0441', 'MEG0411', 'MEG1821', 'MEG0431', 'MEG0421', 'MEG0631', 'MEG0741', 'MEG0711', 'MEG1812', 'MEG1822', 'MEG0432', 'MEG0742', 'MEG0712', 'MEG1813', 'MEG1823', 'MEG0433', 'MEG0423', 'MEG0633', 'MEG0743', 'MEG0713'};
  idx = find(contains(label, chans));
  
  t1 = nearest(time, tsel(1));
  t2 = nearest(time, tsel(2));
  t_short = time(t1:t2);
  tf_short = squeeze(mean(tf(:,idx,:,t1:t2),2));
  
  % subjects
  for s=1:numel(sub)
    figure(s)
    subplot(1,3,i)
    plot_TF(t_short, freq, squeeze(tf_short(s,:,:)), 'Fourier sensor')
  end
  
  % Group
  figure(34);
  subplot(1,3,i)
  plot_TF(t_short, freq, squeeze(mean(tf_short,1)),'Fourier sensor')
end






