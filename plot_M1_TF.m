datadir = '/Volumes/T5_OHBA/analysis/HMM-gamma/TF/M1/';
datarank = 55;
% do group statistics on the M1 data

filename = 'group_TF_fieldtrip.mat';
load([datadir, filename]);
p = parcellation('dk_full');

cfg.frequency = [60 90];
cfg.latency = [0 0.5];
dat = ft_selectdata(cfg, dat);
% this is already versus baseline, so we can compare with zero
null = dat;
null.powspctrm = null.powspctrm*0;

cfg=[];
cfg.design = [1:33 1:33; ones(1,33) 2*ones(1,33)];
cfg.ivar = 2;
cfg.uvar = 1;
cfg.method = 'analytic';
cfg.statistic = 'depsamplesT';
stat = ft_freqstatistics(cfg, dat, null);
s=stat.stat;
s=mean(mean(s,3),2);

[s_sort, idx_sort] = sort(s,1, 'descend');
M1_idx = idx_sort(1:datarank);

% plot the average T-value for each location
steps = 6;
cmap = flipud(brewermap(2*steps, 'RdBu'));
cmap = cmap(steps+1:end,:);
xx = linspace(0, max(s(:)), steps);
for j=1:92
  B = xx-s(j);
  B(B<0)=inf;
  [~, ixtmp] = min(B);
  ix(j) = ixtmp;
end
for j=1:92
  c(j,:) = cmap(ix(j),:);
end
figure; hold on
sc=scatter3(p.template_coordinates(:,1),p.template_coordinates(:,2),p.template_coordinates(:,3),8,'MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5]);
sc.MarkerFaceAlpha = 0.1;sc.MarkerEdgeAlpha = 0.1;
for j=1:92
  if any(M1_idx==j)
    MarkerSize = 30;
  else
    MarkerSize = 15;
  end
  scatter3(p.roi_coordinates{23}(j,1),p.roi_coordinates{23}(j,2),p.roi_coordinates{23}(j,3),MarkerSize,'MarkerEdgeColor',c(j,:),'MarkerFaceColor',c(j,:));
end



%%
%{
p = parcellation('dk_full');
tf = D(:,:,:);
t = -1.8:.05:1.8-0.05;
f = 5:1:145;
f0 = nearest(f,60);
f1 = nearest(f,90);
t0 = nearest(t,0);
t1 = nearest(t,0.5);
coor = p.roi_coordinates{23};
res = 8;
gmin = min(coor);
gmax = max(coor);
[x y z] = meshgrid(gmin(1):res:gmax(1), gmin(2):res:gmax(2), gmin(3):res:gmax(3));
[~, index] = ismember( coor, [x(:) y(:) z(:)], 'rows'); % find indices in grid
[s1, s2, s3] = size(x);
TF = nan(s1*s2*s3,length(f), length(t));
TF(index,:,:) = tf;
TF = reshape(TF, [s1,s2,s3,length(f), length(t)]);
TFsmall = TF(:,:,:, f0:f1, t0:t1);
%}
if ~exist('remove_parc', 'var'), remove_parc = 1; end
PATH  = '/Volumes/T5_OHBA/analysis/HMM-gamma/TF/M1/';
if remove_parc
  d = dir([PATH, 'rmptf_sel*.mat']);
else
  d = dir([PATH, 'rmptf_efd*.mat']);
end

if ~exist('subs', 'var'), subs = 1:33; end
p = parcellation('dk_full');
time_bl     = [-1 -0.5];

for s=subs
D = spm_eeg_load([PATH, d(s).name]);


tf = D(:,:,:);
t = D.time;
f = D.frequencies;
f0 = nearest(f,60);
f1 = nearest(f,90);
fS = f(f0:f1);
t0 = nearest(t,0);
t1 = nearest(t,0.5);
tS = t(t0:t1);
tfS = tf(:,f0:f1, t0:t1);
tfS_avT = nanmean(tfS,3);

% max voxel
[voxel,col]=find(tfS_avT==max(tfS_avT(:)));
relpower = tfS_avT(voxel,col);
peakfreq = fS(col);
tf_2d = squeeze(tf(voxel, f0:f1,:));

if remove_parc
  datarank = 55;
  dip_index_sorted = [36,38,52,54,26,68,50,67,37,39,55,71,42,51,57,27,79,53,70,49,78,41,72,83,66,58,69,64,43,82,73,89,56,28,77,84,40,17,81,92,65,63,85,80,88,75,91,59,74,90,76,18,44,87,29,31,86,11,20,12,9,61,30,21,32,46,19,15,33,22,24,47,62,13,34,60,45,35,2,25,10,16,23,48,4,14,7,3,1,5,8,6];
  dip_index = sort(dip_index_sorted(1:datarank));
else
  dip_index = 1:92;
end

tmp = extractBetween(d(s).name, 'efd_', '_go');
if remove_parc
  filename = [PATH,'peaks_60_90Hz/', tmp{1},'_sel'];
else
  filename = [PATH,'peaks_60_90Hz/', tmp{1}];
end
save([filename, '.mat'],'voxel','peakfreq','relpower', 'dip_index', 'tf_2d');

h = figure; h.WindowState = 'maximized';
% subplot(4,2,[1 3]);hold on;
subplot(2,2,1); hold on
plot(fS,tfS_avT,'color',[0.4 0.4 0.4]);
plot(fS,tfS_avT(voxel,:)','r','Linewidth',2);
plot([peakfreq peakfreq],[min(tfS_avT(:)) max(tfS_avT(:))],'m--','Linewidth',2);
grid on;axis tight;
xticks([f(1):1:f(end)]);
title([{d(s).name,'Power per frequency per M1 Voxel between movement on and -offset'}],'interp','none');
ylabel('Power');xlabel('Frequency (Hz)');
set(gca,'Fontsize',13);

cmap = flipud(brewermap(64, 'RdBu'));
colormap(cmap)
% subplot(4,2,[2 4]);hold on;
subplot(2,2,[3 4]); hold on;
imagesc(t,fS,tf_2d);
rectangle('Position',[time_bl(1) fS(1) time_bl(2)-time_bl(1) fS(end)-fS(1)],'Linestyle','--','Linewidth',2);

plot([0 0],[fS(1) fS(end)],'k--','Linewidth',2);
plot([t(1) t(end)],[peakfreq peakfreq],'m--','Linewidth',2);
set(gca,'YDir','normal');
axis tight;
cb=colorbar;cb.Position=[0.92 0.15 0.01 0.3];
cc = max(abs([min(min(tf_2d)),max(max(tf_2d))]));caxis([-cc cc]);
xlabel('Time (s)');ylabel('Frequency (Hz)');
title(['Avg of all epochs for the selected voxel'],'interp','none');
set(gca,'Fontsize',13);

% subplot(4,2,6);
% % hist(MT/1000);
% xlim([D.time(1) D.time(end)]);
% grid on;xlabel('Time (s)');ylabel('MT counts');
% set(gca,'Fontsize',13);

% subplot(4,2,[5 7]);hold on;
subplot(2,2,2); hold on
steps = 6;
cmap = flipud(brewermap(2*steps, 'RdBu'));
cmap = cmap(steps+1:end,:);
S = mean(tfS_avT,2);
xx = linspace(0, max(S(:)), steps);
for j=1:numel(dip_index)
  B = xx-S(j);
  B(B<0)=inf;
  [~, ixtmp] = min(B);
  ix(j) = ixtmp;
end
for j=1:numel(dip_index)
  c(j,:) = cmap(ix(j),:);
end

sc=scatter3(p.template_coordinates(:,1),p.template_coordinates(:,2),p.template_coordinates(:,3),8,'MarkerEdgeColor',[.4 .4 .4],'MarkerFaceColor',[.4 .4 .4]);
sc.MarkerFaceAlpha = 0.1;sc.MarkerEdgeAlpha = 0.1;
for j=1:numel(dip_index)
  MarkerSize = 40;
  scatter3(p.roi_coordinates{23}(dip_index(j),1),p.roi_coordinates{23}(dip_index(j),2),p.roi_coordinates{23}(dip_index(j),3),MarkerSize,'MarkerEdgeColor',c(j,:),'MarkerFaceColor',c(j,:));
  if j==voxel
    scatter3(p.roi_coordinates{23}(dip_index(voxel),1),p.roi_coordinates{23}(dip_index(voxel),2),p.roi_coordinates{23}(dip_index(voxel),3),70,'go');
  end
end
view(-90,0);
axis square;axis tight;
zlabel('I to S');ylabel('A to P');xlabel('L to R');
title('green encircled voxel is the maximum voxel')
%p.plot;%gives glassbrain, but cant be a subplot

print('-dpng',[filename '.png']);pause(1);
close;
end

