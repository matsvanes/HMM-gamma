PATH_BASE = '/Volumes/T5_OHBA/analysis/HMM-gamma/';
PATH  = [PATH_BASE 'TF/M1/optimised/'];
PATH_TIMING = [PATH_BASE, 'timing/'];

if ~exist('remove_parc', 'var'), remove_parc = 0;                            end
if ~exist('dosave', 'var'),      dosave      = 1;                            end
if ~exist('bltype', 'var'),      bltype      = 'rel';                        end
if ~exist('subs', 'var'),        subs        = 1:33;                         end
if ~exist('doplot', 'var'),      doplot      = 1;                            end
if ~exist('prefix', 'var'),      prefix      = input('what is the prefix?'); end
subinfo;

p = parcellation('dk_full');
time_bl     = [-.5 -.1];%[-1 -0.5];
time_peri   = [0 0.5];
datarank = 55;

% make subselection of dipole locations
if remove_parc
  sel = 'sel_';
  load([PATH, prefix, '_dip_index_sorted.mat']);
  dip_index = sort(dip_index_sorted(1:datarank));
else
  sel=[];
  dip_index = 1:92;
end

tf_avT_avE_group = nan(33, length(dip_index), 31);
for s=subs
  d = dir([PATH, sprintf('ptf_%s%s*%s*.mat', sel, prefix, sub(s).id)]);
  s
  D = spm_eeg_load([PATH, d.name]);
  t = D.time;
  t0 = nearest(t, -1.8); 
  t1 = nearest(t, 1.75);
  f = D.frequencies;
  f0 = nearest(f,60);
  f1 = nearest(f,90);
  fS = f(f0:f1);
  tf = D(:,f0:f1,t0:t1,:);
  chanlabels = D.chanlabels;
  
  % load timing of individual trials
  tmp = extractBetween(d.name, 'case_', '_go');
  tmp = dir([PATH_TIMING, 'emgf*', sprintf('%s*', tmp{1})]);
  load([PATH_TIMING tmp.name, '/Timings.mat']);
  
  % base line correct
  t1 = nearest(t, time_bl(1));
  t2 = nearest(t,time_bl(2));
  avB =       repmat(squeeze(mean(mean(tf(:,:,t1:t2,:),3),4)), [1 1 size(tf,3), size(tf,4)]);
  switch bltype
    case 'diff'
      tf_blc = tf-avB;
    case 'rel'
      tf_blc = (tf./avB)-1;
  end
  
  [s1,s2,s3,s4] = size(tf_blc);
  tf_blc_avT = nan(s1,s2,s4);
  for e=1:size(tf_blc,4)
    clear MT_start_ind MT_end_ind
    MT_start_ind = nearest(t,(time_peri(1))); %movement onset
    MT_end_ind = nearest(t, (time_peri(1)+(MT(e)/1000))); % movement offset
    tf_blc_avT(:,:,e)=squeeze(mean(tf_blc(:,:,MT_start_ind:MT_end_ind,e),3)); % av over time (avT)
    tf_avT(:,:,e)=squeeze(mean(tf(:,:,MT_start_ind:MT_end_ind,e),3)); % av over time (avT)
  end
  tf_blc_avT_avE=squeeze(nanmean(tf_blc_avT,3)); % av over epochs
  tf_avT_avE_group(s,:,:)=tf_blc_avT_avE; % av over epochs
  
  % max voxel in spectrospatial space
  [voxel,col]=find(tf_blc_avT_avE==max(tf_blc_avT_avE(:)));
  power = tf_blc_avT_avE(voxel,col);
  peakfreq = fS(col);
  tf_2d = squeeze(mean(tf_blc(voxel,:,:,:),4));
  tf_3d = squeeze(tf_blc(voxel,:,:,:));
  tf_2d_group(:,:,s) = tf_2d;
  if doplot
  tmp = extractBetween(d.name, sprintf('%s_', prefix), '_go');
  if remove_parc
    filename = [PATH,'peaks_60_90Hz/', prefix, '_', tmp{1},'_sel'];
  else
    filename = [PATH,'peaks_60_90Hz/', prefix,'_', tmp{1}];
  end
  if dosave
    save([filename, '.mat'],'voxel','peakfreq','power', 'dip_index', 'tf_2d', 'tf_3d');
  end
  
  h = figure;% h.WindowState = 'maximized';
  % subplot(4,2,[1 3]);hold on;
  subplot(2,2,1); hold on
  plot(fS,tf_blc_avT_avE,'color',[0.4 0.4 0.4]);
  plot(fS,tf_blc_avT_avE(voxel,:)','r','Linewidth',2);
  plot([peakfreq peakfreq],[min(tf_blc_avT_avE(:)) max(tf_blc_avT_avE(:))],'m--','Linewidth',2);
  grid on;axis tight;
  xticks([f(1):1:f(end)]);
  title([{d.name,'Power per frequency per M1 Voxel between movement on and -offset'}],'interp','none');
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
  cc = max(abs(tf_2d(:)));caxis([-cc cc]);
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
  tf_blc_avF_avT_avE = mean(tf_blc_avT_avE,2);
  xx = linspace(0, max(tf_blc_avF_avT_avE(:)), steps);
  for j=1:numel(dip_index)
    B = xx-tf_blc_avF_avT_avE(j);
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
  if dosave
    saveas(h, [filename], 'png');pause(10);
  end
end
  close; clear tf_blc_avT_avE h avB tf_blc tf_avT tf_blc_avT tf_2d tf_3d tf_blc_avF_avT_avE tf stat D
end
%%
% plot group stats
dat = [];
dat.powspctrm = tf_avT_avE_group;
dat.dimord = 'rpt_chan_freq';
dat.freq = fS;
dat.label = chanlabels;

datB = rmfield(dat, 'powspctrm');
datB.powspctrm = dat.powspctrm*0;

osl_shutdown;
addpath('/Volumes/T5_OHBA/software/fieldtrip/');

cfg=[];
cfg.frequency = [60 90];
cfg.design = [1:33 1:33; ones(1,33) 2*ones(1,33)];
cfg.ivar = 2;
cfg.uvar = 1;
cfg.method = 'analytic';
cfg.statistic = 'depsamplesT';
stat = ft_freqstatistics(cfg, dat, datB);

osl_startup()
T=stat.stat;
T=max(T,[],2); % select maximum frequency for each location

[T_sort, dip_index_sorted] = sort(T,1, 'descend');
if remove_parc
  filename = [PATH, sprintf('%s_',prefix), 'dip_index_sorted_sel.mat'];
else
  filename = [PATH, sprintf('%s_',prefix), 'dip_index_sorted.mat'];
end
if dosave
  save(filename, 'dip_index_sorted', 'T_sort', 'stat', 'tf_2d_group')
end

M1_idx = dip_index_sorted(1:datarank);

figure;
subplot(1,2,2);
tf_2d_groupav = mean(tf_2d_group,3);
imagesc(t, fS,tf_2d_groupav');
cmap = flipud(brewermap(64, 'RdBu'));
colormap(cmap)
title('group average TF for subject specific locations')
rectangle('Position',[time_bl(1) fS(1) time_bl(2)-time_bl(1) fS(end)-fS(1)],'Linestyle','--','Linewidth',2);

set(gca,'YDir','normal');
axis tight;
cb=colorbar;cb.Position=[0.92 0.15 0.01 0.3];
cc = max(abs(tf_2d_groupav(:)));caxis([-cc cc]);
xlabel('Time (s)');ylabel('Frequency (Hz)');

% plot the average T-value for each location
steps = 6;
cmap = flipud(brewermap(2*steps, 'RdBu'));
cmap = cmap(steps+1:end,:);
xx = linspace(0, max(T(:)), steps);
for j=1:numel(dip_index)
  B = xx-T(j);
  B(B<0)=inf;
  [~, ixtmp] = min(B);
  ix(j) = ixtmp;
end
for j=1:numel(dip_index)
  c(j,:) = cmap(ix(j),:);
end
subplot(1,2,1); hold on
sc=scatter3(p.template_coordinates(:,1),p.template_coordinates(:,2),p.template_coordinates(:,3),8,'MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5]);
sc.MarkerFaceAlpha = 0.1;sc.MarkerEdgeAlpha = 0.1;
for j=1:numel(dip_index)
  if any(M1_idx==j)
    MarkerSize = 30;
  else
    MarkerSize = 15;
  end
  scatter3(p.roi_coordinates{23}(dip_index(j),1),p.roi_coordinates{23}(dip_index(j),2),p.roi_coordinates{23}(dip_index(j),3),MarkerSize,'MarkerEdgeColor',c(j,:),'MarkerFaceColor',c(j,:));
end
view(-90,0);
axis square;axis tight;
zlabel('I to S');ylabel('A to P');xlabel('L to R');
title('maximum T per selected per location - group level spatiospectral stat')

if remove_parc
  filename = [PATH, sprintf('%s_groupT_dipole_select_sel', prefix), '.png'];
else
  filename = [PATH, sprintf('%s_groupT_dipole_select', prefix), '.png'];
end
if dosave
  print('-dpng',filename);pause(1);
end

