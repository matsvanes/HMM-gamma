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
scatter3(p.roi_coordinates{23}(j,1),p.roi_coordinates{23}(j,2),p.roi_coordinates{23}(j,3),30,'MarkerEdgeColor',c(j,:),'MarkerFaceColor',c(j,:));
end

[s_sort, idx_sort] = sort(s,1, 'descend');
M1_idx = idx_sort(1:datarank);

