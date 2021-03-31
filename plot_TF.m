function h = plot_TF(t, f, tf, titl, clim, cmap)

if numel(t)==size(tf,1) && numel(f)==size(tf,2)
  tf=tf';
end

if ~exist('cmap', 'var') || isempty(cmap)
  cmap = flipud(brewermap(64, 'RdBu'));
end
if ~exist('titl', 'var'), titl=[]; end

if exist('clim', 'var') && ~isempty(clim)
  imagesc(t, f, tf, clim);
else
  imagesc(t, f, tf);
end
colormap(cmap)
set(gca,'YDir','normal');
xlabel('Time (s)');ylabel('Frequency (Hz)');
title(titl)
colorbar


