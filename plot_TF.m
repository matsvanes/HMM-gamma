function h = plot_TF(t, f, tf, titl)

cmap = flipud(brewermap(64, 'RdBu'));
if ~exist('title', 'var'), titl=[]; end

imagesc(t, f, tf);
colormap(cmap)
set(gca,'YDir','normal');
xlabel('Time (s)');ylabel('Frequency (Hz)');
title(titl)
colorbar


