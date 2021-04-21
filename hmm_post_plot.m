function hmm_post_plot(options, T, time, F0, MLGamma, spectramt, tfmt_avg, time_tf)
    fmt=spectramt.state(1).f;
    figure;
    subplot(2,3,[1 2])
    
    timeidx = options.order+numel(options.embeddedlags);
    time_short = time(timeidx:end); 
    plot(time_short, F0); 
    xlabel('time'), ylabel('occupancy'), title('Fractional occupancy')
    
    subplot(2,3,3)
    imagesc(time, 1:size(MLGamma,2), MLGamma')
    set(gca,'YDir','normal'), xlabel('Time (s)'), ylabel('Trial #'), title('Maximum Likelihood Gamma'), hold on
    for k=1:numel(T)
      hline(numel(cat(1,T{1:k})), 'r');
    end
    xlim([-1 1])

    subplot(2,3,4)
    plot(fmt, cat(2,spectramt.state(:).psd))
    xlabel('Frequency (Hz)'), ylabel('Power'), title('mt State Spectra')
    
    if isfield(options, 'blwindow') && ~isempty(options.blwindow)
      t1 = nearest(time_tf, options.blwindow(1));
      t2 = nearest(time_tf, options.blwindow(2));
      tfmt_avg = tfmt_avg./repmat(nanmean(tfmt_avg(t1:t2,:),1), [size(tfmt_avg,1) 1]);
    end
    subplot(2,3,5)
    plot_TF(time_tf,fmt,tfmt_avg, 'mt TF mtmconvol - full')
    ylim([1 100]), xlim([-1 1])
    
    subplot(2,3,6)
    f1 = nearest(fmt, 60); f2 = nearest(fmt, 90);
    plot_TF(time_tf,fmt(f1:f2),tfmt_avg(:,f1:f2)', 'mt TF mtmconvol - gamma')
    xlim([-1 1])

    colormap(inferno(64))
    suptitle(sprintf('HMM, order %d, %d states, lags %s', options.order, options.K, num2str(options.embeddedlags)))

  end