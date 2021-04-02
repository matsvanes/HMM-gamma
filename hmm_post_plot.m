function hmm_post_plot(options, T, time, F0, MLGamma, spectramt, tfmt_avg)
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

    subplot(2,3,4)
    plot(fmt, cat(2,spectramt.state(:).psd))
    xlabel('Frequency (Hz)'), ylabel('Power'), title('mt State Spectra')
    
    subplot(2,3,5)
    imagesc(time,fmt,tfmt_avg')
    set(gca,'YDir','normal'), xlabel('Time (s)'), ylabel('Frequency (Hz)'), title('mt TF mtmconvol - full'),  ylim([1 100])
    
    subplot(2,3,6)
    f1 = nearest(fmt, 60); f2 = nearest(fmt, 90);
    imagesc(time,fmt(f1:f2),tfmt_avg(:,f1:f2)')
    set(gca,'YDir','normal'), xlabel('Time (s)'), ylabel('Frequency (Hz)'), title('mt TF mtmconvol - gamma')
    
    
    suptitle(sprintf('HMM, order %d, %d states, lags %s', options.order, options.K, num2str(options.embeddedlags)))

  end