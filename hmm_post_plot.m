function hmm_post_plot(options, T, time, F0, MLGamma, spectra, tf_avgT_blC_avgS, tfconvol)
    f=spectra.state(1).f;
    
    figure;
    subplot(2,3,1)
    
    timeidx = options.order+numel(options.embeddedlags);
    time_short = time(timeidx:end); 
    plot(time_short, F0); 
    xlabel('time'), ylabel('occupancy'), title('Fractional occupancy')
    
    subplot(2,3,2)
    imagesc(time, 1:size(MLGamma,2), MLGamma')
    set(gca,'YDir','normal'), xlabel('Time (s)'), ylabel('Trial #'), title('Maximum Likelihood Gamma'), hold on
    for k=1:numel(T)
      hline(numel(cat(1,T{1:k})), 'r');
    end
    
    subplot(2,3,4)
    plot(spectra.state(1).f, cat(2,spectra.state(:).psd))
    xlabel('Frequency (Hz)'), ylabel('Power'), title('State Spectra')
    
    subplot(2,3,5)
    imagesc(time,f,tf_avgT_blC_avgS')
    set(gca,'YDir','normal'), xlabel('Time (s)'), ylabel('Frequency (Hz)'), title('State TF spectra'),  ylim([1 100])
    
    subplot(2,3,6)
    imagesc(tfconvol.time,tfconvol.freq,nanmean(cat(3,tfconvol.tf{:}),3)')
    set(gca,'YDir','normal'), xlabel('Time (s)'), ylabel('Frequency (Hz)'), title('State TF mtmconvol'),  ylim([1 100])
    
    suptitle(sprintf('HMM, order %d, %d states, lags %s', options.order, options.K, num2str(options.embeddedlags)))

  end