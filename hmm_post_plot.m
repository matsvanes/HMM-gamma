function hmm_post_plot(options, T, time, F0, spectra, tf_avgT_blC_avgS, MLGamma)
    f=spectra.state(1).f;
    
    figure;
    subplot(2,2,1)
    
    time_short = time; time_short(1:options.order)=[];
    plot(time_short, F0-mean(F0)); 
    xlabel('time'), ylabel('occupancy (demeaned)'), title('Relative fractional occupancy')
    
    subplot(2,2,2)
    plot(spectra.state(1).f, cat(2,spectra.state(:).psd))
    xlabel('Frequency (Hz)'), ylabel('Power'), title('State Spectra')
    
    subplot(2,2,3)
    imagesc(time, 1:size(MLGamma,2), MLGamma')
    set(gca,'YDir','normal'), xlabel('Time (s)'), ylabel('Trial #'), title('Trial time courses'), hold on
    for k=1:numel(T)
      hline(numel(cat(1,T{1:k})), 'r');
    end
    
    subplot(2,2,4)
    imagesc(time,f,tf_avgT_blC_avgS')
    set(gca,'YDir','normal'), xlabel('Time (s)'), ylabel('Frequency (Hz)'), title('State TF')
    
    suptitle(sprintf('HMM, order %d, %d states', options.order, options.K))

  end