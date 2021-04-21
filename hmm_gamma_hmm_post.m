function [Gamma, MLGamma, dynamics, spectramt, tfmt, tfmt_avg, time_tf] = hmm_gamma_hmm_post(X, Gamma, hmm, T, time, options, doplot, filename)
if ~exist('doplot', 'var'), doplot=1; end
if contains(filename, '.mat'), filename = extractBefore(filename, '/mat'); end
if ~isfield(options, 'keeptrials'), options.keeptrials = false; end

% Post process HMM
if iscell(hmm),   hmm  = hmm{1};  end
if iscell(time),  time = time{1}; end
if ~isfield(options, 'embeddedlags') || isempty(options.embeddedlags), options.embeddedlags=0; end

% Pad Gamma until original length
if ndims(Gamma)>2
  Gamma = reshape(Gamma, [], options.K);
end
Gamma_pad = padGamma(Gamma,T,options);

% Measures of state dynamics
dynamics=[];
dynamics.F0 = getFractionalOccupancy(Gamma, T, options, 1);% state fractional occupancies per session
dynamics.LifeTimes = getStateLifeTimes(Gamma,T,hmm.train); % state life times
dynamics.Intervals = getStateIntervalTimes (Gamma,T,hmm.train); % interval times between state visits
dynamics.SwitchingRate =  getSwitchingRate(Gamma,T,hmm.train); % rate of switching between states

% Trial time courses
MLGamma = gammaRasterPlot(Gamma_pad, cat(1,T{:}));
MLGamma = reshape(MLGamma,T{1}(1),[]);

%% spectral
switch options.freqmethod
  case 'hmm_mt'
    % State spectra
    spectramt = hmmspectramt(cat(1,X{:}), cat(1,T{:}), Gamma, options);

    % reshape Gamma into ntime*ntrl*nstates
    Gamma = reshape(Gamma,T{1}(1)-options.order-numel(options.embeddedlags)+1,[],options.K);

    % Time-frequency decomposition
    % tf = hmmtimefreq(spectra, Gamma_pad, 0); % based on HMM parameters
    tfmt = hmmtimefreq(spectramt, Gamma_pad, 0); % based on multitapers

    % psd per subject
    [tfmt, tfmt_avg] = tftosub(tfmt, T, options.keeptrials);
    time_tf = time;
  case 'ft_mt'
    [spectramt, tfmt, tfmt_avg, time_tf] = hmmspectramt_ft(X, T, time, Gamma_pad, options, options.keeptrials);
end

%% plot
if doplot
  hmm_post_plot(options, T, time, dynamics.F0, MLGamma, spectramt, tfmt_avg, time_tf);
  if ~exist('filename', 'var') || isempty(filename)
    filename = 'HMM';
  end
  print('-dpng',filename);pause(1);
end
end
