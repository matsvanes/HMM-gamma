function [Gamma, MLGamma, dynamics, spectra, spectramt, tf, tfmt, tf_avg, tfmt_avg] = hmm_gamma_hmm_post(X, Gamma, hmm, T, time, options, doplot, filename)
if ~exist('doplot', 'var'), doplot=1; end
if contains(filename, '.mat'), filename = extractBefore(filename, '/mat'); end

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

try
% State spectra
if options.order>0
  spectra = hmmspectramar(cat(1,X{:}), cat(1,T{:}), hmm, Gamma, options);
else
  spectra = hmmspectratde(hmm, options);
end
spectramt = hmmspectramt(cat(1,X{:}), cat(1,T{:}), Gamma, options);

% reshape Gamma into ntime*ntrl*nstates
Gamma = reshape(Gamma,T{1}(1)-options.order-numel(options.embeddedlags)+1,[],options.K);

% Time-frequency decomposition
% tf = hmmtimefreq(spectra, Gamma_pad, 0); % based on HMM parameters
tfmt = hmmtimefreq(spectramt, Gamma_pad, 0); % based on multitapers

% psd per subject
[tf, tf_avg] = tftosub(tf, T, false);
[tfmt, tfmt_avg] = tftosub(tfmt, T, false);


if doplot
  hmm_post_plot(options, T, time, dynamics.F0, MLGamma, spectra, spectramt, tf_avg, tfmt_avg);
  if ~exist('filename', 'var') || isempty(filename)
    filename = 'HMM';
  end
  print('-dpng',filename);pause(1);
end

catch
  Gamma = reshape(Gamma,T{1}(1)-options.order-numel(options.embeddedlags)+1,[],options.K);
  spectra=[]; spectramt=[]; tf=[]; tfmt=[]; tf_avg=[]; tfmt_avg=[];
end
end