function [MLGamma, dynamics, spectra, tf] = hmm_gamma_hmm_post(X, Gamma, hmm, T, options, doplot, filename)
if ~exist('doplot', 'var'), doplot=1; end

% Post process HMM
if iscell(hmm)
  hmm = hmm{1};
end

% Pad Gamma until original length
Gamma_pad = padGamma(Gamma,T,options);
time = -2:1/500:2;

% Measures of state dynamics
dynamics=[];
dynamics.F0 = getFractionalOccupancy(Gamma, T, options, 1);% state fractional occupancies per session
dynamics.LifeTimes = getStateLifeTimes(Gamma,T,hmm.train); % state life times
dynamics.Intervals = getStateIntervalTimes (Gamma,T,hmm.train); % interval times between state visits
dynamics.SwitchingRate =  getSwitchingRate(Gamma,T,hmm.train); % rate of switching between states

% Trial time courses
MLGamma = gammaRasterPlot(Gamma_pad, cat(1,T{:}));
MLGamma = reshape(MLGamma,T{1}(1),[]);

% State spectra
spectra = hmmspectramar(cat(1,X{:}),cat(1,T{:}),hmm, Gamma,options);

% HMM based time-frequency representation
tf = hmmtimefreq(spectra, Gamma_pad, 1);

% psd per subject
for ii=1:numel(T)
  t0 = sum(cat(1,T{1:ii}))-sum(cat(1,T{ii}));
  tf_tmp{ii} = zeros(max(cat(1,T{:})), size(tf,2), length(T{ii}));
  for jj=1:length(T{ii})
    tf_tmp{ii}(:,:,jj)=tf(t0+sum(T{ii}(1:jj))-sum(T{ii}(jj)) +1 : t0+ sum(T{ii}(1:jj)),:);
  end
  tf_avgTr{ii} = nanmean(tf_tmp{ii},3);
end
tf = tf_tmp;
tf_avgTr = permute(cat(3, tf_avgTr{:}), [3 1 2]);

% group TF
tf_avgT_avgS = squeeze(nanmean(tf_avgTr,1));
t1=nearest(time,-0.6); t2=nearest(time,-0.2);

if doplot
  hmm_post_plot(options, T, time, dynamics.F0, spectra, tf_avgT_avgS, MLGamma);
  if ~exist('filename', 'var') || isempty(filename)
    filename = 'HMM';
  end
  print('-dpng',filename);pause(1);
end
% factorisation into frequency bands
%{
if 1 % manual bands
  bands = [0 4; 4 8; 8 13; 14 30; 30 60; 60 inf];
  sp_fit = spectbands(spectra,bands);
else
  % or if using data-driven bands: FIXME:MIGHT NOT WORK? --> use cov
  Gamma = padGamma(Gamma,T,options);
  nsubjects=numel(T);
  sp_fit_subj = cell(nsubjects,1); Tacc = 0;
  for jj = 1:nsubjects
    ind = Tacc + (1:sum(T{jj})); Tacc = Tacc + length(ind);
    sp_fit_subj{jj} = hmmspectramt(X{jj},T{jj},Gamma(ind,:),options);
  end
  params_fac = struct();
  params_fac.Method = 'NNMF';
  params_fac.Ncomp = 4;
  params_fac.Base = 'psd';
  [spectral_factors,spectral_profiles,spectral_factors_subj] = spectdecompose(sp_fit_subj,params_fac);
end
%}


  

end