function [tf, tf_avg] = tftosub(tf, T, keeptrials)
% intended to restructure the HMM's tf spectra (time * freq) to a
% time*freq*trial matrix per subject
%
% NOTE: does not yet work for multi-channel input
%
% Input
% tf: matrix time * freq, output of hmmtimefreq
% T: cell (#sub) of amount of samples per trial
% keeptrials: boolean, whether to keep individual trials or subject average
%
% Output
% tf: cell (#sub), matrix of time*freq(*trial)
% tf_avg: matrix time*freq group average tf

if ndims(tf)>2, error('this function does not yet work on multichannel data'), end
if nargin<3 || isempty('keeptrials'), keeptrials=false; end

for ii=1:numel(T)
  t0 = sum(cat(1,T{1:ii}))-sum(cat(1,T{ii}));
  tf_tmp{ii} = zeros(max(cat(1,T{:})), size(tf,2), length(T{ii}));
  for jj=1:length(T{ii})
    tf_tmp{ii}(:,:,jj)=tf(t0+sum(T{ii}(1:jj))-sum(T{ii}(jj)) +1 : t0+ sum(T{ii}(1:jj)),:);
  end
  tf_avgTr{ii} = nanmean(tf_tmp{ii},3);
end
if keeptrials
  tf = tf_tmp;
else
  tf = tf_avgTr;
end
tf_avgTr = permute(cat(3, tf_avgTr{:}), [3 1 2]);

% group TF
tf_avg = squeeze(nanmean(tf_avgTr,1));