function tf_blc_avT_avE = tf_trialspecific_timeavg(tf_blC, timing, t)
% In a TF spectrum of dim (chan)*freq*time*trial, average over the movement
% time (trial specific) and then over trials. TF spectrum should be
% baseline corrected.

time_peri   = [0 0.5];
if ndims(tf_blC)==4
  [nchan, nfreq, ntim, ntrl] = size(tf_blC);
  tf_blC = reshape(tf_blC, [nchan*nfreq, ntim, ntrl]);
  chanflag=1;
  dim1 = nchan*nfreq;
else
  [dim1, ntim, ntrl] = size(tf_blC);
  chanflag=0;
end


tf_avT= zeros(dim1,ntrl);
for e=1:length(timing)
  clear MT_start_ind MT_end_ind
  start_ind = nearest(t,(time_peri(1))); %movement onset
  end_ind = nearest(t, (time_peri(1)+(timing(e)/1000))); % movement offset
  tf_avT(:,e)=squeeze(mean(tf_blC(:,start_ind:end_ind,e),2)); % av over time (avT)
end
  tf_blc_avT_avE = nanmean(tf_avT,2);
  
if chanflag
  tf_blc_avT_avE = reshape(tf_blc_avT_avE, [nchan,nfreq]);
end