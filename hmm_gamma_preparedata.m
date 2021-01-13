function [D, dat] = hmm_gamma_preparedata(PATH, filename, roi, remove_parc, p, dipole_group, timewin)
if ~exist('timewin', 'var'), timewin = []; end
S = [];
S.D=spm_eeg_load(fullfile(PATH.DATA, filename));
if ~isempty(timewin)
  S.timewin = 1000*timewin;
  S.prefix = 'tmp';
  D = spm_eeg_crop(S);
  
  S = [];
  S.D = D;
end
if (strcmp(roi, 'parc') || strcmp(roi, 'M1')) && remove_parc
  S.outfile = [PATH.TARGET, 'sel_', filename];
else
  S.outfile = fullfile(PATH.TARGET, filename);
end
D = spm_eeg_copy(S);
if ~isempty(timewin)
  delete([fullfile(PATH.DATA, ['tmp', filename(1:end-4)]), '*'])
end

if strcmp(roi, 'M1') || strcmp(roi, 'parc')
  % orthogonalise
  switch roi
    case 'parc'
      use_montage = 3;
      D=D.montage('switch',use_montage);
      dat_org = D(:,:,:);
      dat = reshape(dat_org,size(dat_org,1),[]);
      datarank = rank(dat);
      if remove_parc % or come up with a way to reduce the number of parcels to the rank of the data
        parcels_to_be_del = [1 4 5 6 11 15 17 18 19 34 35 36 37 38 39 40 41 1+41 4+41 5+41 6+41 11+41 15+41 17+41 18+41 19+41 34+41 35+41 36+41 37+41 38+41 39+41 40+41 41+41];
        dat(parcels_to_be_del,:) = [];
      end
    case 'M1'
      use_montage = 2;
      D = D.montage('switch', use_montage);
      M1idx = find(contains(p.labels, 'Left Precentral'));
      dip_index = p.parcelflag;
      dip_index = find(dip_index(:,M1idx)); % find the dipole locations in the M1 parcel
      dat_org = D(dip_index,:,:);
      dat = reshape(dat_org,size(dat_org,1),[]);
      datarank = rank(dat);
      if remove_parc
        % dip_index_sorted is based on the group level T-stats in the [60
        % 90] Hz range and [0 0.5]s window.
        dip_index_sorted = sort(dipole_group.dip_index_sorted(1:datarank));
        dat = dat(dip_index_sorted,:);
      end      
  end
  % remove source leakage. Only possible if certain parcels/dipoles were
  % removed.
  if size(dat,1)==datarank
    dat = ROInets.remove_source_leakage(dat,'symmetric');
  end
  dat = reshape(dat,size(dat,1),D.nsamples,D.ntrials);
  outfile = fullfile(PATH.TARGET,D.fname);
  Dnode = clone(montage(D,'switch',0),outfile,[size(dat,1),D.nsamples,D.ntrials]);
  Dnode = chantype(Dnode,1:Dnode.nchannels,'VE');
  Dnode(:,:,:)=dat;
  D = Dnode;
else
  use_montage = 0;
  D = D.montage('switch', use_montage);
  dat = D(:,:,:);
  Dtmp = D.montage('remove', 1:D.montage('getnumber'));
  Dtmp(:,:,:)=dat;
  D=Dtmp;
end
D.save;