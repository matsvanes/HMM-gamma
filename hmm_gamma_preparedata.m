function [D, POI, dat] = hmm_gamma_preparedata(source_path, target_path, filename, roi, remove_parc)
S = [];
S.D=spm_eeg_load(fullfile(source_path, filename));
if (strcmp(roi, 'parc') || strcmp(roi, 'M1')) && remove_parc
  S.outfile = [target_path, 'sel_', filename];
else
  S.outfile = fullfile(target_path, filename);
end
D = spm_eeg_copy(S);
if strcmp(roi, 'M1') || strcmp(roi, 'parc')
  % orthogonalise
  switch roi
    case 'parc'
      use_montage = 6;
      D=D.montage('switch',use_montage);
      dat_org = D(:,:,:);
      dat = reshape(dat_org,size(dat_org,1),[]);
      datarank = rank(dat);
      if remove_parc % or come up with a way to reduce the number of parcels to the rank of the data
        parcels_to_be_del = [1 4 5 6 11 15 17 18 19 34 35 36 37 38 39 40 41 1+41 4+41 5+41 6+41 11+41 15+41 17+41 18+41 19+41 34+41 35+41 36+41 37+41 38+41 39+41 40+41 41+41];
        dat(parcels_to_be_del,:) = [];
        POI = 14;
      else
        POI = 23;
      end
    case 'M1'
      use_montage = 5;
      D = D.montage('switch', use_montage);
      p = parcellation('dk_full');
      M1idx = find(contains(p.labels, 'Left Precentral'));
      dip_index = p.parcelflag;
      dip_index = find(dip_index(:,M1idx)); % find the dipole locations in the M1 parcel
      dat_org = D(dip_index,:,:);
      dat = reshape(dat_org,size(dat_org,1),[]);
      datarank = rank(dat);
      if remove_parc
        % dip_index_sorted is based on the group level T-stats in the [60
        % 90] Hz range and [0 0.5]s window.
        load([target_path, 'dip_index_sorted.mat']);
        dip_index_sorted = sort(dip_index_sorted(1:datarank));
        dat = dat(dip_index_sorted,:);
      end
      warning('POI for M1 not yet supported: first the maximum gamma location should be selected')
      POI = [];
  end
  % remove source leakage. Only possible if certain parcels/dipoles were
  % removed.
  if size(dat,1)==datarank
    dat = ROInets.remove_source_leakage(dat,'symmetric');
  end
  dat = reshape(dat,size(dat,1),D.nsamples,D.ntrials);
  outfile = fullfile(target_path,D.fname);
  Dnode = clone(montage(D,'switch',0),outfile,[size(dat,1),D.nsamples,D.ntrials]);
  Dnode = chantype(Dnode,1:Dnode.nchannels,'VE');
  Dnode(:,:,:)=dat;
  D = Dnode;
else
  use_montage = 2;
  D = D.montage('switch', use_montage);
  dat = D(:,:,:);
  Dtmp = D.montage('remove', 1:6);
  Dtmp(:,:,:)=dat;
  D=Dtmp;
end
D.save;