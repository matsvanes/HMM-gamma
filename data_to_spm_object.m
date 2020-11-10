function D_out = data_to_spm_object(D0, D_in, filename)
% move data into spm meeg object without any montages

if strcmp(class(D_in), 'meeg')
  if ndims(D_in)==2
    D_in = D_in(:,:);
  elseif ndims(D_in)==3
    D_in = D_in(:,:,:);
  elseif ndims(D_in)==4
    D_in = D_in(:,:,:,:);
  end
end
if ndims(D_in)==2
  dim = [size(D_in) 1];
else
  dim=size(D_in);
end
D_out = clone(montage(D0,'switch',0),filename,dim);
D_out = D_out.montage('remove', 1:D_out.montage('getnumber'));

if ndims(D_out)==2
  D_out(:,:) = D_in;
elseif ndims(D_out)==3
  D_out(:,:,:) = D_in;
elseif ndims(D_out)==4
  D_out(:,:,:,:) = D_in;
end