
islocal = isfolder('/Volumes/T5_OHBA/');
if islocal
  projectdir = '/Volumes/T5_OHBA/';
  homedir = '/Volumes/T5_OHBA/';
else
  projectdir = '/ohba/pi/mwoolrich/mvanes/';
  homedir = '/home/mvanes/';
end
  fieldtrippath = [homedir, 'software/fieldtrip/'];