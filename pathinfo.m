
islocal = isfolder('/Volumes/T5_OHBA/');
if islocal
  projectdir = '/Volumes/T5_OHBA/';
  homedir = '/Volumes/T5_OHBA/';
  startupdir = '/Users/matsvanes/Documents/MATLAB/startup.m';
else
  projectdir = '/ohba/pi/mwoolrich/mvanes/';
  homedir = '/home/mvanes/';
  startupdir = [homedir, 'startup.m'];
end
  fieldtrippath = [homedir, 'software/fieldtrip/'];