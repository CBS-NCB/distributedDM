% change protocol to chose specific trials
Dir = 'C:\Users\Mohammad-BL\Google Drive\BL\MPEP_LOGS\M160608_MA\1\2\';
load([Dir 'Protocol.mat']);
nrep = 19;
Protocol.seqnums(:,nrep+1:end)=[];
Protocol.nrepeats=nrep;
save([Dir 'Protocol.mat'],'Protocol')
