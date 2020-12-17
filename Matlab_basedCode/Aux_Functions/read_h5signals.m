function read_h5signals(sigName,sigNameInput)
%read_h5signals(sigName)

if nargin<2
    sigNameInput='signals.h5';
end
    
info     = h5info(sigNameInput);
channels = numel(info.Groups(1).Datasets);
t = h5read (sigNameInput, sprintf('/signal/times'));
[t,I]=sort(t);
for i=1:channels
     A = ( h5read (sigNameInput, sprintf('/signal/channels/%02i',i-1)) )';
     M(:,:,i) = A(I,:);
end
d=diff(diff(t));d(d<1e-5)=0;I=[0;find(d)+1;length(t)];
save(sigName,'t','M','I');
