% prepare a matrix containing ANITA data
% Input: num represents one of 24 available events
% Output: two matrices [observations, impulse_responses], where rows are signals
%	ax will be of size 1561 and aximp is of size 1024
function [ax, aximp] = prepsig(num)


datapath = 'data/waveforms/info_';

impulse_pathname = 'data/notches_imp_resp/';

% naming of the observation signals: ax(i,:)
% naming of the impulse response: aximp(i,:)
ax = load(strcat(datapath,num2str(num),'_sig'))';
antlist = load(strcat(datapath,num2str(num),'_ant'));

aximp = zeros(length(antlist),1024);
for i = 1:length(antlist)
	aximp(i,:) = load(strcat(impulse_pathname, getantenna(antlist(i))))(:,2)';
end
