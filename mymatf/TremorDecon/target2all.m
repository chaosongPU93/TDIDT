function targetind = target2all(datavec,sepran)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % targetind = target2all(datavec,sepran)
  %
  % This function can be regarded as a generalized version of 'srcdistall.m'
  % which tries to obtain distance of each source to all others that are 
  % separated within some time. This code simply looks at one 'datavec', and
  % identify indices of target members for each member whose separation in 
  % 'datavec' is within 'sepran'
  %
  %
  % Chao Song, chaosong@princeton.edu
  % First created date:   2023/06/12 
  % Last modified date:   2023/06/12 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  
  defval('sepran',[]);
  
  nsrc = size(datavec,1);
  
  targetind = cell(nsrc, 1); %differential time, since it is sorted by time, it is always positive

  %loop for every source
  for i = 1: nsrc-1
    othersrc = i+1:nsrc;
    tempdata = datavec(othersrc);
    if ~isempty(sepran)
      ind = find(abs(tempdata-datavec(i))>=sepran(1) & abs(tempdata-datavec(i))<=sepran(2));   %closer in time, if not all
    else
      ind = 1:length(othersrc);
    end
    targetind{i} = othersrc(ind);
    
  end