function impsft=shiftsrctoorigin(imp,nsrc,nbst,off1i,off1iwk,irccrank)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% impsft=shiftsrctoorigin(imp,nsrc,nbst,off1i,off1iwk)
%
% This function to shift the deconvolved impulses which were computed based
% on aligned records back to the origin. For example, the whole-win data 
% detection aligned the whole window and detect around a range within it. So
% we want to shift that alignment back to the origin.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/04/18
% Last modified date:   2024/04/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('off1iwk',[]);
defval('irccrank',[]);

impsft = [];

for i = 1: nbst
  ist = sum(nsrc(1:i-1))+1;
  ied = ist+nsrc(i)-1;
  impi = imp(ist:ied,:);
  off1ii = off1i(i,:);

  impisft = impi;

  if ~isempty(off1iwk) && ~isempty(irccrank)
    off1iw = off1iwk{i};
    irccran = irccrank{i};
    for ip = 1: size(impi,1)
      %find which subwin this ref impulse belongs to
      iwin = findwhichrange(impi(ip,1),irccran);
      %shift the sources grouped relative to the best alignment of that win to the origin
      % impisft(ip,7:8) = impi(ip,7:8) - (off1iw(iwin,2:3) - off1ii(2:3)); %account for prealignment
      impisft(ip,7:8) = impi(ip,7:8) - off1iw(iwin,2:3); %account for prealignment
    end
  else    
    impisft(:,7:8) = impi(:,7:8) - repmat([off1ii(2) off1ii(3)],size(impi,1),1);
  end

  impsft = [impsft; impisft];

end

