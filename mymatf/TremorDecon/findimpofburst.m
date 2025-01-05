function impi = findimpofburst(imp,nsrc,ibst)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % impi = findimpofburst(imp,nsrc,ibst)
  %
  % Imagine you have all sources 'imp' and the num of sources for each burst 
  % 'nsrc', and now you need the sources 'impi' from a particular burst 'ibst'
  % 
  %
  % 
  % By Chao Song, chaosong@princeton.edu
  % First created date:   2024/11/20
  % Last modified date:   2024/11/20
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ist = sum(nsrc(1:ibst-1))+1;
  ied = ist+nsrc(ibst)-1;
  impi = imp(ist: ied, :);
  