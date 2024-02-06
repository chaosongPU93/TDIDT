function diffmat = diffcustom(datamat,nsep,direction)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diffmat = diffcustom(datamat,nsep)
%
% Sometimes you want to know the difference between the sorted elements in a
% data matrix 'datamat', like N - (N-1), then you can use 'diff'. 'diffcustom' aims 
% to solve similar problem but for non-consecutive elements, like N-(N-2); 
% N-(N-3); N-(N-nsep)...
% You could do this along different 'direction', 'forward' [N-(N-nsep)], 
% or 'backward' [(N+nsep)-N].
% The output 'diffmat' has the length of N-nsep, where N is length of the input
% 'datamat'. 
% Note that the opreation is made to each column of 'datamat'.
%
%
% 
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/06/09
% Last modified date:   2022/06/09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('direction','forward');

% datamat = reshape(datamat,[],1);
N = size(datamat,1);
diffmat = zeros(N-nsep, size(datamat,2));

if isequal(direction,'forward')
  for i = 1+nsep: N 
    diffmat(i-nsep, :) = datamat(i, :)-datamat(i-nsep, :);
  end
elseif isequal(direction,'backward')
  for i = 1: N-nsep 
    diffmat(i, :) = datamat(i+nsep, :)-datamat(i, :);
  end
end