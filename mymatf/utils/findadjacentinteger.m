function adjint = findadjacentinteger(N,Ns)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Let's say I already have Ns random integers within [1, N], N>Ns, now I want
% to find all integers that are right next to the existing numbers. For example,
% if I have 1, 4, 8 within [1,10], then the code needs to return 2, 3, 5, 7 and 
% 9, because 2 is next to 1, 3 and 5 are next to 4, 7 and 9 are next to 8.
% --To find all integers that are right next to the existing numbers within a 
% range in MATLAB, you can use the following approach:
% 
%   1. Define your range of integers.
%   2. Identify the numbers that are adjacent to your given numbers.
%   3. Make sure these adjacent numbers fall within the range and are unique.
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2024/07/07
% Last modified date:   2024/07/07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize an empty array to store adjacent numbers
adjint = [];

% Loop through each number in Ns to find adjacent numbers
for i = 1:length(Ns)
    current_number = Ns(i);
    
    % Check if the current number is not the first in the range
    if current_number > 1
        adjint(end+1) = current_number - 1;
    end
    
    % Check if the current number is not the last in the range
    if current_number < N
        adjint(end+1) = current_number + 1;
    end
end

% Remove duplicates from the adjacent_numbers array
adjint = unique(adjint);

% Remove numbers that are in Ns from the adjacent_numbers array
adjint = setdiff(adjint, Ns);

% % Display the result
% disp('The adjacent numbers are:');
% disp(adjint);
