function plotPursuitCurves(P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot pursuit curves. 
% simulates and plots pursuit curves for any collection of starting points (mice).
%   Start by constructing P, a vector of four complex-valued points, so that they
%   form the corners of a square in the complex plane.
%
% EXAMP 1: a square
%   P = [0 1 1+1i 1i];
%   clf
%   plotPursuitCurves(P)
%   axis([-0.1 1.1 -0.1 1.1])
%
% EXAMP 2: a hexagon
%   clf
%   theta = 0:60:330;
%   P = cosd(theta) + 1i*sind(theta);
%   plotPursuitCurves(P)
%   axis([-1.1 1.1 -1.1 1.1])
%
% From Steve Eddins's blog

P = reshape(P,1,[]);
N = length(P);
d = 0.01;
for k = 1:1000
    V = P(end,[(2:N) 1]) - P(end,:);
    P(end+1,:) = P(end,:) + d*V;
end

hold on
Q = P(:,[(2:N) 1]);
for k = 1:10:size(P,1)
    for m = 1:N
        T = [P(k,m) Q(k,m)];
        plot(real(T),imag(T),'LineWidth',0.5,'Color',[.7 .7 .7])
    end
end

for k = 1:N
    plot(real(P(:,k)),imag(P(:,k)),'LineWidth',1.5)
end
hold off
axis equal
end


