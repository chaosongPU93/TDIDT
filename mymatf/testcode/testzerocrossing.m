t = [1:0.01:5];                                                         % Time Vector
y = sin(2*pi*t);                                                        % Signal
zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);                    % Returns Zero-Crossing Indices Of Argument Vector
zx = zci(coefhf);                                                            % Approximate Zero-Crossing Indices
% figure(1)
% plot(t, y, '-r')
% hold on
% plot(t(zx), y(zx), 'bp')
% hold off
% grid
% legend('Signal', 'Approximate Zero-Crossings')