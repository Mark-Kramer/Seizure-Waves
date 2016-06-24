function [ src_dir, speed, ci_dir, ci_sp ] = estimate_wave( delay, position, varargin )
%ESTIMATE_WAVE Summary of this function goes here
%   Detailed explanation goes here

P_THRESH = 0.05;                                                        % Threshold to accept/reject the fit
MIN_RATIO_FINITE = 0.5;                                                 % we want 50% of electrodes to have a defined delay

src_dir = nan;
speed = nan;

[~, center] = min((position(:,1) - mean(position(:,1))).^2 +...         % find the most central electrode
                  (position(:,2) - mean(position(:,2))).^2);
delay = delay(center,:);                                                % we use delays relative to the center

if length(find(isfinite(delay))) > MIN_RATIO_FINITE * size(position,1)  % check enough delay data is not NaN.
    [b,stats] = robustfit(position, delay, 'fair');
    H = [0 1 0; 0 0 1];  
    c = [0 ; 0];
    P0 = linhyptest(b, stats.covb, c, H, stats.dfe);                    % perform F test that last two coefs are both 0.
    
    if P0 < P_THRESH                                                    % if the fit is significant
        src_dir = squeeze(atan2(b(3), b(2)));                           % compute the source direction
        speed = squeeze(1./sqrt(b(2).^2+b(3).^2));                      % compute velocity mm/s
        samples = mvnrnd(b(2:3), stats.covb(2:3,2:3), 1000);            % generate samples from model.
        boot_dir = atan2(samples(:,2), samples(:,1));                   % bootstrap direction
        boot_sp = 1./sqrt(samples(:,1).^2 + samples(:,2).^2);           % bootstrap velocity.
        ci_dir = quantile(boot_dir, [0.025, 0.975]);                    % CI for direction
        ci_sp = quantile(boot_sp, [0.025, 0.975]);                      % CI for speed
    end
end

if isnan(src_dir)
    fprintf('Unable to fit a plane to the delays in space\n');
elseif nargin ==3 && strcmp(varargin{1}, 'plot')
    scatter3(position(:,1), position(:,2), 1000 * delay, 200, 1000 * delay, 'filled');
    hold on    
    plot3(position(center,1), position(center,2), 0, 'X');
    x1fit = linspace(min(position(:,1)), max(position(:,1)), 10);
    x2fit = linspace(min(position(:,2)), max(position(:,2)), 10);
    [X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
    YFIT = 1000*(b(1) + b(2)*X1FIT + b(3)*X2FIT);
    mesh(X1FIT,X2FIT,YFIT, 'FaceColor', 'interp', 'FaceAlpha', 0.8, 'LineStyle', 'none')
    xlabel('Electrode position (mm)');
    ylabel('Electrodes position (mm)');
    zlabel('Delay (ms)')
    axis tight
end

end

