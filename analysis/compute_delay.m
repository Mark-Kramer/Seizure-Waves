function [ delay, delay_ci_lo, delay_ci_up ] = compute_delay( coh, coh_conf, phi, freq )
%COMPUTE_DELAY Summary of this function goes here
%   Detailed explanation goes here

delay = zeros(size(coh,1));
delay_ci_lo = zeros(size(coh,1));
delay_ci_up = zeros(size(coh,1));

for i = 1 : size(coh,2)
    for j = i+1 : size(coh,2)
        fprintf('Compute delay electrodes %d - %d...\n', i, j);
        coh0 = squeeze(coh(i,j,:));
        conf0 = coh_conf(i,j);
        phi0 = squeeze(phi(i,j,:));
        
        phi0(coh0 < conf0) = NaN; % find indices where coherence not sign. and set bad indices to NaN.    
        [delay(i,j), delay_ci_lo(i,j), delay_ci_up(i,j)] = delay_from_phase(phi0, freq);
        delay(j,i) = -delay(i,j);
    end
end

end


function [delay, delaylo, delayhi] = delay_from_phase(phi, freq)
%INPUTS
%  phi = phase of coherence at each freq.
%  freq = frequency axis.

    FREQ_INTERVAL = 3;      % we demand 3 Hz sequence of significant coherence to compute delay.
    PVALUE_THRESH = 0.05;   % threshold of the p-value for the linear fit
    
    df = freq(2)-freq(1);             % define frequency resolution (from coherence calculation).
    N_pts_in_sequence = ceil(FREQ_INTERVAL/df);
    delay = NaN;
    delaylo = NaN;
    delayhi = NaN;
    
    % find sequences of repeated finite values in phi
    v = [false isfinite(phi(:)') false]; % add false the ensure a proper grpS and grpE
    grpS = find(~v(1 : end - 1) & v(2 : end));
    grpE = find(v(1 : end - 1) & ~v(2 : end)) - 1;
    grpN = grpE - grpS + 1;
    [mx, imx] = max(grpN);                              % find longest stretch of 1's (i.e., longest stretch of finite values).
    if mx > N_pts_in_sequence                           % need consecutive run that is long enough in frequency.
        phi_to_fit = unwrap(phi(grpS(imx):grpE(imx)));  % get phase in the longest sequence.
        f_to_fit  =  freq(grpS(imx):grpE(imx));         % get frequency in the longest sequence.
        X = [ones(mx, 1), f_to_fit'];                   % establish predictors (include constant),
        [b,bint,~,~,stats] = regress(phi_to_fit,X);    % perform linear regression,
        delay = b(2) / (2*pi);                          % the delay is the slope of fit, then convert from rad to sec
        delaylo = bint(2,1) / (2*pi);                   % ... lower boundary of 95% confidence interval.
        delayhi = bint(2,2) / (2*pi);                   % ... upper boundary of 95% confidence interval.
        mp = stats(3);                                  % ... p-value of fit.
        
        % find p-values of fit that are too big and set to NaN.
        if mp > PVALUE_THRESH
            delay = NaN;
            delaylo = NaN;
            delayhi = NaN;
        end
    end
    
end
