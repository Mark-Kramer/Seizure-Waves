%% Load the example
% data: SAMPLES x ELECTRODES
% fs: sampling rate
% position: position of the electrodes
load('example_seizure_waves');
time = 1/fs * (0 : 1 : size(data,1) - 1);

figure;
plot(time, data);
xlabel('Time (s)');
ylabel('Voltage (uV)');


%% Compute coherence and phase
CHRONUX_PATH = '/Volumes/Data HDD/Matlab/chronux'; % required toolbox
addpath(genpath(CHRONUX_PATH));

BAND = [1 13];                  % Select a frequency range to analyze
TW = 20;                        % Time-bandwidth product
ntapers = 2*TW-1;               % Choose the # of tapers.
params.tapers = [TW, ntapers];  % ... time-bandwidth product and tapers.
params.Fs = fs;                 % ... sampling rate
params.pad = -1;                % ... no zero padding.
params.fpass = BAND;            % ... freq range to pass
params.err = [1 0.05];          % ... theoretical error bars, p=0.05.

[coh, phi, freq] = compute_coherence(data, params);
save('example_coherence', 'coh', 'phi', 'freq');

figure;
plot(freq, squeeze(10 * log10(coh(1,36,:))));
xlabel('Frequency (Hz)');
ylabel('Coherence');

figure;
subplot(1,2,1)
imagesc(squeeze(coh(:,:, find(freq >= 6, 1))));
c = colorbar;
c.Label.String = 'Coherence at 6 Hz';
xlabel('Electrodes');
ylabel('Electrodes');
subplot(1,2,2)
imagesc(squeeze(phi(:,:, find(freq >= 6, 1))));
c = colorbar;
c.Label.String = 'Phase at 6 Hz';
xlabel('Electrodes');
ylabel('Electrodes');


%% Estimate delays between electrodes
[delay] = compute_delay(coh, phi, freq, BAND);

%% Estimate parameters of the waves
[src_dir, speed] = estimate_wave(delay, position);

%% Plot the result
plot_wave(delay, position, src_dir, speed);
