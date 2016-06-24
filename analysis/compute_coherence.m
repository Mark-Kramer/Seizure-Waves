function [ coh, phi, freq, coh_conf ] = compute_coherence( data, params )
%COMPUTE_COHERENCE Summary of this function goes here
%   Detailed explanation goes here

nb = size(data,1);
nfft = max(2^(nextpow2(nb) + params.pad), nb);
freq = getfgrid(params.Fs, nfft, params.fpass); 

coh = ones([size(data,2) size(data,2) length(freq)]);
phi = zeros([size(data,2) size(data,2) length(freq)]);
coh_conf = zeros([size(data,2) size(data,2) length(freq)]);

for i = 1 : size(data,2)
    d1 = data(:,i);
    d1 = d1 - mean(d1);
    for j = i+1 : size(data,2)
        fprintf('Coherence electrodes %d - %d...\n', i, j);
        d2 = data(:,j);
        d2 = d2 - mean(d2);
        [coh(i,j,:), phi(i,j,:), ~, ~, ~, ~, coh_conf(i,j,:), ~] = coherencyc(d1, d2, params);
        coh(j,i,:) = coh(i,j,:);
        coh_conf(j,i,:) = coh_conf(i,j,:);
        phi(j,i,:) = -phi(i,j,:);
    end
end

end

