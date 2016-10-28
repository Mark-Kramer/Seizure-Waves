%Run the Seizing Cortical Model with a chosen source distribution map.
%
%  When the map is set to "fixed_point_source":
%  For 180 s, the source is set to be either on or off, as follows,
%         0- 40 s: the source is *off*
%        41-140 s: the source is *on*
%       141-180 s: the source is *off*
%
%  When the map is "ictal_wavefront":
%  For 200 s, the source is a slowly evolving ictal wavefront, for which,
%         0- 40 s: the source is *off*
%        41-200 s: the source is *on*
%
%CALLS:
%
% . seizing_cortical_field.m
% . make_map.m
%
%OUTPUTS:
%
% . Cell A: The output of function "seizing_cortical_field.m" is saved after
%   each 1 s simulation.  The output is saved in a .mat file 
%   with name "seizing_cortical_field_" MAP_TYPE "k_X.mat" where "X" is an
%   integer specifying the time step, and MAP_TYPE is either "fixed_point_source"
%   or "ictal_wavefront".  The file is saved in the directory OUTPUT_DIR,
%   specified in the code below.
%
% . Cell B: Plots an image of the excitatory population
%   activity every 2 s. To do so, the code reads in the .mat files
%   generated in this cell and saved in directory OUTPUT_DIR.
%
% . Cell C: Select a 10 s interval of simulated ECOG data, and convert it
%   to a format appropriate for wave analysis. The save .mat file can be
%   analyzed in analysis/main_seizure_wave.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Seizing Cortical Model (a modification of the Waikato Cortical Model)
% Copyright (c) 2016 M. A. Kramer
% Department of Mathematics and Statistics, Boston University, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CELL A.  Run the simulation.

clear; close all

%%%% Set the output directory ---------------------------------------------
OUTPUT_DIR = '../data/';  %<--- update this to make sense locally.

%%%% Define the source map. -----------------------------------------------
k=0;                                            % The initial state.
map_type = 'fixed_point_source'; K=140;         % The map is either this one, 
%map_type = 'ictal_wavefront'; K=200;            % Or the map is this one.
[map,state] = make_map(map_type,k,0);           % Then, make the source map.

%%%% Run the simulation. --------------------------------------------------
T0 = 1;                                         %Simulate 1 s intervals.
load('seizing_cortical_field_IC.mat')           %Load the initial conditions to start.

if strcmp(map_type, 'ictal_wavefront')          %When using the ictal wavefront map,     
    last.dVe = zeros(size(last.dVe)) - 1 ;      %... adjust offset to resting potential of excitatory population (mV)
end

for k=1:K                                       %For K intervals of 1 s duration,
    if k<=40                                    %... for first 40 s,
        source_drive=mean(last.dVe(:));         %... include no source drive,
    else                                        %... beyond 40 s,
        source_drive=3;                         %... increase source drive.
    end
    fprintf(['Running simulation ' map_type ', ' num2str(k) ' ... \n'])
                                                %Run the simulation,
    [NP, EC, time, last] = seizing_cortical_field(source_drive, map, T0, last);
                                                %... save the results of this run.
    save([OUTPUT_DIR 'seizing_cortical_field_' map_type '_k_' num2str(k) '.mat'], ...
        'NP','EC','time','last')
    
    if strcmp(map_type, 'ictal_wavefront')      %... when appropriate, update ictal wavefront location.
        [map,state] = make_map(map_type,k,state);
    end
end

if strcmp(map_type, 'fixed_point_source')       %When using the focal point source map,    
                                                %... run simulation for more time, and remove source driver.
    source_drive = 1.5;                         %Set the source drive equivalent to background,
    for k=141:180                               %...  run the simulation for 40 s more,
        fprintf(['Running simulation ' num2str(k) '... \n'])
        [NP, EC, time, last] = seizing_cortical_field(source_drive, map, T0, last);
                                                %... save the results of this run.
        save([OUTPUT_DIR 'seizing_cortical_field_' map_type '_k_' num2str(k) '.mat'], ...
            'NP','EC','time','last')
    end
end

%% Cell B. Visualize the activity of the excitatory population. -----------
%  NOTE:  This cell should be run after Cell A is complete.

clear; close all

%%%% Specifiy directory where simulation results were saved. --------------
OUTPUT_DIR = '../data/';  %<--- update this to make sense locally.

%%%% Choose the map type, and specify total number of indices. ------------
%map_type = 'fixed_point_source'; K=180;         % The map is either this one, 
map_type = 'ictal_wavefront'; K=200;            % Or the map is this one.

%%%% Visualize the excitatory population activity. ------------------------
figure()
counter=1;
for k=1:2:K                                     %For each 1 s interval, in steps of 2 s,
    fprintf(['Read in ' num2str(k) '\n'])       %... print counter,
                                                %... load the saved data,
    load([OUTPUT_DIR 'seizing_cortical_field_' map_type '_k_' num2str(k) '.mat'])
    Q0 = last.Qe;                               %... get the excitatory activity at last moment of time,
    subplot(10,10,counter)                      %... assign a subplot
    imagesc(Q0, [0 25])                         %... and image the excitatory activity with range [0 25] Hz.
    axis off                                    %... exclude axes for visualization.
    counter=counter+1;                          %... augment the counter.
end

%% Cell C. Create example 10 s of ECoG data to analyze wave dynamics.
%  NOTE:  This cell can be run after the Cell A is complete.

clear; close all

%%%% Specifiy directory where simulation results were saved ---------------
OUTPUT_DIR = '../data/';  %<--- update this to make sense locally.

%%%% Choose the map type, and specify a starting index to load. -----------
map_type = 'fixed_point_source';k0=120;         % The map is either this one, 
%map_type = 'ictal_wavefront'; k0=180;          % Or the map is this one.

%%%% Load 10 s of simulated ECoG and create a new data variable -----------
K=10;                                           %# of intervals to load.
T=5000;                                         %Time indices per interval.
N=9;                                            %# electrodes.

ECOG = zeros(K*T, N);                           %Variable to hold ECoG.
for k=1:K                                       %For each time point,
    fprintf(['Read in ' num2str(k0+k) '\n'])    %... load in simulation, 
    load([OUTPUT_DIR 'seizing_cortical_field_' map_type  '_k_'  num2str(k0+k) '.mat'])
    ECOG(1+(k-1)*T:k*T,:)=EC.Qe;                %... and store the result.
end

%%%% Downsample the simulation data from to 500 Hz. -----------------------
dec = 10;
data = zeros(K*T/10, N);
for k=1:N                                       %For each channel,
    d0 = ECOG(:,k);                             %... get the ECoG,
    d0 = decimate(d0,dec);                      %... downsample it,
    data(:,k) = d0;                             %... and store the result.
end

dt = (time(10)-time(9))*dec;                    %Define sampling interval.
fs = 1/dt;                                      %Define sampling frequency.

%%%% Define the positions. ------------------------------------------------
y = [-12, -12, -12,   0,  0,  0,  12, 12, 12];
x = [-12,   0,  12, -12,  0, 12, -12,  0, 12];
position = [x;y]';

%%%% Save the results. ----------------------------------------------------
% NOTE: This file can be loaded and analyzed in analysis/main_seizure_wave.m
save(['../data/example_simulation_waves_' map_type '.mat'], 'data', 'fs', 'position')