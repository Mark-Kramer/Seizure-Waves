%Run the Seizing Cortical Model for 180 s. During this time, the "source"
%  is set to be either on or off, as follows,
%         0- 40 s: the source is *off*
%        41-140 s: the source is *on*
%       141-180 s: the source is *off*

%CALLS:
%
% . seizing_cortical_field.m
%
%OUTPUTS:
%
% . The output of function "seizing_cortical_field.m" is saved after
%   each 1 s simulation.  The output is saved in a .mat file 
%   with name "seizing_cortical_field_k_X.mat" where "X" is an integer 1 to
%   180.  The file is saved in the directory OUTPUT_DIR, specified in the
%   code below.
%
% . The last cell of this file plots an image of the excitatory population
%   activity every 2 s. To do so, the code reads in the .mat files
%   generated in this cell and saved in directory OUTPUT_DIR.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Seizing Cortical Model (a modification of the Waikato Cortical Model)
% Copyright (c) 2016 M. A. Kramer
% Department of Mathematics and Statistics, Boston University, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all

%%%% Set the output directory ---------------------------------------------
OUTPUT_DIR = '~/research/MULTISCALE/dat/simulations/';  %<--- update this to make sense locally.

%%%% Run the simulation. --------------------------------------------------
T0 = 1;                                         %Simulate 1 s intervals.
load('seizing_cortical_field_IC.mat')           %Load the initial conditions to start.

for k=1:140                                     %For 140 intervals of 1 s duration,
    if k<=40                                    %... for first 40 s,
        source_drive=1;                         %... include no source drive,
    else                                        %... beyond 40 s,
        source_drive=3;                         %... increase source drive.
    end
    fprintf(['Running simulation ' num2str(k) '... \n'])
                                                %Run the simulation,
    [NP, EC, time, last] = seizing_cortical_field(source_drive, T0, last);
                                                %... save the results of this run.
    save([OUTPUT_DIR 'seizing_cortical_field_k_' num2str(k) '.mat'], ...
        'NP','EC','time','last')
end

%%%%  Run simulation for more time, and remove source driver.
source_drive = 1.5;                             %Set the source drive equivalent to background, 
for k=141:180                                   %...  run the simulation for 40 s more,
    fprintf(['Running simulation ' num2str(k) '... \n'])
    [NP, EC, time, last] = seizing_cortical_field(source_drive, T0, last);
                                                %... save the results of this run.
    save([OUTPUT_DIR 'seizing_cortical_field_k_' num2str(k) '.mat'], ...
        'NP','EC','time','last')
end

%% Visualize the activity of the excitatory population. -------------------
%  NOTE:  This cell can be run after the cell above is complete.

clear; close all

%%%% Set the output directory ---------------------------------------------
OUTPUT_DIR = '~/research/MULTISCALE/dat/simulations/';  %<--- update this to make sense locally.

%%%% Visualize the excitatory population activity. ------------------------
figure()
counter=1;
for k=1:2:180                                   %For each 1 s interval, in steps of 2 s,
    fprintf(['Read in ' num2str(k) '\n'])       %... print counter,
                                                %... load the saved data,
    load([OUTPUT_DIR 'seizing_cortical_field_k_' num2str(k) '.mat'])
    Q0 = last.Qe;                               %... get the excitatory activity at last moment of time,
    subplot(9,10,counter)                       %... assign a subplot
    colormap(1-gray)                            %... in grayscale, with black = more active
    imagesc(Q0, [0 25])                         %... and image the excitatory activity with range [0 25] Hz.
    axis off                                    %... exclude axes for visualization.
    counter=counter+1;                          %... augment the counter.
end