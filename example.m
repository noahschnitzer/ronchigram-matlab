%% Generate an aberration with random magnitude, angles
aberration = aberration_generator(1);
    % stores aberration function in terms of:
    %   m: order of aberration
    %   n: degree of aberration
    %   angle: angle of aberration in degrees
    %   mag: magnitude of aberration
    %   unit: unit for magnitude of aberration
    % for each aberration (generated up to 5th order by generator).

%% Set up simulation parameters
imdim = 1024;   % Output Ronchigram size in pixels
simdim = 180;   % Simulation RADIUS in reciprocal space in mrad
                % Full simulation extent is 2*simdim, so scaling factor from 
                % mrad to px is (imdim/(2*simdim))
ap_size = 128;  % Objective aperture semi-angle (RADIUS) in mrad
shifts = [0 0]; % Shifts in Ronchigram center in pixels
%% Simulate Ronchigram
ronch = shifted_ronchigram(aberration,shifts,ap_size,imdim,simdim);
%% Display
figure; imagesc(ronch); colormap gray; axis image;
%% Add defocus and recalculate and display
aberration.mag(1) = 100; %in Angstroms
ronch = shifted_ronchigram(aberration,shifts,ap_size,imdim,simdim);
imagesc(ronch); colormap gray; axis image;
%% Calculate and overlay Strehl aperture
S = strehl_calculator(aberration,imdim,simdim,.8,0); %takes a bit
viscircles([imdim/2,imdim/2],S.*imdim/(2.*simdim));