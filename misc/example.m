%% Generate an aberration with random magnitude, angles
aberration = aberration_generator(1);
    % stores aberration function in terms of:
    %   m: order of aberration
    %   n: degree (rotational symmetry) of aberration
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
%% Simulate Ronchigram and display
ronch = shifted_ronchigram(aberration,shifts,ap_size,imdim,simdim);
figure; imagesc(ronch); colormap gray; axis image;
%% Set defocus and recalculate and display
aberration.mag(1) = 100; %in Angstroms
ronch = shifted_ronchigram(aberration,shifts,ap_size,imdim,simdim);
imagesc(ronch); colormap gray; axis image;
%% Calculate and overlay Strehl aperture
S = strehl_calculator(aberration,imdim,simdim,.8,0); %takes a bit
viscircles([imdim/2,imdim/2],S.*imdim/(2.*simdim),'Color','blue');

%% Calculate and overlay pi/4 total aberration phase shift aperture
p4_ap = pi4_calculator(aberration, imdim, simdim);
viscircles([imdim/2,imdim/2],S.*imdim/(2.*simdim),'Color','yellow');

%% Plot the aberration function phase shift, mask with pi/4 aperture
phase_shift = calculate_aberration_function(aberration,imdim,simdim);
% recommend finding a nice cyclic colormap and plotting phase with that
figure; subplot(121); imagesc(phase_shift); 
subplot(122); imagesc(phase_shift.*aperture_mask(imdim,simdim,p4_ap));
%% Find 50% probe current diameter probe sizes for set of convergence angles and plot
aperture_sizes = 1:25;
probe_sizes = probe_sizer(aberration,imdim,simdim,aperture_sizes); %probe sizes in pixels
probe_sizes = probe_sizes*px_to_ang(simdim); % converting from px to angstroms
figure;
plot(aperture_sizes,probe_sizes);

%% Identify the minimum probe size, visualize the smallest possible probe
% and a more aberrated one (1.5x the optimal CA)
[min_probe_size, min_probe_size_ap] = min(probe_sizes);
probe_opt = calculate_probe(phase_shift, imdim, simdim, min_probe_size_ap, [0,0]);
probe_over = calculate_probe(phase_shift, imdim, simdim, 1.5*min_probe_size_ap, [0,0]);
figure; subplot(121); imagesc(fftshift(probe_opt.*conj(probe_opt)));
subplot(122); imagesc(fftshift(probe_over.*conj(probe_over)));