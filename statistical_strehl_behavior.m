%%% Goal: show rmse(?) of different strehl ratios w/r/t minimum probe size
%%% over a large number of probes

%1 min per, 6pm -> 10 am (16 hrs = 960 mins)

%% params
n_ab = 1000;
imdim = 1024; simdim=180;

%% Building set of aberrations
abs  = aberration_generator(n_ab);

%% Calculating probe sizes as fn of aperture size
aperture_sizes = 1:1:100;
probe_sizes = zeros(n_ab,length(aperture_sizes));
tic

parfor it = 1:n_ab
    ab = abs(it);
    ab.mag(1) = 0;
    ab.mag = .5*ab.mag;
    ress = [];
    ab_probe_sizes = zeros(1,length(aperture_sizes));
    for jt = 1:length(aperture_sizes)
        aperture_size = aperture_sizes(jt);
        [~,~,~,probe,~] = shifted_ronchigram(ab,[0 0], aperture_size, imdim, simdim);
        lambda = 1.97e-12;
        px_size = lambda/(2*simdim*1e-3);
        res = resolution_test(probe,'probe');
        %title([num2str(aperture_size) 'mrad, ' num2str(res) 'px']);
        %ress(end+1) = res;
        ab_probe_sizes(jt) = res;
    end
    %probe_sizes(it,jt) = res;
    probe_sizes(it,:) = ab_probe_sizes;
    it
end
toc

%% Finding minimum probe size
min_probe_sizes = zeros(1,n_ab);
min_probe_indices = cell(1,n_ab);
for it = 1:n_ab
   min_probe_sizes(it) = min(probe_sizes(it,:)); 
   min_probe_indices{it} = find(probe_sizes(it,:)==min_probe_sizes(it));
end
%% Calculating strehl apertures for aberrations for various strehl cutoffs
[~,Ss] = strehl_calculator(abs, imdim, simdim, .8, 1);
%% Calculating rmses for strehl cutoffs
strehl_cutoffs = 0.1:.01:.95;
num_heur = length(strehl_cutoffs)+2;
probe_size_err = zeros(1,num_heur);
probe_size_err_count = zeros(1,num_heur);
for it = 1:n_ab
    for jt = 1:num_heur
        cutoff_aperture = 0;
        if jt <= length(strehl_cutoffs)
            cutoff_aperture = find(Ss(it,:)>strehl_cutoffs(jt));
            cutoff_aperture = cutoff_aperture(end);
        elseif jt == num_heur - 1 % indiv pi4
            cutoff_aperture = indiv_p4_calculator(abs(it), imdim, simdim);
            
        else %moving window pi4
            cutoff_aperture = pi4_calculator(abs(it), imdim, simdim);
            
            it
        end
        cutoff_aperture = ceil(cutoff_aperture);
        
        if(cutoff_aperture > aperture_sizes(end))
           cutoff_aperture = aperture_sizes(end); 
        end
        
        cutoff_probe_size = probe_sizes(it,cutoff_aperture);
        err = (cutoff_probe_size - min_probe_sizes(it));
        probe_size_err(jt) = probe_size_err(jt) + err;
        probe_size_err_count(jt) = probe_size_err_count(jt) + (err>0);

    end
end
probe_size_err = (probe_size_err./n_ab);
probe_size_err_count = probe_size_err_count/n_ab;

%% Plotting
%figure; plot(aperture_sizes,probe_sizes(1,:)); xlabel('aperture size (mrad)'); ylabel('probe size (px)');
%figure; plot(strehl_cutoffs,probe_size_err,'LineWidth',4); xlabel('Strehl Cutoff'); ylabel('Probe Size Error (px)'); set(gca,'FontSize',20);
%figure; plot(strehl_cutoffs,probe_size_err_count,'LineWidth',4); xlabel('Strehl Cutoff'); ylabel('Probe Size Error Count'); set(gca,'FontSize',20);

figure; 
yyaxis left;
plot([strehl_cutoffs],probe_size_err(1:end-2),'LineWidth',4);  ylabel('Mean Probe Size Error');

line([0 1],[probe_size_err(end-1) probe_size_err(end-1)],'LineStyle','--');
line([0 1],[probe_size_err(end) probe_size_err(end)]);


yyaxis right;
plot([strehl_cutoffs], probe_size_err_count(1:end-2),'LineWidth',4);ylabel('Probe Size Error Count');
set(gca,'FontSize',20);
xlabel('Strehl Cutoff');

line([0 1],[probe_size_err_count(end-1) probe_size_err_count(end-1)],'LineStyle','--');
line([0 1],[probe_size_err_count(end) probe_size_err_count(end)]);
