%Program to loop through different parameters of the lens and save the
%results, make plots etc.
clear variables
tic
addpath('./Functions')
%Set parameters
N_molecules = 1e6; %Number of molecules for each simulation
lens_lengths = 0.534; %lens lengths to be used (m)
lens_radii = 0.02:0.00025:0.022; %lens radii (m)
electrode_radii = 0.0127; %electrode radii (m)
designed_source_to_lens = 0.8052; %The source to lens distance that I've been assuming in the other simulations
source_to_lens_distance =  designed_source_to_lens;
lens_voltage = 30000;
e_no_lens_mu = 2471/1e9; %Efficiency of beamline with no lens
e_no_lens_sigma = 0;


%Open file for saving data
filename = './Data/Loop_over_radius/lens_radius_2.txt';
fileID = fopen(filename,'w');
%Set up headers for data
headers = ["L / m", "R_l / m", "R_e / m", "S to L / m" "Voltage / V", "Efficiency","sigmaE", "Gain","sigmaG", "f_fp", "sigmaFP","n_alive","N_total"];
formatSpec = '';
for i = 1:numel(headers)
    formatSpec = strcat(formatSpec,'%s\t');
end
formatSpec = strcat(formatSpec,'\n');
%Print headers to file
fprintf(fileID,formatSpec,headers);

%Specify format for saving data
formatSpec_data = '';
for i = 1:numel(headers)
    formatSpec_data = strcat(formatSpec_data,'%12.5e\t');
end
formatSpec_data = strcat(formatSpec_data,'\n');

%% Start loop over parameter space
for i = 1:numel(lens_lengths)
    parfor j = 1:numel(lens_radii)
        for k = 1:numel(electrode_radii)
            %Set the lens parameters for this iteration of the loop
            lens_length = lens_lengths(i);
            lens_radius = lens_radii(j);
            electrode_radius = lens_radius;
            
            %Store the results
            [n_alive, n_fp] = electrostatic_lens_function(lens_length, lens_radius, electrode_radius, N_molecules, 2);
            
            %Calculate efficiency, gain and number of molecules hitting
            %field plates (with uncertainties, assuming Poisson dist. for numbers)
            efficiency_mu = n_alive/N_molecules;
            efficiency_sigma = sqrt(n_alive)/N_molecules;
            
            %Calculate gain
            gain_mu = efficiency_mu/e_no_lens_mu;
            gain_sigma = sqrt((efficiency_sigma/e_no_lens_mu)^2 + (efficiency_mu*e_no_lens_sigma/e_no_lens_mu^2)^2);
            
            %Calculate fraction of molecules hitting field plate
            f_fp_mu = n_fp/N_molecules;
            f_fp_sigma = sqrt(n_fp)/N_molecules;
            
            %data = [lens_length lens_radius electrode_radius efficiency_mu efficiency_sigma gain_mu gain_sigma f_fp_mu f_fp_sigma];
            compiled_data(j,:) = [lens_length lens_radius electrode_radius efficiency_mu efficiency_sigma gain_mu gain_sigma f_fp_mu f_fp_sigma];
                      
        end
    end
end

fprintf(fileID,formatSpec_data,compiled_data.');
fclose(fileID);
type(filename)

toc

% fopen(filename,'r')
% fscan(fileID,formatSpec_data)

%% Import data for plotting %%
delimiter = '\t';
header_lines = 1;
data = importdata(filename, delimiter, header_lines);

%The data columns are: {'L / m'  'R_l / m'  'R_e / m'  'Efficiency'  'sigmaE'  'Gain'  'sigmaG'  'f_fp'  'sigmaFP'}
figure

%Plot efficiency vs lens length
subplot(2,2,1)
hold on
title('Efficiency vs lens radius')
errorbar(data.data(:,2), data.data(:,4),data.data(:,5),'-k')
xlabel(data.colheaders(1))
ylabel(data.colheaders(4))
hold off

%Plot gain vs lens length
subplot(2,2,2)
hold on
title('Gain vs lens radius')
errorbar(data.data(:,2), data.data(:,6),data.data(:,7),'x')
xlabel(data.colheaders(1))
ylabel(data.colheaders(6))
hold off