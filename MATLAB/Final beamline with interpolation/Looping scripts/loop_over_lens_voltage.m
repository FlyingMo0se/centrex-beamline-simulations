%Program to loop through different parameters of the lens (length, radius) and save the
%results, make plots etc.
clear variables
tic
addpath('./Functions')
%Set parameters
N_molecules = 1e7; %Number of molecules for each simulation
J = 2; %Set J value
lens_length = .60; %lens lengths to be used (m)
lens_radius = 0.0222; %lens radii (m)
electrode_radius = lens_radius; %electrode radii (m)
designed_source_to_lens = 0.8052; %The source to lens distance that I've been assuming in the other simulations
source_to_lens_distance =  designed_source_to_lens;
lens_voltages = (25:.5:30)*1e3;
e_no_lens_mu = 2471/1e9; %Efficiency of beamline with no lens
e_no_lens_sigma = 0;

%Open file for saving data
filename = 'Data/loop_over_voltage_06_27_2018_J=2_d=1.75in_v=200_StoL=95cm_N=1e7_RK4.txt';
fileID = fopen(filename,'w');
%Set up headers for data
headers = ["L / m", "R_l / m", "R_e / m", "S to L / m" "Voltage / V", "Efficiency","sigmaE", "Gain","sigmaG", "f_fp", "sigmaFP","n_alive","N_total"];
formatSpec = '';
for i = 1:numel(headers)
    formatSpec = strcat(formatSpec,'%s\t');
end
formatSpec = strcat(formatSpec,'\n');

fprintf(fileID,formatSpec,headers);

%Specify format for saving data
formatSpec_data = '';
for i = 1:numel(headers)
    formatSpec_data = strcat(formatSpec_data,'%12.5e\t');
end
formatSpec_data = strcat(formatSpec_data,'\n');
N = {};

%% Start loop over parameter space
parfor i = 1:numel(lens_voltages)
    %Set the lens parameters for this iteration of the loop
    lens_voltage = lens_voltages(i);
    
    %Store the results
    [n_alive, n_fp, numbers_of_hits] = electrostatic_lens_function(lens_length, lens_radius, electrode_radius,source_to_lens_distance,lens_voltage, J, N_molecules, true);
    
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
    compiled_data(i,:) = [lens_length lens_radius electrode_radius source_to_lens_distance lens_voltage efficiency_mu efficiency_sigma gain_mu gain_sigma f_fp_mu f_fp_sigma n_alive N_molecules];
    
    N{i} = numbers_of_hits;
    
end


fprintf(fileID,formatSpec_data,compiled_data.');
fclose(fileID);
type(filename)

toc


%% Import data for plotting %%
delimiter = '\t';
header_lines = 1;
data = importdata(filename, delimiter, header_lines);

%The data columns are: {'L / m'  'R_l / m'  'R_e / m' "S to L / m" 'Efficiency'  'sigmaE'  'Gain'  'sigmaG'  'f_fp'  'sigmaFP'}
figure

%Plot efficiency vs lens length
subplot(2,2,1)
hold on
title('Efficiency vs source to lens distance')
errorbar(data.data(:,5), data.data(:,6),data.data(:,7),'-k')
xlabel(data.colheaders(5))
ylabel(data.colheaders(6))
hold off

%Plot gain vs lens length
subplot(2,2,2)
hold on
title('Gain vs source to lens distance')
errorbar(data.data(:,5), data.data(:,8),data.data(:,9),'x')
xlabel(data.colheaders(5))
ylabel(data.colheaders(8))
hold off


