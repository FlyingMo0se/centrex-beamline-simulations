%Note:  Saving trajectories is not working very well with parallelization
%so better to use the non-parallelized version if you want to save the
%trajectories
%This function calculates the trajectories of a specified number of TlF
%molecules sent through the given beamline geometry based on the beam
%properties and the geometry of the beamline. The output from the function
%is two structures, one containing the trajectories and the other
%containing the counts for the numbers of molecules that hit each region of
%the beamline

function [trajectories, counts] = trajectories_parallelized(zone_of_freezing, beamline_geometry, lens, beam, molecule, N_molecules,store_trajectories)

%% Set up parameters for the numerical simulation
%rng(1) %Can set seed for rng for use in debugging


g = 9.81; %take acceleration due to gravity into account

%Step size for RK4 integration inside lens
dz = 0.001; %meters

%Calculate acceleration inside lens as a function of radius from centre
%axis
dr = lens.d_1/2/1e4;
[a_values,r_values] = acceleration_table(lens.V, lens.d_1/2,dr, molecule.J,molecule.m);

%% Set up counters for the numbers of molecules that hit each element and a structure to 
%store trajectories
for i = 1:length(beamline_geometry)
    name = beamline_geometry{i}.name;
    counts.(name) = 0;
    trajectories.(name) = [];
end
%Also set up counters for the total numbers of molecules that are detected
%or hit the apparatus
counts.detected = 0;
counts.dead = 0;

%Set up counters in array form so that they can be used as sliced variables
%within the parallelized loop
counts_array = zeros(length(beamline_geometry),1);
% n_max = 2*(length(beamline_geometry)-1) + 1 + round(lens.L/dz);
%trajectories_detected = zeros(1,3,2*(length(beamline_geometry)-1) + 1 + round(lens.L/dz));
% trajectories_detected = [];
% trajectories_lens = [];
% trajectories_fp = [];
% trajectories_dr = [];

n_alive = 0; %Needed for indexing

%% Start loop to propagate through the beamline
%Using a loop here helps with saving memory, as don't need to store such
%large arrays of initial positions and velocities
N_loops = 100; %number of loops
for n_loop = 1:N_loops
    %Generate initial positions
    theta = random('unif',0,2*pi,N_molecules/N_loops,1); %Initial polar angle for position
    r = (random('unif',0,(zone_of_freezing.d/2)^2,N_molecules/N_loops,1)).^(1/2); %distribution is uniform for r^2
    x_ini = [r.*cos(theta) r.*sin(theta) zeros(N_molecules/N_loops,1)]; % Initial x, y and z position (m)
    
    %Point source initial position
    %x_ini = zeros(N_molecules/N_loops,3);
    
    %Generate initial velocities from Gaussian distributions
    v_ini = ... %Initial velocity from normal distribution
        [random('norm',beam.v_t,beam.sigma_v_x,N_molecules/N_loops,1) ...
        random('norm',beam.v_t,beam.sigma_v_y,N_molecules/N_loops,1) ...
        random('norm',beam.v_z,beam.sigma_v_z,N_molecules/N_loops,1)];

    %% Start loop over molecules to get trajectory for each %%
    parfor i = 1:N_molecules/N_loops
        %% Set up initial velocity and position %%
        [x, v] = deal(zeros(3,2*length(beamline_geometry)- 1 + round(lens.L/dz))); %preallocate position, velocity and acceleration arrays
        x(:,1) = x_ini(i,:); % Initial x, y and z position (m)
        x(3,1) = zone_of_freezing.z; %set initial z coordinate to zone of freezing
        v(:,1) = v_ini(i,:); % Get initial velocity from vector generated earlier
        
        alive = 1; %Vector to describe which molecules are still within allowed region
        n = 1; %set counter to 1 initially
        %%%%Start propagating through beamline %%%%
        %Start loop over the beamline geometry to get all the elements of
        %the beamline one-by-one and propagate the molecules through them
        for j = 1:length(beamline_geometry)
           beamline_element = beamline_geometry{j}; %Get the beamline element
           
            if strcmp(beamline_element.type,'lens') %Test if the element is the lens by comparing strings
                [x,v,alive,n] = propagate_through_lens_interpolation(x,v,a_values,r_values,lens,molecule,n,alive,'circle',g);
            else
                [x,v,alive,n] = propagate_through(x,v,beamline_element,n,alive,beamline_element.type,g);
            end
            %Test if molecule is still alive after passing through the
            %current beamline element:
            if alive == 0
                
                %Craziness to make parallelizing work below
                vector = zeros(length(beamline_geometry));
                vector(j) = 1;
                counts_array = counts_array + vector;
                %Craziness ends
                
                
                %What I'd like to do:
                %counts_array(j) = counts_array(j) + vector;
                %but I haven't managed to get this working
                
%                 %Some if statements below to check if the molecule hit the
%                 %field plates, lens or detection region. If it did, the
%                 %trajectory will be saved:
%                 name = beamline_element.name;
%                 if strcmp(name,'lens') && store_trajectories
%                     trajectories_lens = [trajectories_lens; x];
%                 elseif strcmp(name,'field_plate') && store_trajectories
%                     trajectories_fp = [trajectories_fp; x];
%                 elseif strcmp(name,'detection_region') && store_trajectories
%                     trajectories_dr = [trajectories_dr; x];
%                 end
                
                break %If molecule is dead, move on to next molecule
            end
        end
        %Test to see if molecule is alive after for loop (this should
        %be redundant)
        if alive == 1
            n_alive = n_alive + 1;
            %Store the trajectories in an array every 4th line of the array
            %contains information for a new particle
%             if store_trajectories
%                 trajectories_detected = [trajectories_detected; x];
%             end
        end
    end
end


%% Now extract the data from the arrays into the output structures:
%The counts first:
n_dead = 0;
for i = 1:length(beamline_geometry)
    name = beamline_geometry{i}.name;
    counts.(name) = counts_array(i);
    n_dead = n_dead + counts_array(i);
end
counts.dead = n_dead;
counts.detected = n_alive;

%Trajectories next:
% for i = 1:size(trajectories_detected,1)/3
%     trajectories.detected(i,:,:) = trajectories_detected(3*i-2:3*i,:);
% end
% for i = 1:size(trajectories_lens,1)/3
%     trajectories.lens(i,:,:) = trajectories_lens(3*i-2:3*i,:);
% end
% for i = 1:size(trajectories_fp,1)/3
%     trajectories.field_plate(i,:,:) = trajectories_fp(3*i-2:3*i,:);
% end
% for i = 1:size(trajectories_dr,1)/3
%     trajectories.detection_region(i,:,:) = trajectories_dr(3*i-2:3*i,:);
% end

