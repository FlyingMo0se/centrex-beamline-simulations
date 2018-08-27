%This function calculates the trajectories of a specified number of TlF
%molecules sent through the given beamline geometry based on the beam
%properties and the geometry of the beamline. The output from the function
%is two structures, one containing the trajectories and the other
%containing the counts for the numbers of molecules that hit each region of
%the beamline

function [trajectories, counts] = trajectories(zone_of_freezing, beamline_geometry, lens, beam, molecule, N_molecules,store_trajectories)

%% Set up parameters for the numerical simulation
rng(1) %Can set seed for rng for use in debugging


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


%Also initialize storage for the trajectories of the detected molecules
trajectories.detected = [];

%Initialize the position and velocity vectors
[x, v] = deal(NaN(3,2*length(beamline_geometry) - 1 + round(lens.L/dz))); %preallocate position, velocity and acceleration arrays


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
    %x_ini = zeros(N_molecules/N_loops,3);s
    
    %Generate initial velocities from Gaussian distributions
    v_ini = ... %Initial velocity from normal distribution
        [random('norm',beam.v_t,beam.sigma_v_x,N_molecules/N_loops,1) ...
        random('norm',beam.v_t,beam.sigma_v_y,N_molecules/N_loops,1) ...
        random('norm',beam.v_z,beam.sigma_v_z,N_molecules/N_loops,1)];
    
    
    
    
    %% Start loop over molecules to get trajectory for each %%
    for i = 1:N_molecules/N_loops
        %% Set up initial velocity and position %%
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
                name = beamline_element.name;
                counts.dead = counts.dead + 1;
                counts.(name) = counts.(name) + 1;
                
                %If the molecule hit the lens, field plates or detection
                %region aperture, store the trajectory
                if strcmp(name,'lens') || strcmp(name,'field_plate') || strcmp(name,'detection_region') && store_trajectories
                    trajectories.(name)(counts.(name),:,:) = [x(:,1:n) NaN(3,length(x)-n)];
                end
                
                break %If molecule is dead, move on to next molecule
            end
        end
        %Test to see if molecule is alive after for loop (this should
        %be redundant)
        if alive == 1
            counts.detected = counts.detected + 1; %Add to count of molecules that made it
            if store_trajectories
                trajectories.detected(counts.detected,:,:) = x;
            end
        end
    end
end