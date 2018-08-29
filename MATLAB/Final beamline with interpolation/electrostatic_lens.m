%Program for simulating trajectories of TlF molecules in CeNTREX from the
%beam source through the electrostatic lens to the detector region.
%Integrating equation of motion for each molecule analytically except for
%the electrostatic lens where RK4 integration is used.
%Initial conditions for molecules are based on
%Taking beamline to consist of: 
%1. Circular aperture at source
%2. Conical collimator close to source
%3. Collimator close to lens
%4. Electrostatic lens
%5. Field plate aperture
%6. Field plates
%7. Rectangular entrance to detection region
%Using an interpolating function to get acceleration from Stark potential

tic
addpath('../Functions')
clear variables

rng(1) %Can set seed for rng for use in debugging

%Natural constants
g = 9.81; %take acceleration due to gravity into account
%g = 0;  %...or don't take g into account (definitely should)
meters_per_inch = 0.0254; %conversion factor from inches to meters

%Molecule parameters
molecule.m = (204.38+19.00)*1.67e-27; %TlF molecule mass in kg
molecule.J = 2; %Rotational quantum number of molecule (m_J = 0 by assumption)

%Simulation parameters
%Step size for RK4 integration inside lens
dz = 0.001; %meters
N_molecules = 1e6; %number of molecules for simulation

%% Start defining the geometry of the beamline
n_elements = 0; %set counter for number of beamline elements
%Beamline parameters (i.e. positions and dimensions of apertures, lenses etc.) (all in meters):
%Measuring z-positions relative to cell aperture

%Cell aperture and zone of freezing dimensions (zone of freezing is the
%region where the molecules "spawn")
cell_aperture.d = meters_per_inch*0.25; %Cell aperture diameter
cell_aperture.z = 0; %Cell aperture position z
zone_of_freezing.z = 1*cell_aperture.d; %position of zone of freezing
zone_of_freezing.d = .02; %diameter of zone of freezing
n_elements = n_elements + 1;

%4K collimator dimensions (cone)
FourK_collimator.z_1 = 0.325*meters_per_inch; %4K collimator start z position
FourK_collimator.d_1 = 0.232*meters_per_inch; %4K collimator start diameter
FourK_collimator.z_2 = FourK_collimator.z_1 + 1*meters_per_inch; %4K collimator end z position
FourK_collimator.d_2 = 1*meters_per_inch; %4K collimator end diameter
n_elements = n_elements + 1;

%4K shield dimensions
FourK_shield.z_1 = FourK_collimator.z_2+ meters_per_inch * 0.375; %4K shield z position
FourK_shield.z_2 = FourK_shield.z_1 + meters_per_inch * 0.25; %4K shield z position
FourK_shield.d_1 = meters_per_inch*1; %4K shield aperture diameter
FourK_shield.d_2 = meters_per_inch*1; %4K shield aperture diameter

%40K shield dimensions
FortyK_shield.z_1 = FourK_shield.z_2 + meters_per_inch * 1.25; %40K shield z position
FortyK_shield.z_2 = FortyK_shield.z_1 + meters_per_inch * .25; %40K shield z position
FortyK_shield.d_1 = meters_per_inch*1; %40K shield aperture diameter
FortyK_shield.d_2 = meters_per_inch*1; %40K shield aperture diameter

%Beam box exit dimensions
BB_exit.z_1 = FortyK_shield.z_2 + meters_per_inch * 2.5; %z position
BB_exit.z_2 = BB_exit.z_1 + meters_per_inch * 0.75; %z position
BB_exit.d_1 = meters_per_inch*4; %diameter of beam box exit
BB_exit.d_2 = BB_exit.d_1; %diameter of beam box exit

%Aperture dimensions (supposed to be 6 mm high and 3 mm wide)
aperture.z_1 = BB_exit.z_2 + 0.03; % 12.5*meters_per_inch; %Take position of aperture to be 3cm from exit of beamsource for now
aperture.z_2 = aperture.z_1; %Thae thickness to be negligible
aperture.h = 0.006; %Height is 6 mm
aperture.w = aperture.h/2; %Width is 3 mm
aperture.x_1 = -aperture.w/2;
aperture.x_2 = aperture.w/2;
aperture.y_1 = -aperture.h/2;
aperture.y_2 = aperture.h/2;

%Lens dimensions
source_to_lens_distance = BB_exit.z_2 + meters_per_inch*25; %Distance from cell aperture to star of lens
lens.L = 0.6; %length of lens (m)
lens.d_1 = 1.75*meters_per_inch; %Lens bore diameter = allowed position within lens before escaping
lens.d_2 = lens.d_1; %Lens bore diameter = allowed position within lens before escaping

%Calculate the spring constant for this lens
k = spring_constant(lens.d_1/2,3e4);

%Code to calculate where lens needs to be to have source at focal length
L = lens.L; %length of lens
m = 3.7304e-25; %Mass of molecules (kg)
v_z = 200; %longitudinal velocity
p = sqrt(k/m)/v_z; %lens constant
d = (1-cos(p*L))./(p*sin(p*L)); %distance of principal plane from start of lens
f = 1/(p*sin(p*L)); %focal length

lens.z_1 = source_to_lens_distance; %position of lens start
%lens.z_1 = zone_of_freezing.z + f - d; %position of lens start if want the source to be at focal point
lens.z_2 = lens.z_1 + lens.L; %Lens end position

lens.dz = 0.001; %step size within lens
lens.V = 30000; %lens voltage
lens.electrode_d = lens.d_1; %electrode diameter
n_elements = n_elements + 1;

%Field plates
aperture_to_field_plate_distance = 0.05; %Distance from field plate aperture to field plates
field_plate.L = 3.000; %length of interaction region
field_plate.z_1 = 1.7727; %position of interaction region start
field_plate.z_2 = field_plate.z_1 + field_plate.L; %position of interaction region end
field_plate.spacing = 0.02; %distance between field plates
field_plate.x_1 = -field_plate.spacing/2; %y-position of lower FP
field_plate.x_2 = +field_plate.spacing/2; %y-position of upper FP
field_plate.y_1 = -5; %not sure what the x-limits should be
field_plate.y_2 = 5; %not sure what the x-limits should be
n_elements = n_elements + 1;


%Detection region
field_plate_to_detection_region_distance = 12.5*meters_per_inch; %Distance from lens to detection region
detection_region.L = 0.05; %length of detection region
detection_region.z_1 = field_plate.z_2 + field_plate_to_detection_region_distance;
detection_region.z_2 = detection_region.z_1 + detection_region.L;
% detection_region.d_1 = 0.010;
% detection_region.d_2 = 0.010;
detection_region.h = 0.03; %height of detection region
detection_region.w = 0.01; %width of detection region
detection_region.x_1 = -detection_region.w/2;
detection_region.x_2 =  detection_region.w/2;
detection_region.y_1 = -detection_region.h/2;
detection_region.y_2 = detection_region.h/2;
n_elements = n_elements + 1;

%% Beam parameters
beam.v_z = 200; %longitudinal velocity (m/s)
beam.sigma_v_z = 13; %spread in longitudinal velocity (m/s)
beam.v_t = 0; %transverse velocities should be zero on average but with some spread
FWHM = pi/4; %Full width half max of distribution of polar angles
%beam.sigma_v_x = tan(FWHM/2)/sqrt(2*log(2)) * beam.v_z; %spread in x velocity (m/s)
beam.sigma_v_x = 57; %spread in x velocity (m/s)
beam.sigma_v_y = beam.sigma_v_x; %spread in y velocity (m/s)
beam.sigma_v_t = 27.6; %spread in transverse velocity (m/s)
beam.d = 0.01; %initial diameter of beam

%Calculate the spring constant for this lens

%set up counters for how many molecules live and die and where
n_alive = 0;
n_dead = 0;
n_enter_lens = 0; %number of molecules that enter lens
n_fp = 0; %number of molecules hitting field plate
n_FourK_collimator = 0; %Number of molecules hitting the conical collimator
n_FourK_shield = 0; %Number of molecules hitting 4K collimator
n_FortyK = 0; %number of molecules hitting 50K collimator
n_lens = 0; %Number of molecules hitting lens
n_fpa = 0; %Number of molecules hitting field plate aperture
n_dr = 0; %Number of molecules that miss the detection region aperture
n_aperture = 0;
n_cell = 0;


disp('Calculating interpolation table for acceleration')

%Calculate acceleration inside lens as a function of radius from centre
%axis
dr = lens.d_1/2/1e4;
[a_values,r_values] = acceleration_table(lens.V, lens.d_1/2,dr, molecule.J,molecule.m);

disp('Interpolating table created, moving on to simulate trajectories')

%% Generate initial velocities and positions (uniform density within zone of
%freezing)

%Using a loop here helps with saving memory, as don't need to store such
%large arrays
N_loops = 100; %number of loops
for n_loop = 1:N_loops
    
    theta = random('unif',0,2*pi,N_molecules/N_loops,1); %Initial polar angle for position
    r = (random('unif',0,(zone_of_freezing.d/2)^2,N_molecules/N_loops,1)).^(1/2); %distribution is uniform for r^2
    x_ini = [r.*cos(theta) r.*sin(theta) zeros(N_molecules/N_loops,1)]; % Initial x, y and z position (m)
    
    %Preallocation for speed
    [x, v] = deal(zeros(3,n_elements + round(lens.L/dz))); %preallocate position, velocity and acceleration arrays
    
    %Point source initial position
    %x_ini = zeros(N_molecules/N_loops,3);
    
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
        
        %% Propagate molecules to 4K collimator %%
        %Use function 'propagate_through' to pass molecules through each
        %beamline element
        [x,v,alive,n] = propagate_through(x,v,FourK_collimator,n,alive,'circle',g);
        
        if alive == 0
            n_dead = n_dead + 1; %Count number of molecules that didn't make it to detection
            n_FourK_collimator = n_FourK_collimator + 1; %Add to number of molecules hitting 4K collimator
            %x_FourK_collimator{n_FourK_collimator} = x(:,1:n); %Store trajectory
            continue %If molecule is dead, start with next molecule
        end
        
        %% Propagate molecules to 4K shield %%
        [x,v,alive,n] = propagate_through(x,v,FourK_shield,n,alive,'circle',g);
        
        if alive == 0
            n_dead = n_dead + 1; %Count number of molecules that didn't make it to detection
            n_FourK_shield = n_FourK_shield+ 1; %Add to number of molecules hitting 4K collimator
            %x_FourK_shield{n_FourK_shield} = x(:,1:n); %Store trajectory
            continue %If molecule is dead, start with next molecule
        end
        
        %% Propagate molecules to 40K shield %%
        [x,v,alive,n] = propagate_through(x,v,FortyK_shield,n,alive,'circle',g);
        
        if alive == 0
            n_dead = n_dead + 1; %Count number of molecules that didn't make it to detection
            n_FortyK = n_FortyK + 1; %Add to number of molecules hitting 4K collimator
            %x_FourK{n_FourK} = x(:,1:n); %Store trajectory
            continue %If molecule is dead, start with next molecule
        end
        %% Propagate molecules to beambox exit %%
        [x,v,alive,n] = propagate_through(x,v,BB_exit,n,alive,'circle',g);
        
        if alive == 0
            n_dead = n_dead + 1; %Count number of molecules that didn't make it to detection
            continue %If molecule is dead, start with next molecule
        end
        %% Propagate to aperture
        [x,v,alive,n] = propagate_through(x,v,aperture,n,alive,'square',g);
        if alive == 0
            n_dead = n_dead + 1; %Count number of molecules that didn't make it to detection
            n_aperture = n_aperture + 1;
            continue %If molecule is dead, start with next molecule
        end
        %% Propagate to lens %%
        [x,v,alive,n] = propagate_through_lens_interpolation(x,v,a_values,r_values,lens,molecule,n,alive,'circle',g);
        %[x,v,alive,n] = propagate_through_lens_interpolation_euler(x,v,a_values,r_values,lens,molecule,n,alive,'circle',g);

        if alive == 0
            n_dead = n_dead + 1; %Count number of molecules that didn't make it to detection
            n_lens = n_lens + 1;
            x_lens{n_lens} = x(:,1:n); %Store trajectory
            continue %If molecule is dead, start with next molecule
        end
        %% Propagate through field plates %%
        [x,v,alive,n] = propagate_through(x,v,field_plate,n,alive,'fplate',g);
        if alive == 0
            n_dead = n_dead + 1; %Count number of molecules that didn't make it to detection
            n_fp = n_fp+1; %number of molecules hitting field plates
            x_fp{n_fp} = x(:,1:n);
            continue %If molecule is dead, start with next molecule
        end
        %% Propagate to detection region %%
        [x,v,alive,n] = propagate_through(x,v,detection_region,n,alive,'square',g);
        
        if alive == 0
            n_dead = n_dead + 1; %Count number of molecules that didn't make it to detection
            n_dr = n_dr + 1;
            x_dr{n_dr} = x(:,1:n);
            continue %If molecule is dead, start with next molecule
        elseif alive == 1
            n_alive = n_alive+1;
            x_alive{n_alive} = x(:,1:n);
            v_alive{n_alive} = v(:,1:n);
        end
    end
end
efficiency = n_alive/N_molecules;

%hold off
%% Plotting below
%Only plot if number of molecules < 1e6
if N_molecules <= 1e6
    %close all
    %Plot trajectories x vs z
    figure
    subplot(2,1,1)
    hold on
    %Plot cell aperture
    plot([0 0],[cell_aperture.d/2 0.05],'k','linewidth',1)
    plot([0 0],[-cell_aperture.d/2 -0.05],'k','linewidth',1)
    %Plot 4K collimator
    plot([FourK_collimator.z_1 FourK_collimator.z_2],[FourK_collimator.d_1/2 FourK_collimator.d_2/2],'k','linewidth',1)
    plot([FourK_collimator.z_1 FourK_collimator.z_2],[-FourK_collimator.d_1/2 -FourK_collimator.d_2/2],'k','linewidth',1)
    %Plot 4K shield
    rectangle('Position',[FourK_shield.z_1 FourK_shield.d_1/2 FourK_shield.z_2-FourK_shield.z_1 1],'FaceColor',[0.85 0.85 0.85])
    rectangle('Position',[FourK_shield.z_1 -FourK_shield.d_1/2-1 FourK_shield.z_2-FourK_shield.z_1 1],'FaceColor',[0.85 0.85 0.85])
    %     plot([FourK_shield.z_1 FourK_shield.z_2],[FourK_shield.d_1/2 0.05],'k','linewidth',4)
    %     plot([FourK_shield.z_1 FourK_shield.z_2],[-FourK_shield.d_1/2 -0.05],'k','linewidth',4)
    %Plot 40K shield
    rectangle('Position',[FortyK_shield.z_1 FortyK_shield.d_1/2 FortyK_shield.z_2-FortyK_shield.z_1 1],'FaceColor',[0.85 0.85 0.85])
    rectangle('Position',[FortyK_shield.z_1 -FortyK_shield.d_1/2-1 FortyK_shield.z_2-FortyK_shield.z_1 1],'FaceColor',[0.85 0.85 0.85])
    %Plot beambox exit
    rectangle('Position',[BB_exit.z_1 BB_exit.d_1/2 BB_exit.z_2-BB_exit.z_1 1],'FaceColor',[0.85 0.85 0.85])
    rectangle('Position',[BB_exit.z_1 -BB_exit.d_1/2-1 BB_exit.z_2-BB_exit.z_1 1],'FaceColor',[0.85 0.85 0.85])
    %Plot lens
    rectangle('Position',[lens.z_1 lens.d_1/2 lens.L 0.02],'facecolor','b')
    rectangle('Position',[lens.z_1 -lens.d_1/2-0.02 lens.L 0.02],'facecolor','b')
    %Plot field plate aperture
    rectangle('Position',[aperture.z_1 aperture.x_2 0.002 0.02],'facecolor','black')
    rectangle('Position',[aperture.z_1 aperture.x_1-0.02 0.002 0.02],'facecolor','black')
    %     plot([aperture.z_1 aperture.z_2],[aperture.d_1/2 0.05],'k','linewidth',1)
    %     plot([aperture.z_1 aperture.z_2],[-aperture.d_1/2 -0.05],'k','linewidth',1)
    %Plot field plates
    rectangle('Position',[field_plate.z_1 field_plate.x_2 field_plate.L 0.02],'facecolor','y')
    rectangle('Position',[field_plate.z_1 field_plate.x_1-0.02 field_plate.L 0.02],'facecolor','y')
    %Plot detection region entrance
    rectangle('Position',[detection_region.z_1 detection_region.x_2 0.02 0.02],'facecolor','black')
    rectangle('Position',[detection_region.z_1 detection_region.x_1-0.02 0.02 0.02],'facecolor','black')
    for i = 1:n_dr
        x_i = x_dr{i};
        plot(x_i(3,:),x_i(1,:),'k-')
    end
    for i = 1:n_lens
        x_i = x_lens{i};
        plot(x_i(3,:),x_i(1,:),'k-')
    end
    for i = 1:n_fpa
        x_i = x_fpa{i};
        plot(x_i(3,:),x_i(1,:),'k-')
    end
    for i = 1:n_fp
        x_i = x_fp{i};
        plot(x_i(3,:),x_i(1,:),'r-')
    end
    for i = 1:n_alive
        x_i = x_alive{i};
        plot(x_i(3,:),x_i(1,:),'g-')
    end
    for i = 1:n_fpa
        x_i = x_fpa{i};
        plot(x_i(3,:),x_i(1,:),'k-')
    end
    
    ylim([-0.06 0.06])
    xlim([0 detection_region.z_2])
    xlabel('z / m')
    ylabel('x / m')
    hold off
    
    %Plot trajectories x vs z
    %figure(2)
    subplot(2,1,2)
    hold on
    %Plot cell aperture
    plot([0 0],[cell_aperture.d/2 0.05],'k','linewidth',1)
    plot([0 0],[-cell_aperture.d/2 -0.05],'k','linewidth',1)
    %Plot 4K collimator
    plot([FourK_collimator.z_1 FourK_collimator.z_2],[FourK_collimator.d_1/2 FourK_collimator.d_2/2],'k','linewidth',1)
    plot([FourK_collimator.z_1 FourK_collimator.z_2],[-FourK_collimator.d_1/2 -FourK_collimator.d_2/2],'k','linewidth',1)
    %Plot 4K shield
    rectangle('Position',[FourK_shield.z_1 FourK_shield.d_1/2 FourK_shield.z_2-FourK_shield.z_1 1],'FaceColor',[0.85 0.85 0.85])
    rectangle('Position',[FourK_shield.z_1 -FourK_shield.d_1/2-1 FourK_shield.z_2-FourK_shield.z_1 1],'FaceColor',[0.85 0.85 0.85])
    %     plot([FourK_shield.z_1 FourK_shield.z_2],[FourK_shield.d_1/2 0.05],'k','linewidth',4)
    %     plot([FourK_shield.z_1 FourK_shield.z_2],[-FourK_shield.d_1/2 -0.05],'k','linewidth',4)
    %Plot 40K shield
    rectangle('Position',[FortyK_shield.z_1 FortyK_shield.d_1/2 FortyK_shield.z_2-FortyK_shield.z_1 1],'FaceColor',[0.85 0.85 0.85])
    rectangle('Position',[FortyK_shield.z_1 -FortyK_shield.d_1/2-1 FortyK_shield.z_2-FortyK_shield.z_1 1],'FaceColor',[0.85 0.85 0.85])
    %Plot beambox exit
    rectangle('Position',[BB_exit.z_1 BB_exit.d_1/2 BB_exit.z_2-BB_exit.z_1 1],'FaceColor',[0.85 0.85 0.85])
    rectangle('Position',[BB_exit.z_1 -BB_exit.d_1/2-1 BB_exit.z_2-BB_exit.z_1 1],'FaceColor',[0.85 0.85 0.85])
    %Plot lens
    rectangle('Position',[lens.z_1 lens.d_1/2 lens.L 0.02],'facecolor','b')
    rectangle('Position',[lens.z_1 -lens.d_1/2-0.02 lens.L 0.02],'facecolor','b')
    %Plot field plate aperture
    rectangle('Position',[aperture.z_1 aperture.y_2 0.002 0.02],'facecolor','black')
    rectangle('Position',[aperture.z_1 aperture.y_1-0.02 0.002 0.02],'facecolor','black')
    %     plot([aperture.z_1 aperture.z_2],[aperture.d_1/2 0.05],'k','linewidth',1)
    %     plot([aperture.z_1 aperture.z_2],[-aperture.d_1/2 -0.05],'k','linewidth',1)
    %Plot field plates
    rectangle('Position',[field_plate.z_1 field_plate.y_2 field_plate.L 0.02],'facecolor','y')
    rectangle('Position',[field_plate.z_1 field_plate.y_1-0.02 field_plate.L 0.02],'facecolor','y')
    %Plot detection region entrance
    rectangle('Position',[detection_region.z_1 detection_region.y_2 0.02 0.02],'facecolor','black')
    rectangle('Position',[detection_region.z_1 detection_region.y_1-0.02 0.02 0.02],'facecolor','black')
    for i = 1:n_dr
        x_i = x_dr{i};
        plot(x_i(3,:),x_i(2,:),'k-')
    end
    for i = 1:n_lens
        x_i = x_lens{i};
        plot(x_i(3,:),x_i(2,:),'k-')
    end
    for i = 1:n_fpa
        x_i = x_fpa{i};
        plot(x_i(3,:),x_i(2,:),'k-')
    end
    for i = 1:n_fp
        x_i = x_fp{i};
        plot(x_i(3,:),x_i(2,:),'r-')
    end
    for i = 1:n_alive
        x_i = x_alive{i};
        plot(x_i(3,:),x_i(2,:),'g-')
    end
    for i = 1:n_fpa
        x_i = x_fpa{i};
        plot(x_i(3,:),x_i(1,:),'k-')
    end
    ylim([-0.06 0.06])
    xlim([0 detection_region.z_2])
    xlabel('z / m')
    ylabel('y / m')
    hold off
    
    %     %Plot initial positions
    %     figure
    %     n = hist3(x_ini(:,[1 2]),'Nbins',[200 200]);
    %     pcolor(n)
    %
    %     %Plot initial divergence angles
    %     alpha = atan(sqrt((v_ini(:,1).^2+v_ini(:,2).^2))./v_ini(:,3));
    %     figure
    %     histogram(alpha)
    %     xlabel("Polar angle")
    %     ylabel('Frequency')
    %
    %     %Check trajectories in 4K collimator
    %     figure
    %     hold on
    %     %Plot cell aperture
    %     plot([0 0],[cell_aperture.d/2 0.05],'k','linewidth',3)
    %     plot([0 0],[-cell_aperture.d/2 -0.05],'k','linewidth',3)
    %     %Plot 4K collimator
    %     plot([FourK_collimator.z_1 FourK_collimator.z_2],[FourK_collimator.d_1/2 FourK_collimator.d_2/2],'k','linewidth',1)
    %     plot([FourK_collimator.z_1 FourK_collimator.z_2],[-FourK_collimator.d_1/2 -FourK_collimator.d_2/2],'k','linewidth',1)
    %     for i = 1:n_FourK_collimator
    %         x_i = x_FourK_collimator{i};
    %         plot(x_i(3,:),x_i(2,:),'b-')
    %     end
    %     for i = 1:n_FourK_shield
    %         x_i = x_FourK_shield{i};
    %         plot(x_i(3,:),x_i(2,:),'r-')
    %     end
    %     xlim([0 0.05])
    %     ylim([-0.05 0.05])
    
    %Plot acceleration vs r
    figure
    plot(r_values,a_values,':.')
    xlabel('Radius / m')
    ylabel('Acceleration / m/s^2')
    title('Accleration of molecules as function of radius within lens')
end

%Write positions and velocities of molecules that make it to detection
%region into a file
save('velocities.mat','v_alive')
save('positions.mat','x_alive')

toc