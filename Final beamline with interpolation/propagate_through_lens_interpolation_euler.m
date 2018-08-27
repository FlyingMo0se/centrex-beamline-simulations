%Function to propagate electrostatic lens simulation through a lens
%using Runge-Kutta 4th order integration method.
%Pass lens into function as a structure with z_1 and z_2 indicating
%start and end position along z, specified diameter or rectangular cross
%section. n indicates steps within simulation (n=1 => initial state).
%Use option to tell if lens has square or circular cross section
%Acceleration calculated by using interpolation

%Function to propagate electrostatic lens simulation through a lens
%using Runge-Kutta 4th order integration method.
%Pass lens into function as a structure with z_1 and z_2 indicating
%start and end position along z, specified diameter or rectangular cross
%section. n indicates steps within simulation (n=1 => initial state).
%Use option to tell if lens has square or circular cross section

function [x,v,alive,n] = propagate_through_lens_interpolation_euler(x,v,a_values,r_values,lens,molecule,n,alive,option_1,g)
%%Check if molecule is dead or alive
if alive == 0 %If molecule is dead just keep x and v constant
    %Position
    x(:,n+1) = x(:,n);
    x(:,n+2) = x(:,n+1);
    
    %Velocity
    v(:,n+1) = v(:,n);
    v(:,n+2) = v(:,n+1);
    
    %Update n index (up by two)
    n = n+2;
else
    m = molecule.m;
    J = molecule.J;
    if lens.z_2 - lens.z_1 == 0
        return
    end
    %% Move forward to start of lens
    %Calculate change in z from input position to start of lens
    delta_z = lens.z_1 - x(3,n);
    
    %Calculate time taken for each molecule to reach start lens
    t_1 = (delta_z)./v(3,n);
    
    %increase n by 1
    n = n + 1;
    
    %Calculate x and y positions at start of lens
    %Calculate positions at end of element
    x(:,n) = ...
        [x(1,n-1) + v(1,n-1).*t_1 ...                   %x-position
        x(2,n-1) + v(2,n-1).*t_1 - g*t_1.^(2)/2 ...     %y-position
        x(3,n-1) + delta_z];                            %z-position
    
    %Update velocities
    v(:,n) = ...
        [v(1,n-1) ...                           %no change in x-velocity
        v(2,n-1) - g*t_1 ...                   %v_y changed due to gravity
        v(3,n-1)];                             %v_z constant
    
    %Cut molecules outside the allowed region
    if option_1 == 'circle'
        rho = sqrt(x(1,n).^2 + x(2,n).^2); %Calculate radial position
        %molecules stay alive if they are within allowed region and they were previously alive
        alive = alive == 1 & rho < lens.d_1/2;
    elseif option_1 == 'square'
        alive = alive == 1 & x(1,n) > lens.x_1 & x(1,n) < lens.x_2 ...
            & x(2,n) > lens.y_1 & x(2,n) < lens.y_2;
    end
    
    %% Propagate through lens itself using RK4
    %For RK4 have:
    % x(t+dt) = x(t) + dt/6*(k_1+2*k_2+2k_3+k_4)
    % v(t+dt) = v(t) + dt/6*(l_1+2*l_2+2l_3+l_4)
    %Only do this part if lens has finite thickness and molecule didn't hit
    %lens earlier
    if lens.z_2-lens.z_1 > 0 && alive
        N_steps = round((lens.z_2 - lens.z_1)/lens.dz); %calculate number of spatial steps
        dt = lens.dz./v(3,n); %calculate temporal step size
        %Start loop over steps within lens
        for n = n:(n + N_steps)
            if alive == 1 %if molecule is alive update position and velocity
                
                %Calculate x(t+dt) and v(t+dt)
                x(:,n+1) = x(:,n) + dt*v(:,n) + dt^2/2*acceleration_interpolation(x(:,n),a_values,r_values,g);
                v(:,n+1) = v(:,n) + dt*acceleration_interpolation(x(:,n),a_values,r_values,g);
                
                %Cut molecules outside the allowed region
                rho = sqrt(x(1,n+1).^2 + x(2,n+1).^2); %Calculate radial position
                alive = rho < lens.d_2/2 & alive == 1; %molecules stay alive if they are within allowed region and they were previously
            else %if molecule is dead keep position and velocity the same
                x(:,n+1) = x(:,n);
                v(:,n+1) = v(:,n);
            end
        end
    end
end