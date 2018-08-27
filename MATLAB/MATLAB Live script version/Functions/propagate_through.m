%Function to propagate electrostatic lens simulation through an element
%without time-stepping (so need analytic expression for x and v)
%Pass element into function as a structure with z_1 and z_2 indicating
%start and end position along z, specified diameter or rectangular cross
%section. n indicates steps within simulation (n=1 => initial state).
%Use option to tell if element has square or circular cross section
function [x,v,alive,n] = propagate_through(x,v,element,n,alive,option_1,g)
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
    %% Move forward to start of element
    %Calculate change in z from input position to start of element
    delta_z = element.z_1 - x(3,n);
    
    %Calculate time taken for each molecule to reach start of element
    t_1 = (delta_z)./v(3,n);
    
    %increase n by 1
    n = n + 1;
    
    %Calculate x, y and z positions at start of element
    x(:,n) = ...
        [x(1,n-1) + v(1,n-1).*t_1 ...                %x-position
        x(2,n-1) + v(2,n-1).*t_1 - g*t_1.^(2)/2 ... %y-position
        x(3,n-1) + delta_z];                        %z-position
    
    %Update velocities
    v(:,n) = ...
        [v(1,n-1) ...                          %no change in x-velocity
        v(2,n-1) - g*t_1 ...                   %v_y changed due to gravity
        v(3,n-1)];                             %v_z constant
    
    %Cut molecules outside the allowed region
    if option_1 == 'circle'
        rho = sqrt(x(1,n).^2 + x(2,n).^2); %Calculate radial position
        %molecules stay alive if they are within allowed region and they were previously alive
        alive = alive == 1 & rho < element.d_1/2;
    elseif option_1 == 'square'
        alive = alive == 1 & x(1,n) > element.x_1 & x(1,n) < element.x_2 ...
            & x(2,n) > element.y_1 & x(2,n) < element.y_2;
    elseif option_1 == 'fplate'
        alive = alive == 1 && x(1,n) > element.x_1 && x(1,n) < element.x_2;
    end
    
    %% Propagate through element itself
    %Only do this part if element has finite thickness and molecule didn't hit
    %element earlier
    if element.z_2-element.z_1 > 0 && alive
        %calculate change in z within element
        delta_z = element.z_2 - element.z_1;
        
        %Calculate time taken to propagate through element
        t_2 = delta_z./v(3,n);
        
        %increase n by 1
        n = n + 1;
        
        %Calculate positions at end of element
        x(:,n) = ...
            [x(1,n-1) + v(1,n-1).*t_2 ...                %x-position
            x(2,n-1) + v(2,n-1).*t_2 - g*t_2.^(2)/2 ...  %y-position
            x(3,n-1) + delta_z];                         %z-position
        
        
        %Update velocities
        v(:,n) = ...
            [v(1,n-1) ...                          %no change in x-velocity
            v(2,n-1) - g*t_2 ...                   %v_y changed due to gravity
            v(3,n-1)];
        
        %Cut molecules outside the allowed region
        if option_1 == 'circle'
            rho = sqrt(x(1,n).^2 + x(2,n).^2); %Calculate radial position
            %molecules stay alive if they are within allowed region and they were previously alive
            alive = alive == 1 && rho < element.d_2/2;
            if alive == 0 %calculate position where molecule hits element

                %Calculate how the radius of the element changes with z
                %(for cones)
                if element.d_1 - element.d_2 > 0
                    diff_rho = (element.d_2 - element.d_1)/(2*(element.z_2 - element.z_1));
                else
                    diff_rho = 0;
                end
                       
                %Some messy maths for figuring out how long it takes for
                %the particle to hit the circular or conical element
                t_vec(1) = (diff_rho*v(3,n)*element.d_1 - 2*v(1,n)*x(1,n) - 2*v(2,n)*x(2,n) ... 
                        + sqrt((diff_rho*v(3,n)*element.d_1 - 2*v(1,n)*x(1,n) - 2*v(2,n)*x(2,n))^2 ... 
                        -(4*rho^2-element.d_1^2)*(v(1,n)^2 + v(2,n)^2 - diff_rho^2*v(3,n)^2))) ...
                        /(v(1,n)^2 + v(2,n)^2 - diff_rho^2*v(3,n)^2);
                    
                t_vec(2) = (diff_rho*v(3,n)*element.d_1 - 2*v(1,n)*x(1,n) - 2*v(2,n)*x(2,n) ... 
                        - sqrt((diff_rho*v(3,n)*element.d_1 - 2*v(1,n)*x(1,n) - 2*v(2,n)*x(2,n))^2 ... 
                        -(4*rho^2-element.d_1^2)*(v(1,n)^2 + v(2,n)^2 - diff_rho^2*v(3,n)^2))) ...
                        /(v(1,n)^2 + v(2,n)^2 - diff_rho^2*v(3,n)^2);
                
                %cut away solutions that are not real or are negative
                t_vec(abs(imag(t_vec))>0) = nan;
                t_vec(t_vec<0) = nan;
                
                %smaller solution is the time it took to hit element
                delta_t = real(min(t_vec));
                
                %Update final position to be the point where the molecule
                %hits the element
                x(:,n) = ...
                    [x(1,n-1) + v(1,n-1).*delta_t ...                      %x-position
                    x(2,n-1) + v(2,n-1).*delta_t - g*delta_t.^(2)/2 ...    %y-position
                    x(3,n-1) + delta_t*v(3,n-1)];                          %z-position
            end
        elseif option_1 == 'square'
            alive = alive == 1 && x(1,n) > element.x_1 && x(1,n) < element.x_2 ...
                && x(2,n) > element.y_1 && x(2,n) < element.y_2;
            if alive == 0 %calculate position where molecule hits element
                t_x_minus = (element.x_1-x(1,n-1))/v(1,n-1); %time taken to hit element in -x-direction
                t_x_plus = (element.x_2-x(1,n-1))/v(1,n-1); %time taken to hit element in +x-directio
                
                %For y positions may need to take gravity into account:
                if g == 0
                    t_y_minus = (element.y_1-x(2,n-1))/v(2,n-1); %time taken to hit element in -y-direction
                    t_y_plus = (element.y_2-x(2,n-1))/v(2,n-1); %time taken to hit element in +y-direction
                else
                    %time taken to hit element moving in -y-direction
                    t_y_minus = [v(2,n-1)/g*(1 + sqrt(1-2*(element.y_1-x(2,n-1))*g/v(2,n-1)^2)) ...
                        v(2,n-1)/g*(1 - sqrt(1-2*(element.y_1-x(2,n-1))*g/v(2,n-1)^2))]; %time taken to hit element moving in -y-direction
                    
                    %time taken to hit element moving in +y-direction
                    t_y_plus = [v(2,n-1)/g*(1 + sqrt(1-2*(element.y_2-x(2,n-1))*g/v(2,n-1)^2)) ...
                        v(2,n-1)/g*(1 - sqrt(1-2*(element.y_2-x(2,n-1))*g/v(2,n-1)^2))]; %time taken to hit element moving in +y-direction
                end
                
                
                %Compile all the found times in one vector
                t_vec = [t_x_minus t_x_plus t_y_minus t_y_plus];
                
                
                %Get rid of negative and non-real times
                t_vec(abs(imag(t_vec))>0) = nan;
                t_vec(t_vec<0) = nan;
                
                
                
                delta_t = real(min(t_vec)); %calculate time taken to hit element
                
                %test to check that time is not unreasonably long:
                if delta_t > (element.z_2 - element.z_1)/v(3,n-1)
                    disp("Something wrong in propagate_through.m")
                end
                
                x(:,n) = ...
                    [x(1,n-1) + v(1,n-1).*delta_t ...                       %x-position
                    x(2,n-1) + v(2,n-1).*delta_t - g*delta_t.^(2)/2 ...     %y-position
                    x(3,n-1) + delta_t*v(3,n-1)];                           %z-position
            end
        elseif option_1 == 'fplate' %Separate option for field
            alive = alive == 1 && x(1,n) > element.x_1 && x(1,n) < element.x_2;
            if alive == 0 %calculate position where molecule hits element
                if v(1,n-1) < 0
                    delta_t = (element.x_1-x(1,n-1))/v(1,n-1); %time taken to hit element in -x-direction
                elseif v(1,n-1) > 0
                    delta_t = (element.x_2-x(1,n-1))/v(1,n-1); %time taken to hit element in +x-direction
                end
                
                x(:,n) = ...
                    [x(1,n-1) + v(1,n-1).*delta_t ...                       %x-position
                    x(2,n-1) + v(2,n-1).*delta_t - g*delta_t.^(2)/2 ...     %y-position
                    x(3,n-1) + delta_t*v(3,n-1)];                           %z-position
            end
        else %If none of the options is valid display error message
            disp('Error: Invalid option in propagate_through.m')
        end
    end
end
end