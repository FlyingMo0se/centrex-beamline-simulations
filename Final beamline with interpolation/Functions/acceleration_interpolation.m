%Function to calculate accleration at a position within the lens by using
%interpolation

function a = acceleration_interpolation(x,a_values,r_values,g) %inputs are table of accelerations (a_values) at various radii (r_values) and radial position of molecule in lens (r)
a = zeros(3,1);


r = sqrt(x(1).^2 + x(2).^2); %calculate radial position of molecule
%a_r = interp1(r_values,a_values,r); %Calculate radial acceleration by interpolating
a_r = interpolate_custom(r_values,a_values,r); %Calculate radial acceleration by interpolating

if r ~= 0
    %Resolve acceleration into components
    a(1) = a_r.*x(1)./r; %x-direction
    a(2) = a_r.*x(2)./r - g; %y-direction
    a(3) = 0; %z-direction
else
    a = [0 0 0];
end

