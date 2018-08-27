%Custom interpolating function for calculating acceleration given that it
%is ordered with increasing radius
function a_r = interpolate_custom(r_values,a_values,r)

%Find the index to which the current radial position of the molecule
%corresponds to
r_max = r_values(length(r_values));
index = floor(r/r_max*length(r_values)) + 1;

%If the index exceeds the size of the interpolation tabels, the particle
%must be very close to the edge of the lens so set the  index so that the
%acceleration outputted is the acceleration at r = r_max
if index > length(r_values)
    index = length(r_values);
end

%Calculate acceleration by interpolating linearly:
if index < length(r_values)
    a_r = a_values(index) + (a_values(index+1) - a_values(index))/(r_values(index+1)-r_values(index))...
        *(r - r_values(index));
elseif index == length(r_values)
    a_r = a_values(index);
else
    disp('Error in interpolate_custom.m')
end
    