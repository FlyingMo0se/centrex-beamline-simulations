%Function to build a table for an interpolating function for electrostatic lens by
%using calculating the Stark shift at various electric field values and
%then taking gradient
function [a_values,r_values] = acceleration_table(lens_V,lens_R,dr,J,m) %as input take voltage on electrodes, lens radius, step size for radius

r_values = 0:dr:lens_R; %Make an array of radius values for which to calculate the Stark shift   
E_values = 2*lens_V/lens_R^2*r_values; %Convert radius values into electric field values (assuming E = 2*V/R^2 * r within the lens radius)
V_stark = stark_potentials(E_values,J); %Find the Stark shifts for the given electric field values by using a different function
a_values = -gradient(V_stark,dr)/m; %Calculate radial acceleration at each radius based on dV_stark/dr
