function k = spring_constant(lens_r, lens_voltage)

syms E r real

m = 3.7304e-25;       %Molecule mass in kg
h = 6.63e-34;         %Planck constant
h_bar = h/(2*pi);
J = 2;                %Set J
m_J = 0;              %set m_Js
V = lens_voltage;              %Voltage
B = 6.68992e9;    %Rotational constant for TlF in Hz
mu = 4.2282 * 0.393430307 *5.291772e-9/4.135667e-15 * h * 1e-2; %[J/(V/m)]
R = lens_r;          %Radius of lens in meters

W = mu^2*E^2/(2*h*B)*(J*(J+1) - 3*m_J^2)/(J*(J+1)*(2*J-1)*(2*J + 3));

E = 2*V*r/R^2;

W = subs(W);

F = -diff(W,r);

k = abs(double(coeffs(subs(F),r))); %in N/m
