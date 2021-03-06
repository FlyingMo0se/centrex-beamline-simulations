%Function by Adam West that calculates a polynomial fit to the Stark shift of the J = 2 or 3
%m_j = 0 level of TlF 

J = 2;

if J == 3
    Emax=5.9e6;                         % Maximum E-field (V/m)
elseif J == 2
    Emax = 4e6;
end
Estep=1e4;                              % E-field step (V/m)
Evals=0:Estep:Emax;                     % E-field (V/m)
D=4.23;                                 % Molecule frame dipole moment (Debye)
D=D*3.336e-30;                          % Molecule frame dipole moment (C.m)
Br=0.22315;                             % Rotational constant (cm^-1)
Br=Br*100*3e8*6.63e-34;                 % Rotational constant (J)
kb=1.38e-23;                            % Boltzmann's constant (J/K)
Jmax=10;                                % Maximum J value to consider
msize=(Jmax+1)^2;                       % Size of the resulting Hamiltonian matrix (msize x msize)
mass=224*1.67e-27;                      % Mass (kg)

Js=[];                                  % Define list of J values
Ms=[];                                  % Define list of M values
for i=0:Jmax
    Js=[Js i*ones(2*i+1,1)'];           % Compute list of J values for the Hamiltonian basis
    Ms=[Ms -i:i];                       % Compute list of M values for the Hamiltonian basis
end

V=zeros(msize,(Emax/Estep)+1);          % Preallocate the energy spectrum
for n=1:Emax/Estep+1
    E=(n-1)*Estep;                      % Set the value of the E-field
    H=zeros(msize);                     % Preallocate the Hamiltonian
    for i=1:msize           
        for j=1:msize
            J1=Js(i);                   % Retrieve value of J
            M1=Ms(i);                   % Retrieve value of M
            J2=Js(j);                   % Retrieve value of J'
            M2=Ms(j);                   % Retrieve value of M'
            if J1==J2&&M1==M2           % Diagonal elements given by rotational splitting
                H(i,j)=Br*J1*(J1+1);    
            elseif J1+1==J2&&M1==M2     % Elements for J' = J+1
                H(i,j)=-sqrt(((J1-M1+1)*(J1+M1+1))/((2*J1+1)*(2*J1+3)))*D*E;
            elseif J1-1==J2&&M1==M2     % Elements for J' = J-1
                H(i,j)=-sqrt(((J1-M1)*(J1+M1))/((2*J1-1)*(2*J1+1)))*D*E;
            else                        % All other elements are zero
                H(i,j)=0;
            end
        end
    end
    V(:,n)=sort(real(eig(H)));          % Find the eigenvalues and sort them by size
end

if J == 2
    [V_stark_pc, ~, mu] = polyfit(Evals(1:100),V(9,1:100)-V(9,1),2);
    V_array = V(9,:)-V(9,1);
    figure
    hold on
    plot(Evals,V(9,:)-V(9,1))
    plot(linspace(0,Emax,1000),polyval(V_stark_pc,linspace(0,Emax,1000),[],mu),'--')
    xlabel('Electric field / V/m')
    ylabel('Stark shift / J')
    
    %Scale trapping potential to be in units of transverse velocity:
    v = sqrt(2*V_array/mass);
    
    figure
    plot(Evals,v)
    xlabel('Electric field / V/m')
    ylabel('Transverse velocity / m/s')

    
%     figure
%     plot(Evals,V_array)

elseif J == 3
    [V_stark_pc, ~, mu] = polyfit(Evals(1:100),V(16,1:100)-V(16,1),2);
    V_array = V(16,:)-V(16,1);
    figure
    hold on
    plot(Evals,V(16,:)-V(16,1))
    plot(linspace(0,Emax,1000),polyval(V_stark_pc,linspace(0,Emax,1000),[],mu),'--')
    xlabel('Electric field / V/m')
    ylabel('Stark shift / J')
    
    %Scale trapping potential to be in units of transverse velocity:
    v = sqrt(2*V_array/mass);
    
    figure
    plot(Evals,v)
    xlabel('Electric field / V/m')
    ylabel('Transverse velocity / m/s')
end
