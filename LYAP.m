%close all
%clear



%% 1.3

xi_0 = [0; 0];
handle = @(xi) react(xi);
xi_eq = fsolve(handle, xi_0);

%% Lyapunov stability
% grid details
side = 0.5;
step = 0.05;

%grid points coordinates
xi_1_vec = 0 : 0.01 : 1;
xi_1_vec(94)=xi_eq(1);
xi_2_vec = 200 :10 : 340;
xi_2_vec(12) =xi_eq(2);


%Vdot initialization
Vdot = zeros(length(xi_1_vec), length(xi_2_vec));

%Vdot copmutation
for i = 1 : length(xi_1_vec)
    for j = 1 : length(xi_2_vec)
        Vdot(i, j) = V_dot_fun([xi_1_vec(i); xi_2_vec(j)], xi_eq);
    end
end

find(Vdot > 0)

%plotting
[xi_1_mesh, xi_2_mesh] = meshgrid(xi_1_vec, xi_2_vec);
figure
surf(xi_1_mesh, xi_2_mesh, Vdot')
set(gca, 'FontSize', 24)
xlabel('x1', 'interpreter', 'latex')
ylabel('x2', 'interpreter', 'latex')
zlabel('$\dot{V}$', 'interpreter', 'latex')
colorbar

%% functions
function xi_dot = react(x)

% Reactor Parameters
V = 50; % Reactor volume (l)
F= 50; % Inlet volumetric flow rate to the reactor (l/min)
C_A0 = 1; % Feed concentration of component A (mole/l)
K_0 = 7.8e10; % Pre-exponential factor (l/min)
E_over_R = 8567; % Activation energy in the Arrhenius equation (Cal/mole)
R = 1.987; % Universal gas constant (Cal/mole.K), assuming typical value
rho = 900; % Density of the inlet and outlet stream (g/l)
C_p = 0.329; % Heat capacity of inlet and outlet stream (Cal/g.K)
T_in = 350; % Inlet stream temperature (K)
H_r = -5e4; % Heat of reaction (Cal/mole)
UA = 5e4; % Heat transfer term (Cal/min.K)
rho_c = 1000; % Density of the coolant water in the jacket (g/l)
T_amb=293; %Temperatura ambiente (K)

xi_dot = [(F/V) *(C_A0 - x(1)) - K_0 * x(1) * exp(-E_over_R / x(2));...
   (F/V) * (T_in - x(2))- (K_0 *x(1)*  H_r) / (rho * C_p) * exp(-E_over_R / x(2)) - (UA / (rho * V * C_p)) * ( x(2)- (T_amb))];


end

function Vdot = V_dot_fun(xi, xi_eq)



k = 1; 
Q = k * eye(2);

%per la linearizzazione
C_A =   xi_eq(1)   ;         
T =    xi_eq(2)   ;         

% Reactor Parameters
V = 50; % Reactor volume (l)
F= 50; % Inlet volumetric flow rate to the reactor (l/min)
C_A0 = 1; % Feed concentration of component A (mole/l)
K_0 = 7.8e10; % Pre-exponential factor (l/min)
E_over_R = 8567; % Activation energy in the Arrhenius equation (Cal/mole)
R = 1.987; % Universal gas constant (Cal/mole.K), assuming typical value
rho = 900; % Density of the inlet and outlet stream (g/l)
C_p = 0.329; % Heat capacity of inlet and outlet stream (Cal/g.K)
T_in = 350; % Inlet stream temperature (K)
H_r = -5e4; % Heat of reaction (Cal/mole)
UA = 5e4; % Heat transfer term (Cal/min.K)
rho_c = 1000; % Density of the coolant water in the jacket (g/l)
T_amb=293; %Temperatura ambiente (K)


A=[-(F/V)-K_0*exp(-E_over_R / T)  ,- K_0 * C_A * exp(-E_over_R / T)/T^2*(-E_over_R)
       - (K_0 *  H_r) / (rho * C_p) * exp(-E_over_R / T) ,-(F/V) - (K_0 * C_A *  H_r) / (rho * C_p) * exp(-E_over_R / T)/T^2*(-E_over_R)- (UA / (rho * V * C_p))];

% Risolvi l'equazione di Lyapunov
P = lyap(A', Q);


Vdot = (xi-xi_eq)' * P * react(xi);
end

