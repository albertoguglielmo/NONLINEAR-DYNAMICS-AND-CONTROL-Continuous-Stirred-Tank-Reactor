close all
clear

%% parameters
syms x y u

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

k=K_0;
u=0;
c=C_A0
t=T_in
a=1
m=K_0*H_r/(rho)
g=T_amb;
s=UA/(rho*C_p*V);
R=E_over_R;
% Definisce x come variabile simbolica

g=(exp(-(2*R)/y)* k* R* (2* x* (m *x *(R - y) + k *y^2) + exp(R/y) *(-g *s* x* (R - 2* y) - s *x *(R* (u - y) + y* (-2* u + y)) + a *(R *x *(-t + y) + y *(2 *t* x - c*y + x* y)))))/y^4;

%grid points coordinates
xi_1_vec = 0 : 0.05 : 1;
xi_2_vec =300 : 1 : 400;


%Vdot initialization
Vdot1 = zeros(length(xi_1_vec), length(xi_2_vec));
Vdot2=Vdot1;
Vdot3=Vdot1;
%Vdot copmutation
for i = 1 : length(xi_1_vec)
    for j = 1 : length(xi_2_vec)
        Vdot1(i, j) = subs(g,[x, y,u], [xi_1_vec(i), xi_2_vec(j),-50]);
         % Vdot2(i, j) = subs(g,[x, y,u], [xi_1_vec(i), xi_2_vec(j),0]);
         % Vdot3(i, j) = subs(g,[x, y,u], [xi_1_vec(i), xi_2_vec(j),50]);
         if Vdot1(i, j) < 0.1
             Vdot1(i, j)=-2000;
         end
    end
end


%plotting
[xi_1_mesh, xi_2_mesh] = meshgrid(xi_1_vec, xi_2_vec);
figure
surf(xi_1_mesh, xi_2_mesh, Vdot1')  % Prima superficie
hold on
%surf(xi_1_mesh, xi_2_mesh, Vdot2')
%surf(xi_1_mesh, xi_2_mesh, Vdot3')


set(gca, 'FontSize', 24)
xlabel('Ca', 'interpreter', 'latex')
ylabel('T', 'interpreter', 'latex')
zlabel('Lg(Lf(h)', 'interpreter', 'latex')
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
   (F/V) * (T_in - x(2))- (K_0 *x(1)*  H_r) / (rho * C_p) * exp(-E_over_R / x(2)) - (UA / (rho * V * C_p)) * ( x(2)- (T_amb))+3.377*20];


end

function Vdot = V_dot_fun(xi, xi_eq)



P = [1 0; 0 1];

Vdot = [1,1]*react(xi);

end