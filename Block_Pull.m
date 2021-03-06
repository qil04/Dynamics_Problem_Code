% Title: Block Pull
% 
% Purpose:  The purpose of this code is to model a block being pulled by a rope.
% 
% Inputs: Inputs include the Mass of the block, angle of pull, force of pull, coefficients of kinetic and static friction.
% 
% Output: Output includes the acceleration (ax and ay) of the block.
% 
% Notes:  This problem will have different solutions based on the load.  You will have to check for slip and lift.

function Block_Pull()

[ax, ay, str] = vector_acceler_simulate(5, 0.4, 0.3, pi/6, 2);
disp(str);

end

function [ax, ay, str] = vector_acceler_simulate(Pu_force, mu_s, mu_k, angle, mass)

%Lansing Gravity acceleration; Unit system N kg m/s^2; Angle in radian
g = 9.81;

%mass*ax == Pu_force*cos(angle) - friction
%mass*ay == normal_force + Pu_force*sin(angle) - mass*g

%Not slip ax=ay=0
friction_ns = Pu_force*cos(angle);
normal_f_ns = mass*g - Pu_force*sin(angle);

%Slip
syms asx friction_s normal_fs
eqn_1 = mass*asx == Pu_force*cos(angle) - friction_s;
eqn_2 = 0 == normal_fs + Pu_force*sin(angle) - mass*g;
eqn_3 = friction_s == mu_k * normal_fs;
%Another solve system of linear equ without matrix
sol = solve([eqn_1, eqn_2, eqn_3], [asx, friction_s, normal_fs]);
xSol = sol.asx;
ySol = sol.friction_s;
zSol = sol.normal_fs;

if (ySol <= mu_s * zSol)
    str = "Cool! It is slipping!";
    ax = xSol;
    ay = 0;
elseif (friction_ns <= mu_s * normal_f_ns)
    str = "It is not slipping now!";
    ax = 0;
    ay = 0;
end

end
