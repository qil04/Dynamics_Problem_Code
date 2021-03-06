% Title: Refrigerator push
% 
% Purpose:  The purpose of this code is to model a refrigerator being pushed.
% 
% Inputs: Inputs include the Mass (M), dimensions (B,H), Height (h) Angle(theta) and magnitude of the push (P), coefficients of kinetic and static friction.
% 
% Output: Output includes the acceleration (ax, ay, alpha) and Forces (NA, NB, FfA, FfB).
% 
% Notes:  This problem will have different solutions based on the inputs.  You will have to check for slip, tip, and lift.
% 
function Refrigerator_push()
[ax, ay, alpha, NA, NB, FfA, FfB] = ref_notslip(100, 0.6, 1.8, 1.5, 200, 0.5, 0.8)
end

function [ax, ay, alpha, NA, NB, FfA, FfB] = ref_notslip(m, B, H, h, P, mu_k, mu_s)
%Variable Clarification: B, H as refrigerator's 2-D dimension
%               h as the ground to forcing point force higher than centroid
%               mass in unit kg
%               Assume the refrigerator is geometrically symmetric

%No Slip: x1-friction; y1-Na; z1-Nb 
I = (1/12)*m*(B^2+H^2);
g = 9.81;
syms x1 y1 z1
equ_1x = P - x1 == 0;
equ_1y = y1 + z1 - m*g == 0;
equ_1m = -(H/2)*(x1) - (B/2)*y1 + (B/2)*z1 - (h-H/2)*P == 0;
sol_ns = solve([equ_1x, equ_1y, equ_1m], [x1, y1, z1]);

Ff = sol_ns.x1;
NA = sol_ns.y1;
NB = sol_ns.z1;

%Check condition
% NA>0, NB>0, Ff<=mu_s*(NA+NB)
if Ff <= mu_s*(NA+NB) && NA > 0 && NB >0
    disp("It doesn't slip!")
    ax = 0; ay=0; alpha=0;
    FfA = Ff/2; FfB = Ff/2;
elseif NA <= 0
    syms fb2 Nb2 ax2 ay2 apha2
    equ_2x = P - fb2 == m*ax2;
    equ_2y = Nb2 - m*g == m*ay2;
    equ_2m = -(H/2)*(fb2) + (B/2)*Nb2 - (h-H/2)*P == I*apha2;
    equ_2i = ax2 == (-H/2)*apha2;
    equ_2j = ay2 == (-B/2)*apha2;
    sol_tipB = solve([equ_2x, equ_2y, equ_2m, equ_2i, equ_2j], [fb2, Nb2, ax2, ay2, apha2]);
    FfB = sol_tipB.fb2; FfA = 0;
    NB = sol_tipB.Nb2;
    ax = sol_tipB.ax2;
    ay = sol_tipB.ay2;
    alpha = sol_tipB.apha2;
    disp("Oops, B tips!")
elseif NB <=0
    syms fa2 Na2 ax3 ay3 alpha3
    equ_3x = P - fa2 == m*ax3;
    equ_3y = Na2 - m*g == m*ay3;
    equ_3m = -(H/2)*(fa2) - (B/2)*Na2 - (h-H/2)*P == I*alpha3;
    equ_3i = ax3 == (-H/2)*alpha3;
    equ_3j = ay3 == (B/2)*alpha3;
    sol_tipA = solve([equ_3x, equ_3y, equ_3m, equ_3i, equ_3j], [fa2, Na2, ax3, ay3, alpha3]);
    FfB = 0; FfA = sol_tipA.fa2;
    NA = sol_tipA.Na2;
    ax = sol_tipA.ax3;
    ay = sol_tipA.ay3;
    alpha = sol_tipA.alpha3;
    disp("Oops, A tips!")
elseif Ff > mu_s*(NA+NB)
    syms s_fa s_fb s_ax s_Na s_Nb
    equ_4x = P - s_fa - s_fb == s_ax;
    equ_4y = s_Na + s_Nb - m*g == 0;
    equ_4m = -(H/2)*(s_fa+s_fb) - (B/2)*s_Na + (B/2)*s_Nb - (h-H/2)*P == 0;
    equ_fa = mu_k*s_Na;
    equ_fb = mu_k*s_Nb;
    sol_slip = solve([equ_4x, equ_4y, equ_4m, equ_fa, equ_fb], [s_fa, s_fb, s_ax, s_Na, s_Nb]);
    disp("Look, it slips!")
    ax = sol_slip.s_ax; ay = 0; alpha = 0;
    FfA = sol_slip.s_fa; FfB = sol_slip.s_fb;
    NA = sol_slip.s_Na; NB = sol_slip.s_Nb;
end

end