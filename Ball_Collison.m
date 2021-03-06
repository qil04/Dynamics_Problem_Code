function Ball_Collison( )
%angle should be in radian, input the velocity before collison for both A
%and B as magnitude
[Sa, Sb] = Ball_Collison_caluation(5, 5, pi/3, 10, 15, 0.8);
end

function [Sa, Sb]=Ball_Collison_caluation(ma, mb, angle, Va, Vb, e)
%Input mass A, mass B, angle between x and N, coefficient of restitution
%Output speed A and speed B after collison

%Vb along normal direction no tangential velocity
%Conversation of linear momentum in a system
%   ma*Va*cos(angle) + mb*Vb = ma*Vcna + mb*Vcnb
%Conversation of linear momentum of a and b
%   ma*Va*sin(angle) = ma*Vcta
%   Vctb = 0;
%coefficient of resitution
%   e = (Vcna - Vcnb)/(Vb+Va*cos(angle))

%Define the matrix A
%[   V'na,    V'nb,      V'ta,       V'tb;
A = [   1,      -1,        0,          0;...
       ma,      mb,        0,          0;...
        0,       0,       ma,          0;...
        0,       0,        0,         mb];
    
B = [e*(Vb+Va*cos(angle)); -ma*Va*cos(angle) + mb*Vb; ma*Va*sin(angle); 0];

X = A\B;

Sa = norm([X(1),X(3)]);
Sb = norm([X(2),X(4)]);
end

