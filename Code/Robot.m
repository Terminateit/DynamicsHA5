%% Calculate rotation matrixes and Zi
clc;
clear;
close all;
RobotParameters;

T10 = trotz(q1)*transl(0,0, L(1))*trotx(pi/2);
T21 = trotz(q2)*transl(L(2),0, 0)*troty(pi/2);
T32 = transl(0,0,L(3));
T20 = T10*T21;
T30 = T10*T21*T32;

R21 =  T21(1:3, 1:3);
R32 =  T32(1:3, 1:3);

R10 = T10(1:3, 1:3);
R20 =  T20(1:3, 1:3);
R30 =  T30(1:3, 1:3);

Z0 = [0;
      0;
      1];
Z1 = R10(1:3,3);
Z2 = R20(1:3,3);
Z = [Z0, Z1, Z2];
R0 = [R10, R20, R30];
R = [R10, R21, R32];

d_q = [diff(q1,t), diff(q2,t), diff(d3,t)];
dd_q = [diff(q1,t,2), diff(q2,t,2), diff(d3,t,2)];

w = zeros(3,4,'sym');
dw = zeros(3,4,'sym');
ddpc = zeros(3,4,'sym');
ddp = zeros(3,4,'sym');
dp = zeros(3,4,'sym');
dpc = zeros(3,4,'sym');
tau = zeros(1,4,'sym');
centrifugal = zeros(3,4,'sym');
G = zeros(3,4,'sym');
G(1:3,1) = [0 0 -1]'*g;

r1 = T10(1:3,4);
r2 = T21(1:3,4);
r3 = T32(1:3,4);
r = [r1, r2, r3];
%% FORWARD NE 
for i=1:3
    if i ~= 3
        w(1:3,1+i) = R(1:3, 1+(i-1)*3:3+(i-1)*3)'*(w(1:3,i) + ...
           +d_q(i)*Z(1:3,1));
       
        dw(1:3,1+i) = R(1:3, 1+(i-1)*3:3+(i-1)*3)'*(dw(1:3,i) + ...
        + Z(1:3,1)*dd_q(i) + ... 
        + cross(d_q(i)*w(1:3,i),Z(1:3,1)));
    
    ddp(1:3,1+i) = R(1:3, 1+(i-1)*3:3+(i-1)*3)'*ddp(1:3,i) + ... 
        + cross(dw(1:3,1+i), r(1:3,i)) + ...
        + cross(w(1:3,1+i), cross(w(1:3,1+i),r(1:3,i)));
    
    ddpc(1:3,1+i) = ddp(1:3,1+i) + cross(dw(1:3,i+1),-r(1:3,i)/2) + ...
            + cross(w(1:3,1+i), cross(w(1:3,1+i),-r(1:3,i)/2));
    
    dp(1:3,1+i) = R(1:3, 1+(i-1)*3:3+(i-1)*3)'*dp(1:3,i) + cross(w(1:3,1+i), r(1:3,i));
    
    dpc(1:3,1+i) = dp(1:3, i+1) + cross(w(1:3,i+1), -r(1:3,i)/2);
    centrifugal(1:3,i) = cross(w(1:3,1+i), cross(w(1:3,1+i),r(1:3,i)));
    elseif i == 3
        w(1:3,1+i) = R(1:3, 1+(i-1)*3:3+(i-1)*3)'*w(1:3,i);
        
        dw(1:3,1+i) = R(1:3, 1+(i-1)*3:3+(i-1)*3)'*dw(1:3,i);
        
        ddp(1:3,1+i) = R(1:3, 1+(i-1)*3:3+(i-1)*3)'*(ddp(1:3,i) + ... 
        + dd_q(i)*Z(1:3,1)) + ...
        + cross(2*d_q(i)*w(1:3,1+i), R(1:3, 1+(i-1)*3:3+(i-1)*3)'*Z(1:3,1)) + ...
        + cross(dw(1:3,1+i), r(1:3,i)) + ...
        + cross(w(1:3,1+i), cross(w(1:3,1+i),r(1:3,i)));
    
        coriolis = cross(2*d_q(i)*w(1:3,1+i), R(1:3, 1+(i-1)*3:3+(i-1)*3)'*Z(1:3,1));
        centrifugal(1:3,i) = cross(w(1:3,1+i), cross(w(1:3,1+i),r(1:3,i)));
        ddpc(1:3,1+i) = ddp(1:3,1+i) + cross(dw(1:3,i+1),-r(1:3,i)/2) + ...
            + cross(w(1:3,1+i), cross(w(1:3,1+i),-r(1:3,i)/2));
    
        dp(1:3,1+i) = R(1:3, 1+(i-1)*3:3+(i-1)*3)'*(dp(1:3,i) + d_q(i)*Z(1:3,1)) + cross(w(1:3,1+i), r(1:3,i));
        
    
        dpc(1:3,1+i) = dp(1:3, i+1) + cross(w(1:3,i+1), -r(1:3,i)/2);
        
    end
    G(1:3,1+i) = R(1:3, 1+(i-1)*3:3+(i-1)*3)'*G(1:3,i);
end
w = simplify(w)
dw = simplify(dw)
ddp = simplify(ddp)
ddpc = simplify(ddpc)
dp = simplify(dp)
dpc = simplify(dpc)
G = simplify(G)

%% BACKPROPOGATE NE
f = zeros(3,4,'sym');
mu = zeros(3,4,'sym');

for i = 3:-1:1
    if i == 3
        f(1:3,i) = m(i)*ddpc(1:3,i+1) ...
            -m(i)*G(1:3,i+1);
        mu(1:3,i) = -cross(f(1:3,i),r(1:3,i)/2) + I(1:3, 1+(i-1)*3:3+(i-1)*3)*dw(1:3,1+i) + ...
            cross(w(1:3,1+i),I(1:3, 1+(i-1)*3:3+(i-1)*3)*w(1:3,1+i)) ; 
        tau(1,i) = f(1:3,i)'*R(1:3, 1+(i-1)*3:3+(i-1)*3)'*Z(1:3,1);
    else
        f(1:3,i) = R(1:3, 1+(i)*3:3+(i)*3)*f(1:3,1+i) ... 
            + m(i)*ddpc(1:3,i+1) ...
            -m(i)*G(1:3,i+1);
        mu(1:3,i) = -cross(f(1:3,i),r(1:3,i)/2)  + ... 
             + I(1:3, 1+(i-1)*3:3+(i-1)*3)*dw(1:3,1+i) + ...
            R(1:3, 1+(i)*3:3+(i)*3)*mu(1:3,1+i) + ...
            cross(R(1:3, 1+(i)*3:3+(i)*3)*f(1:3,i+1), -r(1:3,i)/2) + ...
            cross(w(1:3,1+i),I(1:3, 1+(i-1)*3:3+(i-1)*3)*w(1:3,1+i));
        tau(1,i) = mu(1:3,i)'*R(1:3, 1+(i-1)*3:3+(i-1)*3)'*Z(1:3,1);
    end
end

mu = simplify(mu)
f = simplify(f)
time = 0:0.01:10;

i = 1;
Plotting(t,time, mu(1,i),mu(2,i),mu(3,i),'$\tau_x$','$\tau_y$', ... 
    '$\tau_z$', 'Torque of 1st Joint', i, 'torque1')
i = i+1;

Plotting(t,time, mu(1,i),mu(2,i),mu(3,i),'$\tau_x$','$\tau_y$', ... 
    '$\tau_z$', 'Torque of 2nd Joint', i, 'torque2')
i = i+1;

Plotting(t,time, mu(1,i),mu(2,i),mu(3,i),'$\tau_x$','$\tau_y$', ... 
    '$\tau_z$', 'Torque of 3rd Joint', i, 'torque3')

i = 1;

Plotting(t,time, f(1,i),f(2,i),f(3,i),'$f_x$','$f_y$', ... 
    '$f_z$', 'Force of 1st Joint', i, 'force1')
i = i +1;

Plotting(t,time, f(1,i),f(2,i),f(3,i),'$f_x$','$f_y$', ... 
    '$f_z$', 'Force of 2nd Joint', i, 'force2')
i = i +1;

Plotting(t,time, f(1,i),f(2,i),f(3,i),'$f_x$','$f_y$', ... 
    '$f_z$', 'Force of 3rd Joint', i, 'force3')

i = 1;

Plotting(t,time, coriolis(1,1),coriolis(2,1),coriolis(3,1),'$Cor_{x}$','$Cor_{y}$', ... 
    '$Cor_{z}$', 'Coriolis Term for 3rd Joint', i, 'Coriolis');

Plotting(t,time, centrifugal(1,i),centrifugal(2,i),centrifugal(3,i),'$Centr_{1x}$','$Centr_{1y}$', ... 
    '$Centr_{1z}$', 'Centrifugal Term for 1st Joint', i, 'Centrifugal1');
i = i + 1;
Plotting(t,time, centrifugal(1,i),centrifugal(2,i),centrifugal(3,i),'$Centr_{2x}$','$Centr_{2y}$', ... 
    '$Centr_{2z}$', 'Centrifugal Term for 2nd Joint', i, 'Centrifugal2');
i = i +1;
Plotting(t,time, centrifugal(1,i),centrifugal(2,i),centrifugal(3,i),'$Centr_{3x}$','$Centr_{3y}$', ... 
    '$Centr_{3z}$', 'Centrifugal Term for 3rd Joint', i, 'Centrifugal3');

i = 2;

Plotting(t,time, G(1,i),G(2,i),G (3,i),'$Grav_{x}$','$Grav_{y}$', ... 
    '$Grav_{z}$', 'Gravity Term for 1st Joint', i, 'Grav1');
i = i + 1;
Plotting(t,time, G(1,i),G(2,i),G(3,i),'$Grav_{x}$','$Grav_{y}$', ... 
    '$Grav_{z}$', 'Gravity Term for 2st Joint', i, 'Grav2');
i = i + 1;
Plotting(t,time, G(1,i),G(2,i),G(3,i),'$Grav_{x}$','$Grav_{y}$', ... 
    '$Grav_{z}$', 'Gravity Term for 3st Joint', i, 'Grav3');


Plotting(t,time, dp(1,4),dp(2,4),dp(3,4),'$\dot{p}_x$','$\dot{p}_y$', ... 
    '$\dot{p}_z$', 'End Effector Velocity of 3rd Joint', i, 'EEVelocity4');

Plotting(t,time, dpc(1,4),dpc(2,4),dpc(3,4),'$\dot{p}_{cx}$','$\dot{p}_{cy}$', ... 
    '$\dot{p}_{cz}$', 'CoM Velocity of 3rd Joint', i, 'CoMVelocity4');

Plotting(t,time, ddp(1,4),ddp(2,4),ddp(3,4),'$\ddot{p}_x$','$\ddot{p}_y$', ... 
    '$\ddot{p}_z$', 'End Effector Acceleration of 3rd Joint', i, 'EEAcceleration4');

Plotting(t,time, ddpc(1,4),ddpc(2,4),ddpc(3,4),'$\ddot{p}_{cx}$','$\ddot{p}_{cy}$', ... 
    '$\ddot{p}_{cz}$', 'CoM Acceleration of 3rd Joint', i, 'CoMAcceleration4');

Plotting(t,time, tau(1,1),tau(1,2),tau(1,3),'$\tau_{1}$','$\tau_{2}$', ... 
    '$\tau_{3}$', 'Torques', i, 'Torques');



