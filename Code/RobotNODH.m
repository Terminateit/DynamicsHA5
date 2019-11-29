%% Calculate rotation matrixes and Zi
clc;
clear all;
RobotParameters;

T10 = trotz(q1)*transl(0,0, L(1));
T21 = trotx(q2)*transl(0,L(2), 0);
T32 = transl(0,L(3),0);
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

d_q = [d_q1, d_q2, d_d3];
dd_q = [dd_q1, dd_q2, dd_d3];

w = zeros(3,4,'sym');
alfa = zeros(3,4,'sym');
alfac = zeros(3,4,'sym');
alfae = zeros(3,4,'sym');

r1 = T10(1:3,4);
r2 = T21(1:3,4);
r3 = T32(1:3,4);
r = [r1, r3, r3];

for i=1:3
    w(1:3,1+i) = R(1:3, 1+(i-1)*3:3+(i-1)*3)*w(1:3,i) + ...
       + R0(1:3, 1+(i-1)*3:3+(i-1)*3)*Z(1:3,i)*d_q(i);
    alfa(1:3,1+i) = R(1:3, 1+(i-1)*3:3+(i-1)*3)*alfa(1:3,i) + ...
        + R0(1:3, 1+(i-1)*3:3+(i-1)*3)*Z(1:3,i)*dd_q(i) + w(1:3,1+i)*d_q(i);
    alfae(1:3,1+i) = R(1:3, 1+(i-1)*3:3+(i-1)*3)*alfae(1:3,i) + ... 
        cross(diff(w(1:3,1+i), t), r(1:3,i)) + ...
        + cross(w(1:3,1+i), cross(w(1:3,1+i),r(1:3,i)));
    alfac(1:3,1+i) = R(1:3, 1+(i-1)*3:3+(i-1)*3)*alfae(1:3,i) + ... 
        cross(diff(w(1:3,1+i), t), r(1:3,i)/2) + ...
        + cross(w(1:3,1+i), cross(w(1:3,1+i),r(1:3,i)/2));
end
w = simplify(w)
alfa = simplify(alfa)
alfae = simplify(alfae)
alfac = simplify(alfac)