clc; clear;

syms q1 q2 l1 m1 m2 I1 I2 g dq1 dq2 ddq1 ddq2 real 

% Get FK of center of mass
O1 = [0 0 0]';
O2 = [cos(q1) * (q2 + 2 * l1) sin(q1) * (q2 + 2 * l1) 0]';

Oc1 = [0 0 0]';
Oc2 = [cos(q1) * (q2 + l1); sin(q1) * (q2 + l1); 0]

% Velocity jacobians
Jv1 = simplify([diff(Oc1, q1), diff(Oc1, q2)])
Jv2 = simplify([diff(Oc2, q1), diff(Oc2, q2)])

Jw1 = [[0; 0; 1] [0; 0; 0]]
Jw2 = [[0; 0; 0] [0; 0; 0]]

% Calculate kinetic energy

% Rotation matrix for Q1
R1 = [
 cos(q1) -sin(q1) 0
 sin(q1) cos(q1) 0
 0 0 1];

% Second joint is prismatic => no rotation
R2 = [
    0 0 0
    0 0 0
    0 0 0];

% Formula from slide 9
% We assume that we have I1 and I2
D1 = m1 * Jv1' * Jv1 + Jw1' * R1 * I1 * R1' * Jw1;
D2 = m2 * Jv2' * Jv2 + Jw2' * R2 * I2 * R2' * Jw2;

D = D1 + D2;
D = simplify(D)

% Calculate gradient of potentional energy 
% From center of mass
P1 = m1 * g * 0 * sin(q1);
P2 = m2 * g * (sin(q1) * (q2 + l1));
P = P1 + P2;

G1 = diff(P, q1);
G2 = diff(P, q2);
G = [G1; G2];

% Calculate Coriolis
q = [q1; q2]; 

% Velocity of joint 1 and 2
dq = [dq1; dq2];

% Accel
ddq = [ddq1; ddq2];

C = simplify(Coriolis(D, q, dq, 2));

tor = D * ddq + C*dq + G;

D(q1, q2) = subs(D,{m1, m2, I1, I2, l1},{2 1 1 1 2})
C(q1, q2, dq1, dq2) = subs(C * dq ,{m1, m2, I1, I2, l1},{2 1 1 1 2})
G(q1, q2) = subs(G ,{m1, m2, I1, I2, l1, g},{2 1 1 1 2 9.81})

anw = D(q1, q2) * 1 + C(q1, q2, dq1, dq2) * 1 + G(q1, q2)
anw1 = subs(anw, {q1, q2, dq1, dq2, ddq1, ddq2}, {0, 10, 1, 1, 1, 1})

% Init pos
q1_0 = 0;
q2_0 = 10;

% Init speed
dq1_0 = 0;
dq2_0 = 10;

dt = 0.01;

% Forces
U = [0; 20];

n = 100;

for i = 1:n
    q1p(i) = q1_0;
    q2p(i) = q2_0;
    
    dq1p(i) = dq1_0;
    dq2p(i) = dq2_0;
    
    ddq = inv(D(q1_0, q2_0))*(U-C(q1_0, q2_0,dq1_0,dq2_0)-G(q1_0,q2_0));

    dq1_0 = dq1p(i) + double(ddq(1) * dt);
    dq2_0 = dq2p(i) + double(ddq(2) * dt);

    q1_0 = q1p(i) + dq1_0 * dt;
    q2_0 = q2p(i) + dq2_0 * dt;
end

t = 0:0.1:(0.1 * (n-1));

figure
% Plot of q1
plot(t,q1p)

figure
% Plot of q2
plot(t,q2p)