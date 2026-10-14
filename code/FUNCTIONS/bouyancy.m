clear all; close all

%% Solving the center of bouyancy/mass equations to determine new drifter length

%f = 0.75:0.01:0.95;
f = 0.95;

% Masses in grams
mdisc = 259.5;
mplug = 112.7;
mtube = 16.27; %16.27*(L^2/2) Is this correct??
mbattery = 146.6;
msensor = 49.6;
mcap = 250;
mant = 195.75;
mballast = 0%10:10:300;

M = mdisc + mballast + mplug + mtube+ mbattery + msensor + mcap + mant;

% Solve for L(f) via A, B and C

A = ((mtube)/2 - (f)/2);

B = mdisc + mplug + (1 - (f/2))*(mbattery + msensor + mcap + mant + mdisc);

C = (mplug*1.25) + (mballast*2.5) - (mbattery*14) - (msensor*10) - (mcap*2) + (mant*15);

L_plus = (-B + sqrt(B.^2 - (4*A*C)))./(2*A);

L_minus = (-B - sqrt(B.^2 - (4*A*C)))./(2*A);


figure;
plot(mballast,L_plus)
hold on
plot(mballast,L_minus)

legend('L_plus','L_minus')

%% TO FLOAT
