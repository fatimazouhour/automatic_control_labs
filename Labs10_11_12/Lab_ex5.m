clear all
close all
s = tf('s');
Gp = 25/(s^3+3.3*s^2+2*s);
Ga = 0.38;
Gs = 2;

% Performance requirements translation

Kd = 4; %requirement 1
Gf = 1/(Kd*Gs);

% Constraint on nu and |Kc| (from S2)

rho_r = 0.15;
R0 = 1
Kp = 12.5;

nu_S2 = 0;
Kc_S2 = Kd^2*R0/(rho_r*Kp*Ga);

% Constraint on nu and |Kc| (from S3)

rho_da = 5.8;
Da0 = 5.5e-3
Kp = 12.5;

nu_S3 = 1;
Kc_S3 = Da0/(rho_da*Ga*Gf*Gs);

% Constraints derived from sinusoidal disturbances and transient
% requirements

wp_plus = 0.02;
ws_minus = 40;

rho_dp = 3.6e-4;
rho_ds = 1.25e-4;

ap = 2*1e-2;
as = 0.1;

M_S_LF_dB = 20*log10(rho_dp/ap);
M_T_HF_dB = 20*log10(rho_ds*Gs/as);

wL = wp_plus*10^(-M_S_LF_dB/40);

wH = ws_minus*10^(M_T_HF_dB/40);

% Transient performance requirements

s_hat_max = 0.09;
tr_max = 2.5;
ts_max = 5;
smorz = abs(log(s_hat_max))/sqrt(pi^2+(log(s_hat_max))^2)
Tp = 1/(2*smorz*sqrt(1-smorz^2));
Sp = (2*smorz*sqrt(2+4*smorz^2+2*sqrt(1+8*smorz^2)))/(sqrt(1+8*smorz^2)+4*smorz^2-1);
wc_tr = (1/tr_max)*1/(sqrt(1-smorz^2))*(pi-acos(smorz))*sqrt(sqrt(1+4*smorz^4)-2*smorz^2);
wc_ts_5 = (1/ts_max)*(-log(0.05)/smorz)*sqrt(sqrt(1+4*smorz^4)-2*smorz^2);

% Summarizing constraints

nu = 1;
Kcin = 0.015; % |Kc|>=0.01 from S3


wc_min = max([wc_tr wc_ts_5 wL*2]);
wc_max = wH/2;

% Let me compute the initial loop function Lin

Lin= Kcin/(s^nu)*Gp*Ga*Gf*Gs;

% Plot Lin in the Nichols chart together with the forbidden region
% (corresponding to constraints on Tp and Sp)

figure(1)
myngridst(Tp,Sp)
nichols(Lin)

figure(2)
bode(Lin)
[num,den]=tfdata(Lin,'v'),
nyquist1(num,den)