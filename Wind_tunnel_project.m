clear
close all
clc

%% Case 1: Circular cross section

% Test section
R1 = 1;
D1 = 2 * R1;
rho = 1.225;
mu = 1.81e-5;
nu = mu / rho;
TS_ratio = 1.8;         % [1 - 2] To be optimized
L = TS_ratio * D1;
Vt = 50;                % Can change
Re = rho * Vt * D1 / mu;

f_fun = @(f) f - (2*log10(Re*sqrt(f)) - 0.8) .^ -2;
df = @(f) 1 - 1.32547 ./ (f * (log(Re*sqrt(f))) - 0.921034).^3;
f = newton(1, 1e6, 1e-6, f_fun, df);

Kt = f * L / D1;

% Diffuser
theta = 2.5;            % [2° - 3.5°] To be optimized
Aratio = 2.5;           % [2 - 3] To be optimized
A1 = pi * R1^2;
A2 = A1 * Aratio;
R2 = sqrt(A2 / pi);
D2 = 2 * R2;
Vd = Vt / Aratio;

Re_d = Vd * D2 / nu;
f_fun_d = @(f) f - (2*log10(Re_d*sqrt(f)) - 0.8) .^ -2;
df_d = @(f) 1 - 1.32547 ./ (f * (log(Re_d*sqrt(f))) - 0.921034).^3;
f_d = newton(1, 1e6, 1e-6, f_fun_d, df_d);

Kf = (1 - 1 / Aratio^2) * f / (8 * sind(theta));
Ke = -0.00090760^4 - 0.000013310^5 + 0.000013450^6;     % To be changed for Case 2 (square TS)
Kexp = Ke * ((Aratio - 1) / Aratio) ^ 2;

Kd = Kf + Kexp;

% First Corner
Vc = Vd;
CornerWidth = D2;
NVanes = 10;                % To be optimized
GapVanes = CornerWidth / NVanes;
GapToChord = 1/3;
cv = GapVanes / GapToChord;
Rec = Vc * cv / nu;

Kc = 0.1 + 4.55 / (log10(Rec)) ^ 2.58;

% Second leg
V_2l = Vc;

Re_2l = V_2l * D2 / nu;
f_fun_2l = @(f) f - (2*log10(Re_2l*sqrt(f)) - 0.8) .^ -2;
df_2l = @(f) 1 - 1.32547 ./ (f * (log(Re_2l*sqrt(f))) - 0.921034).^3;
f_2l = newton(1, 1e6, 1e-6, f_fun_2l, df_2l);

L_2l = 1.1 * L;

K_2l = f_2l * L_2l / D2;


% Honeycombs
LengthToDiameter = 7;   % [6 - 8] To be optimized
beta = 0.8;             % Porosity
lambda = 0.028;         % To be confirmed

Kh = lambda * (LengthToDiameter + 3) * 1 / beta^2 + (1 / beta - 1) ^ 2;     % Should be around 0.5

% % Safety screens
% wm = 0.1 * D2;          % Can be changed
% sigmas = 0.5;            % [50% - 70%] To be optimized
% dw = 0.05 * wm;         % coefficient [0.05 - 1] To be optimized
% betas = 1 - dw / wm;
% Kmesh = 1;              % [1 - 2] To be optimized
% KRn = 1;              % [2 - 1] by increasing Re To be optimized
% 
% Km = Kmesh * KRn * sigmas + (sigmas / betas)^2;

% Second corner
Kc2 = 0.1 + 4.55 / (log10(Rec)) ^ 2.58;

% Fan leg
theta_fan = 1.5;
Aratio_fan = 1.5;

V_fan_in = V_2l;
V_fan_out = V_fan_in / Aratio_fan;
A3 = Aratio_fan * A2;
D3 = 2 * sqrt(A3 / pi);

Re_fan_in = V_fan_in * D2 / nu;
Re_fan_out = V_fan_out * D3 / nu;
Re_fan = (Re_fan_in + Re_fan_out) / 2;

f_fun_fan = @(f) f - (2*log10(Re_fan*sqrt(f)) - 0.8) .^ -2;
df_fan = @(f) 1 - 1.32547 ./ (f * (log(Re_fan*sqrt(f))) - 0.921034).^3;
f_fan = newton(1, 1e6, 1e-6, f_fun_fan, df_fan);

Kf_fan = (1 - 1 / Aratio_fan^2) * f_fan / (8 * sind(theta_fan));
Kexp_fan = Ke * ((Aratio_fan - 1) / Aratio_fan) ^ 2;

Kfan = Kf_fan + Kexp_fan;

% Third corner
Rec_3 = V_fan_out * cv / nu;

Kc3 = 0.1 + 4.55 / (log10(Rec_3)) ^ 2.58;

% Third leg
V_3l = V_fan_out;

Re_3l = V_3l * D3 / nu;
f_fun_3l = @(f) f - (2*log10(Re_3l*sqrt(f)) - 0.8) .^ -2;
df_3l = @(f) 1 - 1.32547 ./ (f * (log(Re_3l*sqrt(f))) - 0.921034).^3;
f_3l = newton(1, 1e6, 1e-6, f_fun_3l, df_3l);

L_3l = L_2l;

K_3l = f_3l * L_3l / D3;

% Fourth corner
Kc4 = 0.1 + 4.55 / (log10(Rec_3)) ^ 2.58;

% Nozzle
Ln = 2 * L;
Vout_n = Vt;
Vin_n = V_fan_out;

ReIn_n = Vin_n * D3 / nu;
ReOut_n = Vout_n * D1 / nu;
Re_n = (ReIn_n + ReOut_n) / 2;
f_fun_n = @(f) f - (2*log10(Re_n*sqrt(f)) - 0.8) .^ -2;
df_n = @(f) 1 - 1.32547 ./ (f * (log(Re_n*sqrt(f))) - 0.921034).^3;
f_n = newton(1, 1e6, 1e-6, f_fun_n, df_n);

Knt = 0.32 * f_n * Ln / D1;      % Should be 3% of all losses in the tunnel


% Straighteners
t_c_ratio = 0.15;

Ks = 0.045 * t_c_ratio + 0.003;

% Fan
qRatio = 8;                     % [2 - 10] To be optimized
Kfs = (Kt + Kd + Kc + K_2l + Kc2 + Kh + Kfan + Kc3 + K_3l + Kc4 + Knt + Ks);
disp(Knt / Kfs * 100)