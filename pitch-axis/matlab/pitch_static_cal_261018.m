%% Pitch 静态刚度/重力前馈标定（实机，需先采集数据）
% 目的：用少量静态保持点把"等效刚度 k"从 10~18 的不确定区间钉死，
%       并给出 pitch_ff_center_rad（零力矩角）。
%
% 采集方法（Ozone / EmberProbe）：
%   0) 先把 pitch_ff_k 设为 0，否则刚度前馈会把保持力矩抵消掉，拟合出的 k 会偏小；
%   1) manual_target_mode=1，让普通串级角度环工作；
%   2) 依次把 Gimbal.pitch.target_pos 设为当前静止角 ±0.03/±0.06 rad；
%   3) 每个点等 0.5s，记录 INS_angle[2] 与 pitch_torque_log（或 Gimbal.pitch.torque_out）；
%   4) 导出 CSV，至少两列：第 1 列 IMU 俯仰角 (rad)，第 2 列保持力矩 (N·m)。
%
% 用法：把 CSV 路径填到 CSV_PATH；没有 CSV 时可直接在 DATA 里填 [theta tau]。
clear; clc; close all;

CSV_PATH = '';   % 例如 'C:/Users/dfhhds/Desktop/pitch_static_261018.csv'
DATA = [];       % 例：[ -0.23 -0.95; -0.20 -0.55; -0.17 0.00; -0.14 0.48; -0.11 0.95 ]

if ~isempty(CSV_PATH) && isfile(CSV_PATH)
    T = readtable(CSV_PATH, 'VariableNamingRule','preserve','TextType','string');
    names = T.Properties.VariableNames;
    iTh = find(contains(names,'INS_angle'),1);
    iTa = find(contains(names,'pitch_torque_log'),1);
    if isempty(iTa), iTa = find(contains(names,'torque_out'),1); end
    if ~isempty(iTh) && ~isempty(iTa)
        theta = double(T{:,iTh}); tau = double(T{:,iTa});
    else
        M = readmatrix(CSV_PATH);
        theta = M(:,1); tau = M(:,2);
    end
elseif ~isempty(DATA)
    theta = DATA(:,1); tau = DATA(:,2);
else
    error('请填写 CSV_PATH 或 DATA');
end

good = isfinite(theta) & isfinite(tau);
theta = theta(good); tau = tau(good);
if numel(theta) < 3, error('至少需要 3 个静态点'); end

% 线性化：tau = k*(theta - center) + c = k*theta + (c-k*center)
A = [theta, ones(size(theta))];
x = A\tau;
k = x(1);
center = -x(2)/k;
tauFit = A*x;
R2 = 1 - sum((tau-tauFit).^2)/sum((tau-mean(tau)).^2);

fprintf('\n静态标定结果（线性化）：\n');
fprintf('  k = %.3f N·m/rad\n', k);
fprintf('  零力矩角 center = %.4f rad (%.2f deg)\n', center, center*180/pi);
fprintf('  R² = %.4f\n', R2);
fprintf('  建议：pitch_ff_k=%.2f, pitch_ff_center_rad=%.4f\n', k, center);

% 余弦重力模型 tau = A*cos(theta-phi)+C（可选，跨行程更准）
% 用线性初值 + fminsearch
x0 = [abs(k)*0.08, center, mean(tau)];
opt = optimset('Display','off','MaxFunEvals',2e5,'MaxIter',2e5);
xc = fminsearch(@(x) sum((x(1)*cos(theta-x(2))+x(3)-tau).^2), x0, opt);
tauC = xc(1)*cos(theta-xc(2))+xc(3);
R2c = 1 - sum((tau-tauC).^2)/sum((tau-mean(tau)).^2);
fprintf('\n余弦模型 tau = A*cos(theta-phi)+C：\n');
fprintf('  A=%.4f N·m  phi=%.4f rad  C=%.4f N·m  R²=%.4f\n', xc(1),xc(2),xc(3),R2c);

figure('Color','w','Position',[100 100 800 500]);
plot(theta*180/pi, tau, 'ko','MarkerFaceColor','k'); hold on;
thg = linspace(min(theta),max(theta),200)';
plot(thg*180/pi, [thg ones(size(thg))]*x, 'b-','LineWidth',1.5);
plot(thg*180/pi, xc(1)*cos(thg-xc(2))+xc(3), 'r--','LineWidth',1.2);
grid on; xlabel('IMU 俯仰角 (deg)'); ylabel('保持力矩 (N·m)');
legend('静态点','线性 k','余弦重力','Location','best');
title(sprintf('Pitch 静态标定 k=%.2f N·m/rad, center=%.3f rad',k,center));
