%% Pitch 串级 PID 1ms 非线性仿真（2026-10-18）
% 对象（IMU 侧，开环力矩 -> 俯仰角）：
%   J*th'' + b0*th' + tau_c*tanh(th'/w_f) + k*th = tau(t - Td)
% 控制器与固件一致：
%   角度环 P -> 目标角速度 -> 速度环 PI -> 力矩
%   速度环 Ki 写固件时要乘 Ts；前馈用参考角速度的平滑饱和。
%
% 目的：比较现有默认增益、重算增益，以及摩擦前馈开/关，
%       检查 0.07rad（4°）角度阶跃的上升/超调/稳定/峰值力矩。

clear; clc; close all;

here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end

% ---- 辨识模型 ----
P.J    = 0.07008;
P.k    = 15.783;
P.b0   = 0.7784;
P.tauC = 0.9391;
P.Td   = 0.01625;
P.wF   = 0.03;          % 植物摩擦的 tanh 平滑宽度 (rad/s)

% ---- 仿真设置 ----
Ts   = 1e-3;
Tend = 2.5;
t    = (0:Ts:Tend)';
N    = numel(t);

ref = zeros(N,1);
ref(t >= 0.30 & t < 1.50) = 0.07;   % 加长保持段看稳态；实际模式 7 为 0.8s

% 增益族：结构体字段 Kpa, Kpw, Kiw(连续 /s), tauFF, wLin, useFF
G(1) = mk_gain('现有(无弹簧FF)', 11.77, 2.50, 10.0, 1.00, 0.30, true, 0.0);
G(2) = mk_gain('现有/无摩擦FF',  11.77, 2.50, 10.0, 1.00, 0.30, false, 0.0);
G(3) = mk_gain('Kpw2.5 kFF1.0',  11.77, 2.50, 10.0, 0.94, 0.30, true, 1.0*P.k);
G(4) = mk_gain('Kpw1.5 kFF0.8',  11.77, 1.50, 10.0, 0.94, 0.30, true, 0.8*P.k);
G(5) = mk_gain('Kpw1.5 kFF1.0',  11.77, 1.50, 15.0, 0.94, 0.30, true, 1.0*P.k);
G(6) = mk_gain('Kpw1.2 kFF1.0',  11.77, 1.20, 15.0, 0.94, 0.30, true, 1.0*P.k);
G(7) = mk_gain('Kpa15 Kpw1.5',   15.00, 1.50, 15.0, 0.94, 0.30, true, 0.8*P.k);
G(8) = mk_gain('Kpw1.8 kFF0.8',  11.77, 1.80, 15.0, 0.94, 0.30, true, 0.8*P.k);
G(9)  = mk_gain('Kpw1.0 Ki6',    11.77, 1.00,  6.0, 0.94, 0.30, true, 1.0*P.k);
G(10) = mk_gain('Kpw0.8 Ki6',    11.77, 0.80,  6.0, 0.94, 0.30, true, 1.0*P.k);
G(11) = mk_gain('Kpw1.0 Ki8',    11.77, 1.00,  8.0, 0.94, 0.30, true, 1.0*P.k);

res = cell(numel(G),1);
for q = 1:numel(G)
    res{q} = simulate_case(t, ref, P, G(q), Ts);
end

fprintf('\n================ Pitch 0.07rad 角度阶跃（1ms 非线性仿真） ================\n');
fprintf('%-20s  %8s %8s %8s %8s %9s %8s\n', ...
    '方案','上升(ms)','超调(%)','稳定(ms)','峰值Nm','末值误差','饱和tick');
for q = 1:numel(G)
    R = res{q};
    fprintf('%-20s  %8.1f %8.1f %8.1f %8.3f %9.1f%% %8d\n', ...
        G(q).name, 1000*R.tr, 100*R.os, 1000*R.ts, R.tauMax, 100*R.ess, R.satTicks);
end

%% ---- 参数扫描：找"无饱和、超调<5%、稳定<250ms"的折中点 ----
fprintf('\n================ 参数扫描（刚度前馈 + 摩擦前馈） ================\n');
cand = {};
for Kpa = [9 11.77 14]
    for Kpw = [1.0 1.2 1.5]
        for Kiw = [10 15 20]
            for kf = [0.9 1.0]
                g = mk_gain(sprintf('Kpa%.1f/Kpw%.1f/Ki%d/kFF%.1f',Kpa,Kpw,Kiw,kf), ...
                    Kpa,Kpw,Kiw,0.94,0.30,true,kf*P.k);
                R = simulate_case(t, ref, P, g, Ts);
                score = 10*R.satTicks + 100*max(0,R.os-0.05) + ...
                        100*max(0,R.tr-0.15) + 10*max(0,R.ts-0.25) + 10*R.ess;
                cand{end+1} = struct('g',g,'R',R,'score',score); %#ok<SAGROW>
            end
        end
    end
end
sc = cellfun(@(x) x.score, cand);
[~,ord] = sort(sc);
fprintf('%-38s %8s %8s %8s %8s %8s\n','方案','上升ms','超调%','稳定ms','峰值Nm','饱和tick');
for q = 1:min(12,numel(ord))
    x = cand{ord(q)};
    fprintf('%-38s %8.1f %8.1f %8.1f %8.3f %8d\n', ...
        x.g.name, 1000*x.R.tr, 100*x.R.os, 1000*x.R.ts, x.R.tauMax, x.R.satTicks);
end

%% ---- 鲁棒性：控制器 kFF 固定 15.78，植物 k / J / Td 变化 ----
Gchosen = mk_gain('chosen',11.77,1.00,6.0,0.94,0.30,true,P.k);
fprintf('\n================ 鲁棒性检查（控制器 kFF=%.2f 固定） ================\n',P.k);
fprintf('%-18s %8s %8s %8s %8s %8s\n','扰动','上升ms','超调%','稳定ms','峰值Nm','饱和tick');
for kp = [8 10 12 15.783 18 20]
    P2 = P; P2.k = kp;
    R = simulate_case(t, ref, P2, Gchosen, Ts);
    fprintf('plant k=%-10.2f %8.1f %8.1f %8.1f %8.3f %8d\n', ...
        kp,1000*R.tr,100*R.os,1000*R.ts,R.tauMax,R.satTicks);
end
for Jp = [0.06 0.083 0.10]
    P2 = P; P2.J = Jp;
    R = simulate_case(t, ref, P2, Gchosen, Ts);
    fprintf('plant J=%-10.3f %8.1f %8.1f %8.1f %8.3f %8d\n', ...
        Jp,1000*R.tr,100*R.os,1000*R.ts,R.tauMax,R.satTicks);
end
for Tdp = [0.008 0.012 0.01625 0.020]
    P2 = P; P2.Td = Tdp;
    R = simulate_case(t, ref, P2, Gchosen, Ts);
    fprintf('plant Td=%-9.1fms %6.1f %8.1f %8.1f %8.3f %8d\n', ...
        1000*Tdp,1000*R.tr,100*R.os,1000*R.ts,R.tauMax,R.satTicks);
end

%% ---- 图 ----
f = figure('Position',[80 80 1150 760],'Color','w');
subplot(2,2,[1 2]);
hold on; grid on;
for q = [1 4 5 7]
    plot(t, res{q}.theta*180/pi, 'LineWidth',1.5);
end
yline(0.07*180/pi,'k--');
xlabel('时间 (s)'); ylabel('俯仰角 (deg)');
title('角度阶跃：刚度前馈 0 / 0.5k / 1.0k');
legend({G([1 4 5 7]).name},'Location','southeast');

subplot(2,2,3);
hold on; grid on;
for q = [1 4 5 7]
    plot(t, res{q}.tau, 'LineWidth',1.1);
end
yline(3.5,'r--'); yline(-3.5,'r--');
xlabel('时间 (s)'); ylabel('力矩 (N·m)'); title('力矩指令');

subplot(2,2,4);
hold on; grid on;
for q = [2 3 6 8]
    plot(t, res{q}.theta*180/pi, 'LineWidth',1.1);
end
yline(0.07*180/pi,'k--');
xlabel('时间 (s)'); ylabel('俯仰角 (deg)');
title('无摩擦前馈 / 加大积分 的对照');
legend({G([2 3 6 8]).name},'Location','southeast');

saveas(f, fullfile(here,'pitch_closed_loop_261018.png'));
save(fullfile(here,'pitch_closed_loop_261018.mat'),'G','res','P','t','ref');
fprintf('\n图与结果已保存到 %s\n', here);

%% ---- 本地函数 ----
function g = mk_gain(name,Kpa,Kpw,Kiw,tauFF,wLin,useFF,kFF)
    g = struct('name',name,'Kpa',Kpa,'Kpw',Kpw,'Kiw',Kiw, ...
               'tauFF',tauFF,'wLin',wLin,'useFF',useFF,'kFF',kFF);
end

function R = simulate_case(t, ref, P, G, Ts)
    N = numel(t);
    th = 0; w = 0;
    Iout = 0;
    tauCmd = 0;
    Ndelay = max(1, round(P.Td/Ts));
    delay = zeros(Ndelay,1);
    tauApply = 0;

    thLog = zeros(N,1); wLog = zeros(N,1); tauLog = zeros(N,1);
    tvLog = zeros(N,1); sat = 0;

    sub = 10; h = Ts/sub;
    Ki_code = G.Kiw * Ts;

    for n = 1:N
        % 反馈用上一拍的状态（真实系统也是先采样再算）
        err = ref(n) - th;
        tv = G.Kpa * err;
        tv = max(-3.0, min(3.0, tv));       % 角度环输出限幅

        ev = tv - w;
        Iout = Iout + Ki_code * ev;
        Iout = max(-2.0, min(2.0, Iout));   % 速度环积分限幅

        tauRaw = G.Kpw * ev + Iout;
        if G.useFF
            s = tv / G.wLin;
            s = max(-1, min(1, s));
            tauRaw = tauRaw + G.tauFF * s;
        end
        tauRaw = tauRaw + G.kFF * ref(n);   % 模型刚度前馈（按参考角）
        if abs(tauRaw) >= 3.5, sat = sat + 1; end
        tauCmd = max(-3.5, min(3.5, tauRaw));

        % 纯延迟：控制器本拍算出的力矩，Ndelay 拍后作用到对象
        tauApply = delay(1);
        delay = [delay(2:end); tauCmd];

        % 对象用 RK4 子步积分
        x = [th; w];
        for s = 1:sub
            k1 = f_plant(x, tauApply, P);
            k2 = f_plant(x + 0.5*h*k1, tauApply, P);
            k3 = f_plant(x + 0.5*h*k2, tauApply, P);
            k4 = f_plant(x + h*k3, tauApply, P);
            x = x + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
        end
        th = x(1); w = x(2);

        thLog(n)=th; wLog(n)=w; tauLog(n)=tauApply; tvLog(n)=tv;
    end

    % 指标只取 0.3~1.5s 的正阶跃段
    idx = t >= 0.30 & t < 1.50;
    tt = t(idx) - 0.30; yy = thLog(idx);
    target = 0.07;
    i10 = find(yy >= 0.1*target, 1);
    i90 = find(yy >= 0.9*target, 1);
    if isempty(i10) || isempty(i90), tr = NaN; else, tr = tt(i90)-tt(i10); end
    os = max(0, (max(yy)-target)/target);
    band = 0.02*target;
    outside = find(abs(yy-target) > band, 1, 'last');
    if isempty(outside), ts = 0; else, ts = tt(outside); end
    ess = abs(yy(end)-target)/target;

    R = struct('theta',thLog,'w',wLog,'tau',tauLog,'tv',tvLog, ...
        'tr',tr,'os',os,'ts',ts,'tauMax',max(abs(tauLog(idx))), ...
        'ess',ess,'satTicks',sat);
end

function dx = f_plant(x, tau, P)
    th = x(1); w = x(2);
    tauF = P.tauC * tanh(w / P.wF);
    dx = [w; (tau - P.b0*w - tauF - P.k*th) / P.J];
end
