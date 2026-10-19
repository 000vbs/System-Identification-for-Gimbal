%% 自动收集并联合拟合所有独立 Pitch 角度阶跃 CSV
% 识别方式：以实测力矩作为植物输入，联合拟合 J、b0、tau_c 和两个工作区的局部刚度。
% 同一段阶跃若因采样变量不同而出现在多个 CSV 中，会自动去重。
% 多幅值回中段因静摩擦死区不参与辨识，避免把零速保持误差算进动态参数。
clear; clc; close all;

here = fileparts(mfilename('fullpath'));
pitchRoot = fileparts(here);
OUT_DIR = fullfile(pitchRoot,'generated');
if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end
ROOTS = {fullfile(pitchRoot,'data')};

%% ---------- 1. 收集 CSV 中的独立阶跃段 ----------
segs = {};
seen = containers.Map('KeyType','char','ValueType','logical');
duplicates = {};

for r = 1:numel(ROOTS)
    files = dir(fullfile(ROOTS{r},'emberprobe-live-*.csv'));
    for f = 1:numel(files)
        csvPath = fullfile(files(f).folder,files(f).name);
        candidates = extract_segments(csvPath);
        for c = 1:numel(candidates)
            S = candidates{c};
            key = sprintf('%s|%.6f|%d',S.startTime,S.nominalAmp,S.n);
            if isKey(seen,key)
                duplicates{end+1} = sprintf('%s (%s amp=%.4f)',S.source,S.kind,S.nominalAmp); %#ok<SAGROW>
                continue;
            end
            seen(key) = true;
            segs{end+1} = S; %#ok<SAGROW>
        end
    end
end

if isempty(segs)
    error('没有找到包含角度阶跃状态列的有效 CSV。');
end

fprintf('\n找到 %d 个独立阶跃段：\n',numel(segs));
for s = 1:numel(segs)
    S = segs{s};
    fprintf('  %d) %s | %s | amp=%.4f rad | n=%d | median dt=%.1f ms\n', ...
        s,S.source,S.kind,S.nominalAmp,S.n,1000*S.dtMed);
end
if ~isempty(duplicates)
    fprintf('去重 %d 个重复段：\n',numel(duplicates));
    for i = 1:numel(duplicates)
        fprintf('  %s\n',duplicates{i});
    end
end

%% ---------- 2. 统一到 1 ms 并联合拟合植物 ----------
Ts = 1e-3;
Td = 0.01625;
for s = 1:numel(segs)
    segs{s} = prepare_segment(segs{s},Ts);
end

x0s = [ ...
    0.044 0.00 2.44  9.2 39.2 0.005;
    0.044 0.50 2.40  9.0 33.0 0.020;
    0.020 0.50 2.00 15.0 40.0 0.030;
    0.070 1.50 0.10 50.0 50.0 0.010];
opt = optimset('Display','off','MaxFunEvals',2500,'MaxIter',2500, ...
    'TolFun',1e-10,'TolX',1e-9);

best = []; bestCost = inf; bestExit = [];
for i = 1:size(x0s,1)
    [xb,~,ef] = fminsearch(@(x) joint_cost(x,segs,Ts,Td),x0s(i,:),opt);
    cb = joint_cost(xb,segs,Ts,Td);
    fprintf('初值 %d: cost=%.5f -> J=%.5f b0=%.4f tauC=%.4f k=%.3f\n', ...
        i,cb,xb(1),xb(2),xb(3),xb(4));
    if cb < bestCost
        best = xb; bestCost = cb; bestExit = ef;
    end
end

J = best(1); b0 = best(2); tauC = best(3); kLow = best(4); kHigh = best(5); tauAct = best(6);
fprintf('\n联合拟合结果: J=%.6f, b0=%.5f, tau_c=%.5f, kLow=%.4f, kHigh=%.4f, tau_act=%.5f s, cost=%.6f, exit=%d\n', ...
    J,b0,tauC,kLow,kHigh,tauAct,bestCost,bestExit);

%% ---------- 3. 逐段误差和指标 ----------
nSeg = numel(segs);
source = strings(nSeg,1); kind = strings(nSeg,1);
nominalAmp = zeros(nSeg,1); nrmsePct = zeros(nSeg,1); r2 = zeros(nSeg,1);
peakErr = zeros(nSeg,1); peakVelMeas = zeros(nSeg,1); peakVelSim = zeros(nSeg,1);
riseMeasMs = nan(nSeg,1); riseSimMs = nan(nSeg,1);
overshootMeasPct = nan(nSeg,1); overshootSimPct = nan(nSeg,1);

for s = 1:nSeg
    S = segs{s};
    [ys,ws] = sim_segment(S,best,Ts,Td);
    ey = ys-S.ym;
    eScale = max(max(S.ym)-min(S.ym),0.02);
    nrmsePct(s) = 100*sqrt(mean(ey.^2))/eScale;
    r2(s) = 1-sum(ey.^2)/sum((S.ym-mean(S.ym)).^2);
    peakErr(s) = max(abs(ey));
    peakVelMeas(s) = max(abs(S.wm));
    peakVelSim(s) = max(abs(ws));
    source(s) = string(S.source);
    kind(s) = string(S.kind);
    nominalAmp(s) = S.nominalAmp;

    % 只对真实阶跃段计算上升和超调；斜坡段仍保留 NRMSE/R²。
    if S.rampMs <= 1
        tm = step_metrics(S.tss,S.ym-S.y0,S.nominalAmp);
        ts = step_metrics(S.tss,ys-S.y0,S.nominalAmp);
        riseMeasMs(s) = 1000*tm.rise; riseSimMs(s) = 1000*ts.rise;
        overshootMeasPct(s) = 100*tm.overshoot;
        overshootSimPct(s) = 100*ts.overshoot;
    end

    fprintf('段 %d %s (%s): NRMSE=%.2f%%, R2=%.3f, 峰值误差=%.4f rad, 峰值速度 实测/仿真=%.2f/%.2f\n', ...
        s,S.source,S.kind,nrmsePct(s),r2(s),peakErr(s),peakVelMeas(s),peakVelSim(s));
end

metrics = table(source,kind,nominalAmp,nrmsePct,r2,peakErr, ...
    peakVelMeas,peakVelSim,riseMeasMs,riseSimMs,overshootMeasPct,overshootSimPct);
metricsPath = fullfile(OUT_DIR,'pitch_step_all_csv_fit_metrics_261019.csv');
writetable(metrics,metricsPath);

save(fullfile(OUT_DIR,'pitch_step_all_csv_joint_fit_261019.mat'), ...
    'best','J','b0','tauC','kLow','kHigh','tauAct','Td','Ts','metrics','segs','duplicates');

%% ---------- 4. 实测与仿真总图 ----------
fig = figure('Color','w','Visible','off','Position',[30 30 1580 1180]);
tiledlayout(nSeg,3,'TileSpacing','compact','Padding','compact');

for s = 1:nSeg
    S = segs{s};
    [ys,ws] = sim_segment(S,best,Ts,Td);

    nexttile;
    plot(S.tss,S.ym*180/pi,'Color',[0.05 0.35 0.85],'LineWidth',1.5); hold on;
    plot(S.tss,ys*180/pi,'--','Color',[0.85 0.25 0.15],'LineWidth',2.0);
    plot(S.tss,S.refm*180/pi,':','Color',[0.15 0.15 0.15],'LineWidth',1.0);
    grid on; ylabel('角度 (deg)');
    title(sprintf('%s | %s | amp %.3f rad',S.source,S.kind,S.nominalAmp),'Interpreter','none');
    if s == 1
        legend('实测','联合模型','参考','Location','southeast');
    end

    nexttile;
    plot(S.tss,S.wm,'Color',[0.05 0.45 0.25],'LineWidth',1.4); hold on;
    plot(S.tss,ws,'--','Color',[0.85 0.25 0.15],'LineWidth',1.8);
    grid on; ylabel('角速度 (rad/s)');
    if s == 1
        legend('实测','联合模型','Location','northeast');
    end

    nexttile;
    plot(S.tss,S.um,'Color',[0.55 0.20 0.65],'LineWidth',1.4);
    grid on; ylabel('实测力矩 (N·m)');
end
sgtitle(sprintf(['Pitch 全部独立角度阶跃：实测 vs 联合模型\n' ...
    'J=%.5f kg·m², b0=%.4f N·m·s/rad, tau_c=%.4f N·m, kLow=%.2f, kHigh=%.2f N·m/rad, tau_act=%.1f ms, Td=%.2f ms'], ...
    J,b0,tauC,kLow,kHigh,1000*tauAct,1000*Td));

figPath = fullfile(OUT_DIR,'pitch_step_all_csv_joint_fit_261019.png');
exportgraphics(fig,figPath,'Resolution',150);
close(fig);

fprintf('\n已保存:\n  %s\n  %s\n  %s\n', ...
    figPath,metricsPath,fullfile(OUT_DIR,'pitch_step_all_csv_joint_fit_261019.mat'));

%% ========================================================================
function candidates = extract_segments(csvPath)
candidates = {};
try
    T = readtable(csvPath,'VariableNamingRule','preserve','TextType','string');
catch
    return;
end
names = T.Properties.VariableNames;
if ~any(strcmp(names,'time')) || ...
        ~any(strcmp(names,'BPRS_log_pos_pitch [f32]')) || ...
        ~any(strcmp(names,'BPRS_log_w_pitch [f32]')) || ...
        ~any(strcmp(names,'BPRS_log_u_pitch [f32]'))
    return;
end

try
    tv = datetime(T.time,'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''','TimeZone','UTC');
catch
    tv = datetime(T.time,'TimeZone','UTC');
end
t = double(seconds(tv-tv(1)));
pos = numcol(T,'BPRS_log_pos_pitch [f32]');
vel = numcol(T,'BPRS_log_w_pitch [f32]');
tau = numcol(T,'BPRS_log_u_pitch [f32]');

[~,stem,ext] = fileparts(csvPath);
sourceName = string([stem ext]);

if any(strcmp(names,'pitch_angle_step_state [u8]')) && ...
        any(strcmp(names,'Gimbal.pitch.target_pos [f32]'))
    st = numcol(T,'pitch_angle_step_state [u8]');
    tgt = numcol(T,'Gimbal.pitch.target_pos [f32]');
    q = find(st == 2);
    if numel(q) >= 8
        q1 = find(st == 1 & t < t(q(1)));
    else
        q1 = [];
    end
    if numel(q) >= 8 && ~isempty(q1)
        i0 = q1(1); i1 = q(end);
        base = q1(max(1,numel(q1)-4):end);
        p0 = median(pos(base),'omitnan');
        uHold = median(tau(base),'omitnan');
        S = make_segment(t(i0:i1),pos(i0:i1),vel(i0:i1),tau(i0:i1), ...
            tgt(i0:i1),p0,uHold,sourceName,"mode7 step",0,char(string(tv(i0))));
        if S.nominalAmp > 0.02
            candidates{end+1} = S; %#ok<AGROW>
        end
    end
end

if any(strcmp(names,'pitch_multi_state [u8]')) && ...
        any(strcmp(names,'pitch_multi_index [u8]')) && ...
        any(strcmp(names,'pitch_multi_target [f32]'))
    st = numcol(T,'pitch_multi_state [u8]');
    idx = numcol(T,'pitch_multi_index [u8]');
    tgt = numcol(T,'pitch_multi_target [f32]');
    for k = 0:1
        % 只用基线和正向阶跃/保持段。回中段含明显静摩擦死区，
        % 普通 tanh 动摩擦模型无法在零速时维持力矩偏置，会把参数拉偏。
        q = find(ismember(st,[3 4]) & idx == k);
        if numel(q) < 8
            continue;
        end
        i0 = q(1); i1 = q(end);
        q4 = find(st == 4 & idx == k);
        if isempty(q4)
            continue;
        end
        base = q(t(q) < t(q4(1)));
        base = base(max(1,numel(base)-4):end);
        if isempty(base)
            base = q(1);
        end
        p0 = median(pos(base),'omitnan');
        uHold = median(tau(base),'omitnan');
        rampMs = 0;
        if k == 1
            rampMs = 300;
        end
        kindName = sprintf("multi idx%d",k);
        S = make_segment(t(i0:i1),pos(i0:i1),vel(i0:i1),tau(i0:i1), ...
            tgt(i0:i1),p0,uHold,sourceName,kindName,rampMs,char(string(tv(i0))));
        if S.nominalAmp > 0.02
            candidates{end+1} = S; %#ok<AGROW>
        end
    end
end
end

function S = make_segment(t,pos,vel,tau,ref,p0,uHold,source,kind,rampMs,startTime)
keep = isfinite(t) & isfinite(pos) & isfinite(vel) & isfinite(tau) & isfinite(ref);
t = t(keep); pos = pos(keep); vel = vel(keep); tau = tau(keep); ref = ref(keep);
[t,ia] = unique(t,'stable');
pos = pos(ia); vel = vel(ia); tau = tau(ia); ref = ref(ia);
t = t-t(1);

S = struct();
S.t = t;
S.y = pos-p0;
S.w = vel;
S.u = tau-uHold;
S.ref = ref-p0;
S.p0 = p0;
S.uHold = uHold;
S.y0 = S.y(1);
S.nominalAmp = max(ref,[],'omitnan')-p0;
S.source = char(source);
S.kind = char(kind);
S.rampMs = rampMs;
S.startTime = startTime;
S.n = numel(t);
S.dtMed = median(diff(t),'omitnan');
end

function S = prepare_segment(S,Ts)
keep = isfinite(S.t) & isfinite(S.y) & isfinite(S.w) & isfinite(S.u) & isfinite(S.ref);
t = S.t(keep); y = S.y(keep); w = S.w(keep); u = S.u(keep); ref = S.ref(keep);
[t,ia] = unique(t,'stable');
y = y(ia); w = w(ia); u = u(ia); ref = ref(ia);

tEnd = t(end);
tss = (0:Ts:tEnd)';
S.tss = tss;
S.ym = interp1(t,y,tss,'linear');
S.wm = interp1(t,w,tss,'linear');
S.um = interp1(t,u,tss,'linear');
S.refm = interp1(t,ref,tss,'linear');
S.w0 = w(1);
S.y0 = S.ym(1);
S.yScale = max(max(S.ym)-min(S.ym),0.02);
S.wScale = max(max(abs(S.wm)),0.20);
end

function c = joint_cost(x,segs,Ts,Td)
J=x(1); b0=x(2); tauC=x(3); kLow=x(4); kHigh=x(5); tauAct=x(6);
if J<0.005 || J>0.30 || b0<0.00 || b0>6 || ...
        tauC<0.02 || tauC>3 || kLow<2 || kLow>100 || kHigh<2 || kHigh>100 || ...
        tauAct<0.002 || tauAct>0.100
    c = 1e7;
    return;
end
c = 0;
for s = 1:numel(segs)
    S = segs{s};
    kk = local_k(S.p0,x);
    if kk<2 || kk>100
        c = 1e7;
        return;
    end
    [ys,ws] = sim_segment(S,x,Ts,Td);
    ey = ys-S.ym;
    ew = ws-S.wm;
    c = c + mean(ey.^2)/S.yScale^2 + ...
        0.03*mean(ew.^2)/S.wScale^2;
end
c = c/numel(segs);
end

function [yLog,wLog] = sim_segment(S,x,Ts,Td)
J=x(1); b0=x(2); tauC=x(3); kLow=x(4); kHigh=x(5); tauAct=x(6);
kk = local_k(S.p0,x);
N = numel(S.tss);
yLog = zeros(N,1); wLog = zeros(N,1);
y = S.y0; w = S.w0; tauA = 0;
sub = 4; h = Ts/sub;
for n = 1:N
    tt = S.tss(n);
    td = tt-Td;
    if td <= 0
        u = 0;
    else
        u = interp1(S.tss,S.um,min(td,S.tss(end)),'linear');
    end
    xx = [y;w;tauA];
    for q = 1:sub
        a1 = plant_fp(xx,u,J,b0,tauC,kk,tauAct);
        a2 = plant_fp(xx+0.5*h*a1,u,J,b0,tauC,kk,tauAct);
        a3 = plant_fp(xx+0.5*h*a2,u,J,b0,tauC,kk,tauAct);
        a4 = plant_fp(xx+h*a3,u,J,b0,tauC,kk,tauAct);
        xx = xx + (h/6)*(a1+2*a2+2*a3+a4);
    end
    y = xx(1); w = xx(2); tauA = xx(3);
    yLog(n) = y; wLog(n) = w;
end
end

function kk = local_k(p0,x)
if p0 < 0.22
    kk = x(4);
else
    kk = x(5);
end
end

function dx = plant_fp(x,tauCmd,J,b0,tauC,k,tauAct)
th = x(1); w = x(2); tauA = x(3);
dx = [w;
      (tauA-b0*w-tauC*tanh(w/0.03)-k*th)/J;
      (tauCmd-tauA)/tauAct];
end

function m = step_metrics(t,y,AMP)
m.rise = NaN; m.overshoot = NaN;
if ~isfinite(AMP) || AMP <= 0
    return;
end
i10 = find(y >= 0.1*AMP,1);
i90 = find(y >= 0.9*AMP,1);
if ~isempty(i10) && ~isempty(i90)
    m.rise = t(i90)-t(i10);
end
m.overshoot = max(0,(max(y)-AMP)/AMP);
end

function x = numcol(T,name)
v = T.(name);
if isnumeric(v)
    x = double(v);
else
    x = str2double(string(v));
end
end
