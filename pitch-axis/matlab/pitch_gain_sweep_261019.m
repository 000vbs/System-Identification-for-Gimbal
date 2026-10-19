%% Pitch 角度阶跃参数扫描（2026-10-19）
% 对象：J*th'' + b0*th' + tau_c*tanh(th'/wF) + k*th = tau(t-Td)
% 控制器（与固件一致）：
%   角度 P -> 目标角速度（±3 rad/s）
%   速度 PI（Kpw, Kiw；Iout 限 ±2 N·m；PID 输出限 ±3.5 N·m）
%   摩擦前馈 tau_c*sat(tv/0.3)
%   刚度前馈 kFF*(ref-center)，限 ±2 N·m
%   最终力矩限幅可选（0=不限）
%
% 目标：0.07 rad 阶跃，超调 <10%，尽量小的上升/稳定时间。
clear; clc; close all;

here = fileparts(mfilename('fullpath'));
pitchRoot = fileparts(here);
outdir = fullfile(pitchRoot,'generated');
if ~exist(outdir,'dir'), mkdir(outdir); end

% 植物：静态标定工作点 k≈33.67；动态宽频拟合 k≈15.78（用于鲁棒性检查）
P.J = 0.07008; P.b0 = 0.7784; P.tauC = 0.9391;
P.Td = 0.01625; P.wF = 0.03;
P.k = 33.67;

Ts = 1e-3; Tend = 3.0;
t = (0:Ts:Tend)';
N = numel(t);
AMP = 0.07;
ref = zeros(N,1);
ref(t >= 0.30 & t < 2.30) = AMP;   % 2s 保持，足够看稳定

% 参数网格
KpaList = [10 12 15 18 22];
KpwList = [1.0 1.5 2.0 2.5 3.0];
KiwList = [6 10 15 20];
kFFList = [0 15 25 33.67];
tauLimList = [0 3.5];

rows = {};
fprintf('扫描 %d 组...\n', numel(KpaList)*numel(KpwList)*numel(KiwList)*numel(kFFList)*numel(tauLimList));
for Kpa = KpaList
    for Kpw = KpwList
        for Kiw = KiwList
            for kFF = kFFList
                for tauLim = tauLimList
                    R = sim_case(P,t,ref,AMP,Kpa,Kpw,Kiw,kFF,tauLim,Ts);
                    rows(end+1,:) = {Kpa,Kpw,Kiw,kFF,tauLim, ...
                        R.tr,R.os,R.ts2,R.ts5,R.tauMax,R.velMax,R.sat,R.pidSat}; %#ok<SAGROW>
                end
            end
        end
    end
end

cols = {'Kpa','Kpw','Kiw','kFF','tauLim','tr','os','ts2','ts5','tauMax','velMax','sat','pidSat'};
S = cell2table(rows,'VariableNames',cols);

% 只保留超调 <10%、无最终力矩削顶的方案
ok = S.os < 0.10 & S.sat == 0 & isfinite(S.tr) & isfinite(S.ts2);
T = S(ok,:);
fprintf('\n满足超调<10%%且无削顶的方案: %d / %d\n', height(T), height(S));

% 评分：上升 + 2%稳定，超调接近 10% 给一点惩罚
T.score = T.tr + T.ts2 + 2.0*T.os;
T = sortrows(T,'score');

fprintf('\n%-6s %-6s %-6s %-7s %-7s %8s %8s %8s %8s %8s %8s %4s %5s\n', ...
    'Kpa','Kpw','Kiw','kFF','tauLim','tr(ms)','os(%)','ts2(ms)','ts5(ms)','tauPk','velPk','sat','pidSat');
for i = 1:min(20,height(T))
    r = T(i,:);
    fprintf('%-6.1f %-6.2f %-6.1f %-7.2f %-7.2f %8.1f %8.1f %8.1f %8.1f %8.2f %8.2f %4d %5d\n', ...
        r.Kpa,r.Kpw,r.Kiw,r.kFF,r.tauLim,1000*r.tr,100*r.os,1000*r.ts2,1000*r.ts5,r.tauMax,r.velMax,r.sat,r.pidSat);
end

% 前 15 名同时跑 k=15.78 与 k=33.67，挑两边超调都 <10% 的稳健解
fprintf('\n前 15 名鲁棒性检查（k=15.78 / k=33.67）：\n');
fprintf('%-6s %-6s %-6s %-7s %-7s | %8s %8s %8s %8s\n', ...
    'Kpa','Kpw','Kiw','kFF','tauLim','tr15','os15','tr33','os33');
robust = [];
for i = 1:min(15,height(T))
    r = T(i,:);
    P2=P; P2.k=15.78; R15=sim_case(P2,t,ref,AMP,r.Kpa,r.Kpw,r.Kiw,r.kFF,r.tauLim,Ts);
    P3=P; P3.k=33.67; R33=sim_case(P3,t,ref,AMP,r.Kpa,r.Kpw,r.Kiw,r.kFF,r.tauLim,Ts);
    fprintf('%-6.1f %-6.2f %-6.1f %-7.2f %-7.2f | %8.0f %8.1f %8.0f %8.1f\n', ...
        r.Kpa,r.Kpw,r.Kiw,r.kFF,r.tauLim,1000*R15.tr,100*R15.os,1000*R33.tr,100*R33.os);
    if R15.os < 0.10 && R33.os < 0.10 && R15.sat==0 && R33.sat==0
        robust(end+1,:) = [r.Kpa,r.Kpw,r.Kiw,r.kFF,r.tauLim, ...
            R15.tr,R15.os,R15.ts2,R33.tr,R33.os,R33.ts2,max(R15.tauMax,R33.tauMax)]; %#ok<SAGROW>
    end
end

if ~isempty(robust)
    fprintf('\n两边都 <10%% 的稳健候选：\n');
    fprintf('%-6s %-6s %-6s %-7s %-7s %8s %8s %8s %8s %8s %8s %8s\n', ...
        'Kpa','Kpw','Kiw','kFF','tauLim','tr15','os15','ts15','tr33','os33','ts33','tauPk');
    for i=1:size(robust,1)
        fprintf('%-6.1f %-6.2f %-6.1f %-7.2f %-7.2f %8.0f %8.1f %8.0f %8.0f %8.1f %8.0f %8.2f\n', ...
            robust(i,1),robust(i,2),robust(i,3),robust(i,4),robust(i,5), ...
            1000*robust(i,6),100*robust(i,7),1000*robust(i,8), ...
            1000*robust(i,9),100*robust(i,10),1000*robust(i,11),robust(i,12));
    end
else
    fprintf('\n前 15 名里没有两边都 <10%% 的候选\n');
end

save(fullfile(outdir,'pitch_gain_sweep_261019.mat'),'S','T','P');

%% 本地函数
function R = sim_case(P,t,ref,AMP,Kpa,Kpw,Kiw,kFF,tauLim,Ts)
    N = numel(t); th=0; w=0; Iout=0;
    Ndelay = max(1,round(P.Td/Ts));
    delay = zeros(Ndelay,1); tauApply=0;
    thLog=zeros(N,1);
    sat=0; pidSat=0; tauMax=0; velMax=0;
    KiCode = Kiw*Ts;
    sub=10; h=Ts/sub;
    ffMax=2.0; pidMax=3.5; iMax=2.0; angVelMax=3.0;

    for n=1:N
        % 角度 P
        tv = Kpa*(ref(n)-th);
        tv = max(-angVelMax,min(angVelMax,tv));
        % 速度 PI
        ev = tv - w;
        Iout = Iout + KiCode*ev;
        Iout = max(-iMax,min(iMax,Iout));
        pidRaw = Kpw*ev + Iout;
        if abs(pidRaw) >= pidMax, pidSat = pidSat + 1; end
        pidOut = max(-pidMax,min(pidMax,pidRaw));
        % 前馈
        s = max(-1,min(1,tv/0.30));
        tauFF = P.tauC*s;
        tauSpring = max(-ffMax,min(ffMax,kFF*ref(n)));
        raw = pidOut + tauFF + tauSpring;
        if tauLim > 0
            if abs(raw) >= tauLim, sat=sat+1; end
            raw = max(-tauLim,min(tauLim,raw));
        end
        if abs(raw) > tauMax, tauMax = abs(raw); end
        % 延迟
        tauApply = delay(1);
        delay = [delay(2:end); raw];
        % 植物 RK4
        x=[th;w];
        for q=1:sub
            k1=fp(x,tauApply,P);
            k2=fp(x+0.5*h*k1,tauApply,P);
            k3=fp(x+0.5*h*k2,tauApply,P);
            k4=fp(x+h*k3,tauApply,P);
            x=x+(h/6)*(k1+2*k2+2*k3+k4);
        end
        th=x(1); w=x(2);
        if abs(w) > velMax, velMax = abs(w); end
        thLog(n)=th;
    end

    idx = t>=0.30 & t<2.30;
    tt = t(idx)-0.30; yy=thLog(idx);
    i10=find(yy>=0.1*AMP,1); i90=find(yy>=0.9*AMP,1);
    if isempty(i10)||isempty(i90), tr=NaN; else, tr=tt(i90)-tt(i10); end
    os=max(0,(max(yy)-AMP)/AMP);
    b2=0.02*AMP; out2=find(abs(yy-AMP)>b2,1,'last');
    if isempty(out2), ts2=0; else, ts2=tt(out2); end
    b5=0.05*AMP; out5=find(abs(yy-AMP)>b5,1,'last');
    if isempty(out5), ts5=0; else, ts5=tt(out5); end
    R=struct('tr',tr,'os',os,'ts2',ts2,'ts5',ts5, ...
        'tauMax',tauMax,'velMax',velMax,'sat',sat,'pidSat',pidSat);
end

function dx=fp(x,tau,P)
    th=x(1); w=x(2);
    dx=[w;(tau-P.b0*w-P.tauC*tanh(w/P.wF)-P.k*th)/P.J];
end
