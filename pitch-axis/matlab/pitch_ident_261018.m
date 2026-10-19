%% Pitch 轴开环正弦扫频：数据质量审查 + 灰箱辨识（2026-10-18）
% 目的：把"1.5~3.2Hz 六点 + 1.0~1.3Hz 加密复测"放在同一套判据下处理，
%       给出可解释的 J / k / 摩擦模型，并明确哪些频点不能用于线性拟合。
%
% 关键修正（相对 Desktop/pitch_greybox_fit.m）：
%   1) 不再用 u 的上升过零点估计频率。u 带 -0.21N·m 直流偏置，
%      且部分频点记录被截短，过零点判据会直接跳掉 1.6/2.0Hz 两个点。
%      改为在标称频率 ±5% 内对 u 做 [1,t,sin,cos] 网格搜索。
%   2) 每段只用"斜坡结束后、记录结束前"的数据，并同时算 R²、三次谐波、
%      正负不对称和净漂移，先判断该点是否粘滑/顶限位。
%   3) 摩擦不当作粘性阻尼：b_eff = b0 + c/|Omega|，c = 4*tau_c/pi。
%      同时给出 tau_c 的描述函数估计与 b_eff-1/|Omega| 回归交叉验证。

clear; clc; close all;

here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
pitchRoot = fileparts(here);
outdir = fullfile(pitchRoot,'generated');
if ~exist(outdir,'dir'), mkdir(outdir); end

CSVS = { ...
    fullfile(pitchRoot,'data','pitch_rec_6pt.csv'), ...
    fullfile(pitchRoot,'data','pitch_rec_low.csv')};
RAMP_CYC = 3.0;
CYC_PER_FQ = 15;

%% ---------------- 读数据并按连续块切段 ----------------
rec = struct('run',{},'idx',{},'fNom',{},'t',{},'u',{},'pMot',{},'pImu',{});
for r = 1:numel(CSVS)
    if ~isfile(CSVS{r})
        warning('找不到 %s，跳过', CSVS{r});
        continue;
    end
    T = readtable(CSVS{r}, 'VariableNamingRule','preserve','TextType','string');
    tv = datetime(T.time, 'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''','TimeZone','UTC');
    idx = double(T.("sweep_ol_freq_index [u8]"));
    u   = double(T.("sweep_ol_u [f32]"));
    pM  = double(T.("sweep_ol_unwrap_pos [f32]"));
    pI  = double(T.("sweep_ol_imu_pitch [f32]"));
    tAll = seconds(tv - tv(1));

    % 连续块边界：频点索引变化或时间倒退都切开
    brk = [1; find(diff(idx) ~= 0 | diff(tAll) <= 0) + 1; numel(idx)+1];
    for b = 1:numel(brk)-1
        rows = brk(b):brk(b+1)-1;
        k = idx(rows(1));
        rec(end+1) = struct('run',r,'idx',k,'fNom',NaN, ...
            't',tAll(rows)-tAll(rows(1)),'u',u(rows),'pMot',pM(rows),'pImu',pI(rows)); %#ok<SAGROW>
    end
end

%% ---------------- 逐段：频率精估 + 复增益 + 数据质量 ----------------
rows = {};
for n = 1:numel(rec)
    R = rec(n);
    t = R.t; u = R.u; pM = R.pMot; pI = R.pImu;

    % 第一步：不依赖索引名，先用 u 的宽带残差搜索找主频。
    % 6pt 文件是旧表 1.5~7.5Hz，low 文件是新表 1.0~3.2Hz；CSV 里的索引
    % 在两版固件里含义不同，直接用索引查表会把频点整体标错。
    good0 = isfinite(t) & isfinite(u) & isfinite(pM) & isfinite(pI) & t > 0.2;
    if nnz(good0) < 40
        fprintf('run%d idx%d 点数不足，跳过\n', R.run, R.idx);
        continue;
    end
    t0 = t(good0); u0 = u(good0);
    fGrid0 = linspace(0.5, 8.5, 1601);
    sse0 = inf(size(fGrid0));
    for q = 1:numel(fGrid0)
        A0 = [ones(numel(t0),1), t0, sin(2*pi*fGrid0(q)*t0), cos(2*pi*fGrid0(q)*t0)];
        r0 = u0 - A0*(A0\u0);
        sse0(q) = r0'*r0;
    end
    [~,iq0] = min(sse0);
    f0 = fGrid0(iq0);

    % 只保留斜坡结束后、最多 CYC_PER_FQ 个周期的数据
    tRamp = RAMP_CYC / f0;
    tEnd  = CYC_PER_FQ / f0;
    good  = isfinite(t) & isfinite(u) & isfinite(pM) & isfinite(pI) & ...
            t >= tRamp & t <= tEnd;
    if nnz(good) < 40
        fprintf('run%d idx%d 有效点不足(%d)，跳过\n', R.run, R.idx, nnz(good));
        continue;
    end
    t = t(good); u = u(good); pM = pM(good); pI = pI(good);
    t = t - t(1);
    N = numel(t);

    % --- 频率：u 的残差网格搜索 ---
    fGrid = linspace(0.98*f0, 1.02*f0, 401);
    sse = inf(size(fGrid));
    for q = 1:numel(fGrid)
        A = [ones(N,1), t, sin(2*pi*fGrid(q)*t), cos(2*pi*fGrid(q)*t)];
        ru = u - A*(A\u);
        sse(q) = ru'*ru;
    end
    [~,iq] = min(sse);
    fEst = fGrid(iq);
    w = 2*pi*fEst;

    % --- 基波拟合 ---
    A = [ones(N,1), t, sin(w*t), cos(w*t)];
    cu = A\u; cpM = A\pM; cpI = A\pI;
    Gu  = cu(3) + 1i*cu(4);
    GM  = (cpM(3) + 1i*cpM(4))/Gu;
    GI  = (cpI(3) + 1i*cpI(4))/Gu;

    % --- 拟合质量 ---
    fitM = A*cpM; fitI = A*cpI;
    R2M = 1 - sum((pM-fitM).^2)/sum((pM-mean(pM)).^2);
    R2I = 1 - sum((pI-fitI).^2)/sum((pI-mean(pI)).^2);

    % 三次谐波（把基波+趋势+直流去掉后，再拟 3w）
    rM = pM - fitM; rI = pI - fitI;
    A3 = [ones(N,1), sin(3*w*t), cos(3*w*t)];
    c3M = A3\rM; c3I = A3\rI;
    h3M = hypot(c3M(2),c3M(3)) / hypot(cpM(3),cpM(4));
    h3I = hypot(c3I(2),c3I(3)) / hypot(cpI(3),cpI(4));

    % 正负半周幅值不对称（按基波相位分半周，取去趋势后的峰谷）
    ph = angle(cpI(3) + 1i*cpI(4));
    s = sin(w*t + ph);
    posAmp = max(rI(s >= 0)); negAmp = min(rI(s < 0));
    asym = (posAmp + negAmp) / max(posAmp - negAmp, eps);

    rows(end+1,:) = {R.run, R.idx, fEst, f0, N, abs(Gu), ...
        abs(GM), angle(GM)*180/pi, abs(GI), angle(GI)*180/pi, ...
        R2M, R2I, h3M, h3I, asym, cpI(2), cpM(2)}; %#ok<SAGROW>
end

if isempty(rows)
    error('没有任何可用数据段');
end
cols = {'run','idx','fEst','fNom','N','uAmp','Gmot','phMot','Gimu','phImu', ...
        'R2mot','R2imu','h3mot','h3imu','asym','driftImu','driftMot'};
P = cell2table(rows, 'VariableNames', cols);
P.valid = P.R2imu > 0.85 & P.h3imu <= 0.20;

fprintf('\n================ 逐段质量审查 ================\n');
fprintf('run idx  fEst   uAmp  |Gm|    phM     |Gi|    phI    R2i   h3i   asym  drift  有效\n');
for n = 1:height(P)
    fprintf('%3d %3d %5.3f %5.2f %7.5f %+7.2f %7.5f %+7.2f %5.3f %5.1f%% %+5.1f%% %+6.3f  %d\n', ...
        P.run(n),P.idx(n),P.fEst(n),P.uAmp(n),P.Gmot(n),P.phMot(n),P.Gimu(n),P.phImu(n), ...
        P.R2imu(n),100*P.h3imu(n),100*P.asym(n),P.driftImu(n),P.valid(n));
end

% 同一轮里索引才可以比较；不同固件版本的同名索引对应不同频率。
% 因此这里按"实测频率"聚类，而不是按 idx 去重。
fprintf('\n重复测量（频率相差 <2%% 视为同一点）：\n');
fk = P.fEst;
used = false(height(P),1);
for n = 1:height(P)
    if used(n) || ~P.valid(n), continue; end
    q = find(P.valid & abs(fk-fk(n))./fk(n) < 0.02);
    if numel(q) >= 2
        fprintf('  %.3f Hz: |Gi| %s  相位 %s\n', fk(n), ...
            mat2str(P.Gimu(q)',4), mat2str(P.phImu(q)',5));
    end
    used(q) = true;
end

%% ---------------- 选点：所有有效频率点，按频率排序 ----------------
P.phaseUnwrap = rad2deg(unwrap(deg2rad(P.phImu)));
sel = P.valid;
Q = sortrows(P(sel,:), 'fEst');
fprintf('\n用于灰箱拟合的频点：\n');
disp(Q(:,{'run','idx','fEst','uAmp','Gimu','phImu','R2imu','h3imu'}));

if height(Q) < 3
    error('有效频点不足 3 个，无法拟合');
end

w = 2*pi*Q.fEst;
GI = Q.Gimu .* exp(1i*deg2rad(Q.phImu));
GM = Q.Gmot .* exp(1i*deg2rad(Q.phMot));
Omega = w .* abs(GI) .* Q.uAmp;       % IMU 侧速度幅值 (rad/s)
Z  = 1 ./ GI;                         % 位置通道: Re(Z)=k-Jw^2, Im(Z)=w*b
Zv = 1 ./ (1i*w.*GI);                 % 速度通道: Re(Zv)=b_eff, Im(Zv)=w*J-k/w

%% ---------------- 模型 A：纯 LTI 二阶（IMU 侧） ----------------
fprintf('\n================ 模型 A: 1/(J s^2 + b s + k)，IMU 侧 ================\n');
AA = [-w.^2, zeros(size(w)), ones(size(w)); ...
      zeros(size(w)), w, zeros(size(w))];
bb = [real(Z); imag(Z)];
x = AA\bb;
JA = x(1); bA = x(2); kA = x(3);
GA = 1./(-JA*w.^2 + 1i*w*bA + kA);
errA = abs(GA-GI)./abs(GI);
fprintf('  J=%.5f  b=%.4f  k=%.3f   固有频率 %.3f Hz\n', JA,bA,kA,sqrt(kA/JA)/(2*pi));
fprintf('  逐点相对误差: %s%%  中位 %.1f%%\n', mat2str(100*errA',3), 100*median(errA));

%% ---------------- 模型 B：库仑摩擦描述函数 ----------------
% Re(1/G) = b0 + c/Omega ; Im(1/G) = w*J - k/w
fprintf('\n================ 模型 B: b_eff = b0 + c/|Omega|，c = 4*tau_c/pi ================\n');
beff = real(Zv);
AB = [ones(numel(Omega),1), 1./Omega(:)];
xb = AB\beff;
b0 = xb(1); cC = xb(2); tauC = pi*cC/4;
% J,k 用虚部线性解（b 不再混进虚部）
AB2 = [w(:), -1./w(:)];
xJk = AB2\imag(Zv);
JB = xJk(1); kB = xJk(2);
GB = 1 ./ (-JB*w.^2 + 1i*w.*(b0 + cC./Omega(:)) + kB);
errB = abs(GB-GI)./abs(GI);
fprintf('  J=%.5f  k=%.3f  b0=%.4f  c=%.4f  tau_c=%.4f N*m\n', JB,kB,b0,cC,tauC);
fprintf('  b_eff 回归 R^2 = %.4f\n', 1-sum((beff-AB*xb).^2)/sum((beff-mean(beff)).^2));
fprintf('  逐点相对误差: %s%%  中位 %.1f%%\n', mat2str(100*errB',3), 100*median(errB));

%% ---------------- 模型 B+延迟：摩擦描述函数 + 纯延迟 ----------------
% 先扫 Td，每个 Td 下去掉延迟后用线性最小二乘解 (J,k,b0,c)，
% 避免把"结构失配"全都塞给一个非物理的负粘性阻尼。
fprintf('\n================ 模型 B+延迟: 每 Td 重解 (J,k,b0,tau_c) ================\n');
bestTdB = 0; bestErrB = inf; bestXB = [];
for Td = 0:0.00025:0.020
    Gc = GI .* exp(1i*w*Td);
    Zvc = 1 ./ (1i*w.*Gc);
    xbq = [ones(numel(Omega),1), 1./Omega(:)] \ real(Zvc);
    xjq = [w(:), -1./w(:)] \ imag(Zvc);
    Gq = exp(-1i*w*Td) ./ (-xjq(1)*w.^2 + 1i*w.*(xbq(1)+xbq(2)./Omega(:)) + xjq(2));
    eq = median(abs(Gq-GI)./abs(GI));
    if eq < bestErrB
        bestErrB = eq; bestTdB = Td;
        bestXB = [xjq(1), xjq(2), xbq(1), xbq(2)];
    end
end
JBd=bestXB(1); kBd=bestXB(2); b0d=bestXB(3); cCd=bestXB(4); tauCd=pi*cCd/4;
GBd = exp(-1i*w*bestTdB) ./ (-JBd*w.^2 + 1i*w.*(b0d + cCd./Omega(:)) + kBd);
errBd = abs(GBd-GI)./abs(GI);
fprintf('  Td=%.2f ms  J=%.5f  k=%.3f  b0=%.4f  tau_c=%.4f N*m\n', ...
    1000*bestTdB, JBd,kBd,b0d,tauCd);
fprintf('  逐点相对误差: %s%%  中位 %.1f%%\n', mat2str(100*errBd',3), 100*median(errBd));

%% ---------------- 模型 C：双惯量（电机-弹簧-负载） ----------------
fprintf('\n================ 模型 C: 双惯量，电机侧 + IMU 侧联合 ================\n');
opt = optimset('MaxFunEvals',1e5,'MaxIter',1e5,'Display','off','TolFun',1e-12,'TolX',1e-12);
starts = {log([0.08 0.10 8 0.3 0.3]), log([0.05 0.20 3 1.0 0.5]), ...
          log([0.15 0.15 20 0.2 0.2]), log([0.10 0.05 2 0.5 1.0])};
best = []; bestc = inf;
for q = 1:numel(starts)
    xq = fminsearch(@(x) cost_dual(x,w,GM,GI), starts{q}, opt);
    cq = cost_dual(xq,w,GM,GI);
    if cq < bestc, best = xq; bestc = cq; end
end
pp = exp(best);
Jm=pp(1); Jl=pp(2); ks=pp(3); cs=pp(4); bl=pp(5);
GCm = zeros(size(w)); GCl = zeros(size(w));
for q = 1:numel(w)
    s = 1i*w(q); K = ks + cs*s;
    den = (Jm*s^2+K)*(Jl*s^2+bl*s+K) - K^2;
    GCm(q) = (Jl*s^2+bl*s+K)/den;
    GCl(q) = K/den;
end
errC = [abs(GCm-GM)./abs(GM); abs(GCl-GI)./abs(GI)];
fprintf('  Jm=%.5f  Jl=%.5f  ks=%.3f  cs=%.4f  bl=%.4f\n',Jm,Jl,ks,cs,bl);
fprintf('  反共振 %.3f Hz   共振 %.3f Hz\n', ...
    sqrt(ks/Jl)/(2*pi), sqrt(ks*(Jm+Jl)/(Jm*Jl))/(2*pi));
fprintf('  联合中位误差 %.1f%%\n',100*median(errC));

%% ---------------- 延迟扫描 ----------------
fprintf('\n================ 延迟扫描（0~20ms） ================\n');
bestTd = 0; bestEc = inf;
for Td = 0:0.0005:0.020
    Gm = GI .* exp(1i*w*Td);  Zt = 1./Gm;
    xt = AA\ [real(Zt); imag(Zt)];
    Gt = 1./(-xt(1)*w.^2 + 1i*w*xt(2) + xt(3));
    ec = median(abs(Gt-GI)./abs(GI));
    if ec < bestEc, bestEc = ec; bestTd = Td; end
end
fprintf('  最佳 Td = %.2f ms，对应中位误差 %.1f%%（Td=0 时 %.1f%%）\n', ...
    1000*bestTd, 100*bestEc, 100*median(errA));

%% ---------------- 交叉验证：留一 ----------------
fprintf('\n================ 留一交叉验证（模型 B） ================\n');
errLOO = zeros(height(Q),1);
for q = 1:height(Q)
    keep = true(height(Q),1); keep(q) = false;
    Ab = [ones(nnz(keep),1), 1./Omega(keep)];
    xbq = Ab\beff(keep);
    Jkq = [w(keep), -1./w(keep)]\imag(Zv(keep));
    Gq = 1/(-Jkq(1)*w(q)^2 + 1i*w(q)*(xbq(1)+xbq(2)/Omega(q)) + Jkq(2));
    errLOO(q) = abs(Gq-GI(q))/abs(GI(q));
end
fprintf('  留一误差: %s%%  中位 %.1f%%\n', mat2str(100*errLOO',3),100*median(errLOO));

% 模型 B+延迟的留一：每个点都用其余点重新扫 Td 并重解参数
errLOOd = zeros(height(Q),1);
for q = 1:height(Q)
    keep = true(height(Q),1); keep(q) = false;
    bestq = inf; bestpar = [];
    for Td = 0:0.00025:0.020
        Gc = GI(keep) .* exp(1i*w(keep)*Td);
        Zvq = 1 ./ (1i*w(keep).*Gc);
        xbq = [ones(nnz(keep),1), 1./Omega(keep)] \ real(Zvq);
        xjq = [w(keep), -1./w(keep)] \ imag(Zvq);
        Gfit = exp(-1i*w(keep)*Td) ./ (-xjq(1)*w(keep).^2 + ...
            1i*w(keep).*(xbq(1)+xbq(2)./Omega(keep)) + xjq(2));
        eq = median(abs(Gfit-Gc)./abs(Gc));
        if eq < bestq
            bestq = eq; bestpar = [xjq(1), xjq(2), xbq(1), xbq(2), Td];
        end
    end
    Gq = exp(-1i*w(q)*bestpar(5)) ./ (-bestpar(1)*w(q)^2 + ...
        1i*w(q)*(bestpar(3)+bestpar(4)/Omega(q)) + bestpar(2));
    errLOOd(q) = abs(Gq-GI(q))/abs(GI(q));
end
fprintf('  B+延迟 留一误差: %s%%  中位 %.1f%%\n', mat2str(100*errLOOd',3),100*median(errLOOd));

%% ---------------- 图 ----------------
f = figure('Position',[80 80 1100 760],'Color','w');

subplot(2,2,[1 2]);
semilogx(w/(2*pi),20*log10(abs(GI)),'ko','MarkerFaceColor','k'); hold on;
semilogx(w/(2*pi),20*log10(abs(GBd)),'b-','LineWidth',1.6);
semilogx(w/(2*pi),20*log10(abs(GB)),'b:','LineWidth',1.0);
semilogx(w/(2*pi),20*log10(abs(GA)),'r--','LineWidth',1.2);
grid on; ylabel('|G| (dB)'); xlabel('频率 (Hz)');
legend('实测 IMU','库仑摩擦+延迟','库仑摩擦(无延迟)','LTI 模型','Location','southwest');
title('Pitch 开环力矩 \rightarrow IMU 角度');

subplot(2,2,3);
semilogx(w/(2*pi),Q.phaseUnwrap,'ko','MarkerFaceColor','k'); hold on;
semilogx(w/(2*pi),unwrap(angle(GBd))*180/pi,'b-','LineWidth',1.6);
semilogx(w/(2*pi),unwrap(angle(GB))*180/pi,'b:','LineWidth',1.0);
semilogx(w/(2*pi),unwrap(angle(GA))*180/pi,'r--','LineWidth',1.2);
grid on; ylabel('相位 (deg)'); xlabel('频率 (Hz)');
legend('实测','库仑摩擦+延迟','库仑摩擦','LTI','Location','southwest');

subplot(2,2,4);
beff_d = real(1 ./ (1i*w.*GI.*exp(1i*w*bestTdB)));
ABd = [ones(numel(Omega),1), 1./Omega(:)];
xbd = ABd\beff_d;
plot(1./Omega,beff_d,'ko','MarkerFaceColor','k'); hold on;
plot(1./Omega,ABd*xbd,'b-','LineWidth',1.5);
grid on; xlabel('1/|Omega| (s/rad)'); ylabel('Re(1/G) = b_{eff}');
title(sprintf('去延迟后：b_0=%.3f, \\tau_c=%.3f N\\cdotm',xbd(1),pi*xbd(2)/4));

saveas(f, fullfile(outdir,'pitch_ident_261018_bode.png'));

% 数据质量图
f2 = figure('Position',[120 120 1000 600],'Color','w');
subplot(2,1,1);
bar(P.idx + 0.15*(P.run-1), 100*P.h3imu);
grid on; ylabel('IMU 三次谐波 (%)'); xlabel('频点索引');
title('粘滑判据：三次谐波/基波（>12% 不可用）');
subplot(2,1,2);
bar(P.idx + 0.15*(P.run-1), P.driftImu);
grid on; ylabel('IMU 净漂移 (rad/s)'); xlabel('频点索引');
title('低频蠕动/顶限位判据');
saveas(f2, fullfile(outdir,'pitch_ident_261018_quality.png'));

save(fullfile(outdir,'pitch_ident_261018.mat'),'P','Q','JA','bA','kA','JB','kB','b0','cC','tauC', ...
    'Jm','Jl','ks','cs','bl','bestTd','errA','errB','errLOO');
writetable(P, fullfile(outdir,'pitch_ident_261018_points.csv'));

fprintf('\n结果已保存到 %s\n', outdir);

%% ---------------- 本地函数 ----------------
function J = cost_dual(x, w, gm, gl)
    p = exp(x);
    Jm=p(1); Jl=p(2); ks=p(3); cs=p(4); bl=p(5);
    J = 0;
    for n = 1:numel(w)
        s = 1i*w(n); K = ks + cs*s;
        den = (Jm*s^2+K)*(Jl*s^2+bl*s+K) - K^2;
        gmh = (Jl*s^2+bl*s+K)/den;
        glh = K/den;
        d1 = log(gmh/gm(n)); d2 = log(glh/gl(n));
        J = J + real(d1)^2 + imag(d1)^2 + real(d2)^2 + imag(d2)^2;
    end
end
