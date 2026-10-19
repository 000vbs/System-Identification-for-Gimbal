%% 分析 Pitch 闭环角度扫频 CSV：目标正弦 -> 实际角度的幅频/相频
clear; clc; close all;

here = fileparts(mfilename('fullpath'));
pitchRoot = fileparts(here);
outdir = fullfile(pitchRoot,'generated');
if ~exist(outdir,'dir'), mkdir(outdir); end
CSV = fullfile(pitchRoot,'data','emberprobe-live-20261019-002358.csv');
if ~isfile(CSV)
    error('请把 CSV 路径改成实际的 EmberProbe 导出文件');
end
T = readtable(CSV,'VariableNamingRule','preserve','TextType','string');
tv = datetime(T.time,'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''','TimeZone','UTC');
t = seconds(tv-tv(1));
st = double(T.("pitch_sweep_state [u8]"));
flt = double(T.("pitch_sweep_fault [u8]"));
idx = double(T.("pitch_sweep_index [u8]"));
ref = double(T.("pitch_sweep_ref [f32]"));
pos = double(T.("BPRS_log_pos_pitch [f32]"));
vel = double(T.("BPRS_log_w_pitch [f32]"));
tau = double(T.("BPRS_log_u_pitch [f32]"));
fcol = double(T.("pitch_sweep_freq [f32]"));

NOM = [1.00 2.00 3.00 5.00 7.00 10.00 14.00 20.00]; % 固件默认表（优先用 CSV 的 pitch_sweep_freq）

% CSV 可能混入之前几轮历史数据：只取最后一次 state=5 之后的扫频
q5 = find(st==5);
if isempty(q5), t0 = -inf; else, t0 = t(q5(end)); end

fprintf('采样 %d 点，时长 %.1f s，平均 %.1f Hz\n',numel(t),t(end),(numel(t)-1)/t(end));
fprintf(' state 最终 %d，fault 最终 %d\n', st(end), flt(end));

freq=[]; mag=[]; ph=[]; r2a=[]; r2r=[]; ampA=[];
for k=unique(idx)'
    q = find(st==2 & idx==k & t>=t0);
    if numel(q) < 12, continue; end
    tt = t(q)-t(q(1));
    rr = ref(q); pp = pos(q);
    f0 = median(fcol(q));
    if ~isfinite(f0) || f0 <= 0
        f0 = NOM(min(max(k+1,1),numel(NOM)));
    end
    % 去掉第一周期，避免起始瞬态
    keep = tt >= 1/f0;
    tt=tt(keep); rr=rr(keep); pp=pp(keep);
    good = isfinite(tt)&isfinite(rr)&isfinite(pp);
    tt=tt(good); rr=rr(good); pp=pp(good);
    if numel(tt)<8, continue; end
    % 频率细搜：用 ref 残差在标称 ±5% 内找
    fg=linspace(0.95*f0,1.05*f0,201); sse=inf(size(fg));
    for m=1:numel(fg)
        A=[ones(numel(tt),1),tt,sin(2*pi*fg(m)*tt),cos(2*pi*fg(m)*tt)];
        e=rr-A*(A\rr); sse(m)=e'*e;
    end
    [~,im]=min(sse); f=fg(im); w=2*pi*f;
    A=[ones(numel(tt),1),tt,sin(w*tt),cos(w*tt)];
    cr=A\rr; cp=A\pp;
    Gr=cr(3)+1i*cr(4); Gp=cp(3)+1i*cp(4);
    G=Gp/Gr;
    fit=A*cp;
    r2pos=1-sum((pp-fit).^2)/sum((pp-mean(pp)).^2);
    fitr=A*cr;
    r2ref=1-sum((rr-fitr).^2)/sum((rr-mean(rr)).^2);
    freq(end+1)=f; mag(end+1)=abs(G); ph(end+1)=angle(G)*180/pi; %#ok<SAGROW>
    r2a(end+1)=r2pos; r2r(end+1)=r2ref; %#ok<SAGROW>
    ampA(end+1)=abs(Gp); %#ok<SAGROW>
end

fprintf('\n f(Hz)   |G|      |G|(dB)   相位(deg)  R2ref  R2pos  Apos\n');
for i=1:numel(freq)
    fprintf('%6.3f  %8.4f  %8.2f  %+9.2f  %6.3f %6.3f  %.4f\n', ...
        freq(i),mag(i),20*log10(mag(i)),ph(i),r2r(i),r2a(i),ampA(i));
end

usable = ampA > 1e-3 & r2a > 0.80;
figure('Color','w','Position',[80 80 1000 700]);
subplot(2,1,1);
semilogx(freq(usable),20*log10(mag(usable)),'o-','LineWidth',1.5); hold on;
semilogx(freq(~usable),20*log10(mag(~usable)),'o','MarkerSize',6); grid on;
ylabel('|G| (dB)'); title('Pitch 闭环角度响应：目标角 \rightarrow 实际角');
subplot(2,1,2);
semilogx(freq(usable),ph(usable),'o-','LineWidth',1.5); hold on;
semilogx(freq(~usable),ph(~usable),'o','MarkerSize',6); grid on;
xlabel('频率 (Hz)'); ylabel('相位 (deg)');

saveas(gcf,fullfile(outdir,'pitch_angle_sweep_bode_261019.png'));
fprintf('\nBode 图已保存: pitch_angle_sweep_bode_261019.png\n');
