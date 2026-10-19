%% 用 0.1 阶跃 + 0.3 斜坡实测拟合 0.15 工作点局部植物模型
% 输入用实测力矩（含控制器前馈/PID），只拟合植物 J,b0,tau_c,k
clear; clc; close all;

here = fileparts(mfilename('fullpath'));
pitchRoot = fileparts(here);
outdir = fullfile(pitchRoot,'generated');
if ~exist(outdir,'dir'), mkdir(outdir); end
CSV=fullfile(pitchRoot,'data','emberprobe-live-20261019-210710.csv');
T=readtable(CSV,'VariableNamingRule','preserve','TextType','string');
tv=datetime(T.time,'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''','TimeZone','UTC');
t=seconds(tv-tv(1));
st=str2double(string(T.("pitch_multi_state [u8]")));
idx=str2double(string(T.("pitch_multi_index [u8]")));
pos=str2double(string(T.("BPRS_log_pos_pitch [f32]")));
vel=str2double(string(T.("BPRS_log_w_pitch [f32]")));
tau=str2double(string(T.("BPRS_log_u_pitch [f32]")));
fault=str2double(string(T.("pitch_multi_fault [u8]")));

segs={};
for k=0:1
    q=find(st==4 & idx==k); i0=q(1); i1=q(end);
    S=struct('tm',t(i0:i1)-t(i0),'pm',pos(i0:i1)-pos(i0), ...
        'wm',vel(i0:i1),'um',tau(i0:i1),'p0',pos(i0),'w0',vel(i0));
    segs{end+1}=S; %#ok<SAGROW>
end

Ts=1e-3; Td=0.01625;
x0s=[0.07 0.8 0.15 10; 0.09 1.5 0.05 30; 0.07 0.8 0.94 15.78];
opt=optimset('MaxFunEvals',800,'MaxIter',800,'Display','off','TolFun',1e-8,'TolX',1e-8);
best=[];bc=inf;
for s=1:size(x0s,1)
    xb=fminsearch(@(x) cost(x,segs,Ts,Td),x0s(s,:),opt);
    cb=cost(xb,segs,Ts,Td);
    if cb<bc,best=xb;bc=cb;end
end
J=best(1);b0=best(2);tauC=best(3);k=best(4);
fprintf('\n局部模型: J=%.5f b0=%.4f tau_c=%.4f k=%.3f\n',J,b0,tauC,k);

for s=1:numel(segs)
    S=segs{s};
    [thS,wS]=sim_plant(S,J,b0,tauC,k,Ts,Td);
    tss=(0:Ts:S.tm(end))';
    thAt=interp1(tss,thS,S.tm); wAt=interp1(tss,wS,S.tm);
    AMP=max(abs(thS));
    nrmse=sqrt(mean((thAt-S.pm).^2))/AMP;
    fprintf(' idx%d: 角度 NRMSE %.2f%%，峰值实际 %.4f/仿真 %.4f\n', ...
        s-1,100*nrmse,max(abs(S.pm)),AMP);
    segs{s}.thS=thS; segs{s}.wS=wS; segs{s}.tss=tss;
end

figure('Color','w','Position',[70 50 1100 800]);
for s=1:numel(segs)
    S=segs{s};
    subplot(2,2,s);
    plot(S.tm,S.pm*180/pi,'o-','Color',[0.05 0.35 0.85],'LineWidth',1.2,'MarkerSize',4); hold on;
    plot(S.tss,S.thS*180/pi,'-','Color',[0.85 0.25 0.15],'LineWidth',1.5);
    grid on; ylabel('角度 (deg)'); if s==1,title('0.1 rad 阶跃：实测 vs 局部模型');else,title('0.3 rad 斜坡：实测 vs 局部模型');end
    legend('实测','局部模型','Location','southeast'); xlim([0 S.tm(end)]);
    subplot(2,2,s+2);
    plot(S.tm,S.wm,'o-','Color',[0.05 0.35 0.85],'LineWidth',1.2,'MarkerSize',4); hold on;
    plot(S.tss,S.wS,'-','Color',[0.85 0.25 0.15],'LineWidth',1.5);
    grid on; xlabel('段内时间 (s)'); ylabel('角速度 (rad/s)');
    legend('实测','局部模型','Location','northeast'); xlim([0 S.tm(end)]);
end
sgtitle(sprintf('Pitch 0.1/0.3 多幅值阶跃：实测 vs 0.15 工作点局部模型（fault=%d）',fault(end)));
saveas(gcf,fullfile(outdir,'pitch_multi_fit_261019.png'));
fprintf('图已保存: pitch_multi_fit_261019.png\n');

%% ---------- 本地函数 ----------
function c=cost(x,segs,Ts,Td)
    J=x(1);b0=x(2);tauC=x(3);k=x(4);
    if J<0.02||J>0.3||b0<0.05||b0>6||tauC<0.02||tauC>3||k<2||k>60
        c=1e6;return;
    end
    c=0;
    for s=1:numel(segs)
        S=segs{s};
        [th,w]=sim_plant(S,J,b0,tauC,k,Ts,Td);
        tss=(0:Ts:S.tm(end))';
        thAt=interp1(tss,th,S.tm); wAt=interp1(tss,w,S.tm);
        c=c+mean((thAt-S.pm).^2)/var(S.pm) + 0.05*mean((wAt-S.wm).^2)/var(S.wm);
    end
end

function [thLog,wLog]=sim_plant(S,J,b0,tauC,k,Ts,Td)
    N=round(S.tm(end)/Ts)+1; th=0; w=S.w0;
    thLog=zeros(N,1);wLog=zeros(N,1);
    sub=1;h=Ts/sub;
    for n=1:N
        tt=(n-1)*Ts;
        td=tt-Td;
        if td<=S.tm(1), u=S.um(1); else, u=interp1(S.tm,S.um,min(td,S.tm(end)),'linear'); end
        xx=[th;w];
        for q=1:sub
            k1=fp(xx,u,J,b0,tauC,k);
            k2=fp(xx+0.5*h*k1,u,J,b0,tauC,k);
            k3=fp(xx+0.5*h*k2,u,J,b0,tauC,k);
            k4=fp(xx+h*k3,u,J,b0,tauC,k);
            xx=xx+(h/6)*(k1+2*k2+2*k3+k4);
        end
        th=xx(1);w=xx(2);thLog(n)=th;wLog(n)=w;
    end
end

function dx=fp(x,tau,J,b0,tauC,k)
    th=x(1);w=x(2);
    dx=[w;(tau-b0*w-tauC*tanh(w/0.03)-k*th)/J];
end
