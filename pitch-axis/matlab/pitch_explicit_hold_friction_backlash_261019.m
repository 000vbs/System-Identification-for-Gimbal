%% Pitch 显式模型：位置保持力矩 + 静摩擦死区 + 回差/滞回
% 状态：
%   theta     关节角
%   w         角速度
%   tauA      执行器一阶滞后后的实际电机力矩
%   thetaH    保持力矩的有效位置（带回差）
%   stick     静摩擦粘滞状态
%
% 植物：
%   J*w' = tauA - kHold*(H(thetaH)-H(theta0))
%          - b0*w - Fc*tanh(w/wS)
%   |tauA - kHold*(H(thetaH)-H(theta0))| <= Fs 且 w≈0 时保持不动
%
% H(theta) 来自全行程静态标定表，先做单调化，避免原始测量凹坑造成负刚度。
clear; clc; close all;

here = fileparts(mfilename('fullpath'));
pitchRoot = fileparts(here);
OUT_DIR = fullfile(pitchRoot,'generated');
if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end
DATA_DIR = fullfile(pitchRoot,'data');
TABLE_CSV = fullfile(DATA_DIR,'pitch_hold_table_261019.csv');

S{1} = load_mode7(fullfile(DATA_DIR,'emberprobe-live-20261019-000638.csv'),'0.07 rad step (post-fix)');
S{2} = load_multi(fullfile(DATA_DIR,'emberprobe-live-20261019-210710.csv'),0,'0.1 rad step');
S{3} = load_multi(fullfile(DATA_DIR,'emberprobe-live-20261019-210710.csv'),1,'0.3 rad / 300 ms ramp');
S{4} = load_mode7(fullfile(DATA_DIR,'emberprobe-live-20261018-235835.csv'),'0.07 rad step (pre-fix)');

T = readtable(TABLE_CSV,'VariableNamingRule','preserve');
holdTheta = T.theta;
holdTau = T.tau_monotonic;

% J、Td 来自开环扫频；拟合 b0、动能摩擦、静摩擦额外量、平滑宽度、
% 回差宽度、保持力矩比例。
base.J = 0.07008;
base.Td = 0.01625;
base.tauAct = 0.00545;

x0s = [ ...
    0.78 0.94 0.50 0.05 0.002 1.00;
    0.20 1.50 0.50 0.05 0.005 1.00;
    1.00 2.20 0.30 0.03 0.010 0.90;
    0.50 1.00 1.50 0.08 0.002 1.10];
opt = optimset('Display','off','MaxFunEvals',1800,'MaxIter',1800, ...
    'TolFun',1e-9,'TolX',1e-8);

best=[]; bestCost=inf;
for i=1:size(x0s,1)
    [xb,~] = fminsearch(@(x) cost_fit(x,S(1:3),holdTheta,holdTau,base),x0s(i,:),opt);
    cb=cost_fit(xb,S(1:3),holdTheta,holdTau,base);
    fprintf('初值 %d: cost=%.5f b0=%.4f Fc=%.4f Fs=%.4f wS=%.4f bH=%.5f kH=%.3f\n', ...
        i,cb,xb(1),xb(2),xb(2)+xb(3),xb(4),xb(5),xb(6));
    if cb<bestCost, best=xb; bestCost=cb; end
end

b0=best(1); Fc=best(2); Fs=Fc+best(3); wS=best(4); bH=best(5); kHold=best(6);
fprintf('\n显式模型参数:\n');
fprintf('  J=%.6f kg*m^2\n',base.J);
fprintf('  b0=%.5f N*m*s/rad\n',b0);
fprintf('  Fc=%.5f N*m, Fs=%.5f N*m, wS=%.5f rad/s\n',Fc,Fs,wS);
fprintf('  backlash=%.6f rad (%.3f deg)\n',bH,bH*180/pi);
fprintf('  kHold=%.5f\n',kHold);
fprintf('  tauAct=%.5f s, Td=%.5f s\n',base.tauAct,base.Td);

%% 评估全部数据
n=numel(S);
metrics=zeros(n,5);
for k=1:n
    R=sim_case(S{k},best,holdTheta,holdTau,base);
    metrics(k,:)=[R.nrmseAngle,R.nrmseVel,R.nrmseTorque,R.peakErr,R.finalErr];
    fprintf('%s: angle %.1f%%, vel %.1f%%, torque %.1f%%, peakErr %.4f, finalErr %.4f\n', ...
        S{k}.name,100*metrics(k,1),100*metrics(k,2),100*metrics(k,3),metrics(k,4),metrics(k,5));
end

%% 三工况 3x3 图（0.07 post、0.1、0.3）
fig=figure('Color','w','Visible','off','Position',[30 30 1580 1100]);
for k=1:3
    R=sim_case(S{k},best,holdTheta,holdTau,base);
    subplot(3,3,k);
    plot(S{k}.t,S{k}.pos*180/pi,'o-','Color',[0.05 0.35 0.85],'LineWidth',1.3,'MarkerSize',4); hold on;
    plot(R.t,R.x*180/pi,'--','Color',[0.85 0.25 0.15],'LineWidth',2.0);
    plot(R.t,R.ref*180/pi,':','Color',[0.1 0.1 0.1],'LineWidth',1.0);
    grid on; ylabel('角度 (deg)'); title(S{k}.name);
    if k==1, legend('实测','显式模型','参考','Location','southeast'); end

    subplot(3,3,3+k);
    plot(S{k}.t,S{k}.vel,'o-','Color',[0.05 0.45 0.25],'LineWidth',1.3,'MarkerSize',4); hold on;
    plot(R.t,R.w,'--','Color',[0.85 0.25 0.15],'LineWidth',2.0);
    grid on; ylabel('角速度 (rad/s)');
    if k==1, legend('实测','显式模型','Location','northeast'); end

    subplot(3,3,6+k);
    plot(S{k}.t,S{k}.tau,'o-','Color',[0.55 0.20 0.65],'LineWidth',1.3,'MarkerSize',4); hold on;
    plot(R.t,R.tau,'--','Color',[0.85 0.25 0.15],'LineWidth',2.0);
    yline(0,'k:'); grid on; xlabel('阶跃后时间 (s)'); ylabel('力矩 (N·m)');
    if k==1, legend('实测','显式模型','Location','northeast'); end
end
sgtitle(sprintf(['Pitch 显式模型：保持表 + 静摩擦死区 + 回差/滞回\n' ...
    'b0=%.3f, Fc=%.3f, Fs=%.3f, backlash=%.3f deg, kHold=%.3f'], ...
    b0,Fc,Fs,bH*180/pi,kHold),'Interpreter','none');
exportgraphics(fig,fullfile(OUT_DIR,'pitch_explicit_hold_friction_backlash_261019.png'),'Resolution',160);
close(fig);

save(fullfile(OUT_DIR,'pitch_explicit_hold_friction_backlash_261019.mat'), ...
    'best','b0','Fc','Fs','wS','bH','kHold','base','metrics','S');
fprintf('\n图已保存: pitch_explicit_hold_friction_backlash_261019.png\n');

%% ---------- 数据提取 ----------
function S=load_mode7(csvPath,name)
    T=readtable(csvPath,'VariableNamingRule','preserve','TextType','string');
    tv=datetime(T.time,'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''','TimeZone','UTC');
    t=double(seconds(tv-tv(1)));
    st=numcol(T,'pitch_angle_step_state [u8]');
    q=find(st==2); i0=q(1); i1=q(end);
    S=make_seg(T,t,i0,i1,name,'Gimbal.pitch.target_pos [f32]');
end

function S=load_multi(csvPath,idxWanted,name)
    T=readtable(csvPath,'VariableNamingRule','preserve','TextType','string');
    tv=datetime(T.time,'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''','TimeZone','UTC');
    t=double(seconds(tv-tv(1)));
    st=numcol(T,'pitch_multi_state [u8]');
    idx=numcol(T,'pitch_multi_index [u8]');
    q=find(st==4 & idx==idxWanted); i0=q(1); i1=q(end);
    S=make_seg(T,t,i0,i1,name,'pitch_multi_target [f32]');
end

function S=make_seg(T,t,i0,i1,name,refCol)
    keep=i0:i1;
    S.name=name;
    S.t=t(keep)-t(i0);
    pos=numcol(T,'BPRS_log_pos_pitch [f32]');
    vel=numcol(T,'BPRS_log_w_pitch [f32]');
    tau=numcol(T,'BPRS_log_u_pitch [f32]');
    ref=numcol(T,refCol);
    S.p0=pos(i0);
    S.w0=vel(i0);
    S.pos=pos(keep)-S.p0;
    S.vel=vel(keep);
    S.tau=tau(keep);
    S.ref=ref(keep)-S.p0;
end

%% ---------- 仿真 ----------
function R=sim_case(S,x,holdTheta,holdTau,base)
    b0=x(1); Fc=x(2); Fs=Fc+x(3); wS=x(4); bH=x(5); kHold=x(6);
    Ts=1e-3;
    t=(0:Ts:S.t(end))';
    N=numel(t);
    ref=interp1(S.t,S.ref,t,'linear','extrap');
    J=base.J; Td=base.Td; tauAct=base.tauAct;
    Kpa=15; Kpw=3; Kiw=20;
    tauFF=0.94; wLin=0.30; kFF=25; ffMax=2.0;
    pidMax=3.5; iMax=2.0; angVelMax=3.0;

    xAng=0; w=S.w0; tauA=0; Iout=0; thetaH=0; stick=true;
    Ndelay=max(1,round(Td/Ts)); delay=zeros(Ndelay,1);
    xLog=zeros(N,1); wLog=zeros(N,1); tauLog=zeros(N,1); stickLog=false(N,1);
    h0=hold_eval(S.p0,holdTheta,holdTau);

    for n=1:N
        tv=Kpa*(ref(n)-xAng);
        tv=max(-angVelMax,min(angVelMax,tv));
        ev=tv-w;
        Iout=Iout+Kiw*Ts*ev;
        Iout=max(-iMax,min(iMax,Iout));
        pidOut=Kpw*ev+Iout;
        pidOut=max(-pidMax,min(pidMax,pidOut));
        s=max(-1,min(1,tv/wLin));
        springFF=max(-ffMax,min(ffMax,kFF*ref(n)));
        raw=pidOut+tauFF*s+springFF;

        tauDelay=delay(1);
        delay=[delay(2:end);raw];
        tauA=tauA+(Ts/tauAct)*(tauDelay-tauA);

        h=hold_eval(S.p0+thetaH,holdTheta,holdTau)-h0;
        tauHold=kHold*h;
        tauDrive=tauA-tauHold;

        if stick
            if abs(tauDrive)<=Fs
                a=0; wNew=0;
            else
                stick=false;
                a=(tauDrive-b0*w-Fc*tanh(w/wS))/J;
                wNew=w+a*Ts;
            end
        else
            a=(tauDrive-b0*w-Fc*tanh(w/wS))/J;
            wNew=w+a*Ts;
        end

        xNew=xAng+0.5*(w+wNew)*Ts;

        if bH>0
            if xNew>thetaH+bH, thetaH=xNew-bH;
            elseif xNew<thetaH-bH, thetaH=xNew+bH;
            end
        else
            thetaH=xNew;
        end

        if ~stick && w*wNew<=0
            hNew=hold_eval(S.p0+thetaH,holdTheta,holdTau)-h0;
            if abs(tauA-kHold*hNew)<=Fs
                wNew=0; stick=true;
            end
        end

        xAng=xNew; w=wNew;
        xLog(n)=xAng; wLog(n)=w; tauLog(n)=raw; stickLog(n)=stick;
    end

    xm=interp1(t,xLog,S.t,'linear','extrap');
    wm=interp1(t,wLog,S.t,'linear','extrap');
    tm=interp1(t,tauLog,S.t,'linear','extrap');
    R.t=t; R.x=xLog; R.w=wLog; R.tau=tauLog; R.ref=ref; R.stick=stickLog;
    R.nrmseAngle=sqrt(mean((xm-S.pos).^2))/max(max(S.pos)-min(S.pos),0.02);
    R.nrmseVel=sqrt(mean((wm-S.vel).^2))/max(max(abs(S.vel)),0.2);
    R.nrmseTorque=sqrt(mean((tm-S.tau).^2))/max(max(abs(S.tau)),0.5);
    R.peakErr=max(abs(xm-S.pos));
    R.finalErr=xm(end)-S.pos(end);
end

function h=hold_eval(theta,holdTheta,holdTau)
    th=min(max(theta,holdTheta(1)),holdTheta(end));
    h=interp1(holdTheta,holdTau,th,'pchip');
end

function c=cost_fit(x,S,holdTheta,holdTau,base)
    b0=x(1); Fc=x(2); Fs=Fc+x(3); wS=x(4); bH=x(5); kHold=x(6);
    if b0<0 || b0>4 || Fc<0.02 || Fc>5 || x(3)<0 || Fs>7 || ...
            wS<0.002 || wS>0.25 || bH<0 || bH>0.04 || kHold<0.2 || kHold>2.0
        c=1e8; return;
    end
    c=0;
    for k=1:numel(S)
        R=sim_case(S{k},x,holdTheta,holdTau,base);
        c=c+R.nrmseAngle^2+0.2*R.nrmseVel^2+0.2*R.nrmseTorque^2;
    end
    c=c/numel(S);
end

function x=numcol(T,name)
    v=T.(name);
    if isnumeric(v), x=double(v); else, x=str2double(string(v)); end
end
