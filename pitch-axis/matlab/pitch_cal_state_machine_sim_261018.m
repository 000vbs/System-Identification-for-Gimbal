%% Pitch 静态标定状态机离线复现（与固件 Gimbal_Pitch_Static_Cal 同逻辑）
% 目的：先确认 0/+0.03/+0.06/... 序列在 0.94N·m 摩擦下能采到可用点，
%       以及固件内 tau=k*theta+b 拟合能不能恢复 k。
clear; clc; close all;

% 植物
J=0.07008; kTrue=15.783; b0=0.7784; tauC=0.9391; Td=0.01625; wF=0.03;
Ts=1e-3;

% 状态机参数（与固件默认一致）
amp=0.06; step=0.03; rampMs=300; settleMs=1000; avgMs=300;
tauLim=3.5; excLim=0.15; minR2=0.90;

% 控制器
Kpa=11.77; Kpw=1.00; Kiw=6.0; KiCode=Kiw*Ts;

offsets=[0 step amp step 0 -step -amp -step 0];
N=9;
theta0=0; th=0; w=0;
Iout=0; delay=zeros(round(Td/Ts),1); tauApply=0;

% 状态：1=move 2=settle 3=average
state=2; point=0; phaseT=0; prev=theta0; goal=theta0; target=theta0;
sumTh=0; sumTau=0; sumN=0;
thetaPts=zeros(1,N); tauPts=zeros(1,N); valid=false(1,N); stuck=false(1,N);
tLog=[]; thLog=[]; tauLog=[]; targetLog=[];

maxT=20; nMax=round(maxT/Ts);
done=false;
for n=1:nMax
    t=(n-1)*Ts;

    if state==1
        frac=min(1,(t-phaseT)/ (rampMs/1000));
        target=prev+(goal-prev)*frac;
        if frac>=1, state=2; phaseT=t; end
    elseif state==2
        target=goal;
        if t-phaseT >= settleMs/1000
            state=3; phaseT=t; sumTh=0; sumTau=0; sumN=0;
        end
    elseif state==3
        target=goal;
        if t-phaseT >= avgMs/1000 && sumN>0
            thA=sumTh/sumN; tauA=sumTau/sumN;
            thetaPts(point+1)=thA; tauPts(point+1)=tauA;
            if abs(thA-goal)<0.01, valid(point+1)=true; else, stuck(point+1)=true; end
            point=point+1;
            if point>=N, done=true; break; end
            prev=goal; goal=theta0+offsets(point+1);
            phaseT=t; state=1; sumTh=0; sumTau=0; sumN=0;
        end
    end

    % 控制器
    tv=Kpa*(target-th);
    tv=max(-3,min(3,tv));
    ev=tv-w;
    Iout=Iout+KiCode*ev;
    Iout=max(-2,min(2,Iout));
    raw=Kpw*ev+Iout;
    s=max(-1,min(1,tv/0.30));
    raw=raw+tauC*s;
    raw=max(-tauLim,min(tauLim,raw));

    tauApply=delay(1);
    delay=[delay(2:end); raw];

    sub=10; h=Ts/sub; x=[th;w];
    for q=1:sub
        k1=fp(x,tauApply,J,kTrue,b0,tauC,wF);
        k2=fp(x+0.5*h*k1,tauApply,J,kTrue,b0,tauC,wF);
        k3=fp(x+0.5*h*k2,tauApply,J,kTrue,b0,tauC,wF);
        k4=fp(x+h*k3,tauApply,J,kTrue,b0,tauC,wF);
        x=x+(h/6)*(k1+2*k2+2*k3+k4);
    end
    th=x(1); w=x(2);

    if state==3
        sumTh=sumTh+th; sumTau=sumTau+raw; sumN=sumN+1;
    end
    tLog(end+1)=t; thLog(end+1)=th; tauLog(end+1)=raw; targetLog(end+1)=target; %#ok<SAGROW>
end

% 固件内同样的拟合
n=sum(valid);
if n<3
    fprintf('有效点不足: n=%d, stuck=%s\n',n,mat2str(find(stuck)));
else
    mth=mean(thetaPts(valid)); mtau=mean(tauPts(valid));
    sxx=sum((thetaPts(valid)-mth).^2);
    sxy=sum((thetaPts(valid)-mth).*(tauPts(valid)-mtau));
    kFit=sxy/sxx; bFit=mtau-kFit*mth; cFit=-bFit/kFit;
    fit=kFit*thetaPts(valid)+bFit;
    r2=1-sum((tauPts(valid)-fit).^2)/sum((tauPts(valid)-mtau).^2);
    fprintf('\n标定结果: k=%.3f (真值 %.3f), center=%.4f, R2=%.4f, 有效点=%d\n', ...
        kFit,kTrue,cFit,r2,n);
end
fprintf('点序号 目标     实际角   力矩    有效\n');
for i=1:N
    fprintf('  %d   %+7.4f  %+7.4f  %+7.3f   %d\n', ...
        i-1,theta0+offsets(i),thetaPts(i),tauPts(i),valid(i));
end

figure('Color','w','Position',[100 100 1100 700]);
subplot(2,1,1);
plot(tLog,targetLog*180/pi,'k--'); hold on; plot(tLog,thLog*180/pi,'b');
grid on; xlabel('时间 (s)'); ylabel('角度 (deg)'); legend('目标','实际');
title('静态标定状态机：目标与实际角');
subplot(2,1,2);
plot(tLog,tauLog,'b'); grid on; xlabel('时间 (s)'); ylabel('力矩 (N·m)');
title('对应的保持力矩');

figure('Color','w');
plot(thetaPts*180/pi,tauPts,'ko','MarkerFaceColor','k'); hold on;
if n>=3
    thg=linspace(min(thetaPts(valid)),max(thetaPts(valid)),100);
    plot(thg*180/pi,kFit*thg+bFit,'b-','LineWidth',1.5);
end
grid on; xlabel('角度 (deg)'); ylabel('力矩 (N·m)');
title('九个静态点与拟合');

function dx=fp(x,tau,J,k,b0,tauC,wF)
    th=x(1); w=x(2);
    dx=[w;(tau-b0*w-tauC*tanh(w/wF)-k*th)/J];
end
