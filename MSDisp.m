%Main function of MSD calculation
%%
%initialization
DC=TrRead('J:\S2CellGrowth\30thOctober2010\Analysis-Max\1-cellAnalysis\RoiSet(1)11.txt'); %DC=Data Cell
%DC=TrRead('2635.txt'); %DC=Data Cell
%DC=TrRead('1umglycerol.txt'); %DC=Data Cell
NumOfFrames=110; %Number of Frames
dT=10; %time interval
NumOfParticles=size(DC,2);
cf=1/3.5; %conversion factor from pixel to um, 
%%
MSDx=zeros(1,NumOfFrames);
MSDy=MSDx;

for m=1:NumOfFrames
    N=0;
    for n=1:NumOfParticles
        for p=1:size(DC{n},1)
        X0=DC{n}(p,2);
        Y0=DC{n}(p,3);
        idx=(DC{n}(:,1)-DC{n}(p,1)==m-1);
        if sum(idx)
            Xm=DC{n}(idx,2);%2 x
            Ym=DC{n}(idx,3);%3 y
            MSDx(m)=MSDx(m)+(Xm-X0)^2;
            MSDy(m)=MSDy(m)+(Ym-Y0)^2;
            N=N+1;
        end
        end
    end
    if N~=0
       MSDx(m)=MSDx(m)/N;
       MSDy(m)=MSDy(m)/N;
       M(m)=N;
    end
end
MSD=MSDx+MSDy;
%%
MSDx=MSDx*cf^2;
MSDy=MSDy*cf^2;
MSD=MSD*cf^2;
%%
TimeLag = 0:dT:dT*(NumOfFrames-1);
plot(0:dT:dT*(NumOfFrames-1),MSD,0:dT:dT*(NumOfFrames-1), MSDx, 0:dT:dT*(NumOfFrames-1), MSDy)
xlabel('Time (second)')
ylabel('MSD (\mum^2)')
figure
loglog(0:dT:dT*(NumOfFrames-1),MSD,0:dT:dT*(NumOfFrames-1), MSDx, 0:dT:dT*(NumOfFrames-1), MSDy);
%%
% X=datacell{:}(:,2);
% Y=report(:,3);
% m0=report(:,4);%intensity moment 0
% m2=report(:,5);%intensity moment 2
% %%
% for n=1:size(X,1)-1
%      MSD(n)=mean((X(n+1:n:end)-X(1:n:end-n)).^2+(Y(n+1:n:end)-Y(n+1:n:end)).^2);
% end
% %MSD=(X-X(1)).^2+(Y-Y(1)).^2;
% %%

%%
% dim=2;
% for n=3:3
%     aa=DC{n};
%     plot(aa(:,1),aa(:,2)-mean(aa(:,2)),aa(:,1),aa(:,3)-mean(aa(:,3)))
%     axis([0 255 -3 3])
% end
% n=0;
% aa(:,dim)=-aa(:,dim);
% for m=1:234
%     for qq=1:(aa(m+1,1)-aa(m,1))
% plot(aa(1:m,1)*dT,aa(1:m,dim)-mean(aa(:,dim))); 
% %axis off
% axis([0 255*dT -3 3])
% xlabel('Time (sec)')
% ylabel('Displacement (\mum)')
% %print('-djpeg', '-r67',strcat('L:\cell dispersion analysis\beads analysis\tracking\',num2str(aa(m)+qq),'.jpg')); %-djpeg
% print('-djpeg', '-r67',strcat('M:\Fluorescence Beads\experimental data\1642\',num2str(aa(m)+qq),'.jpg')); %-djpeg
%     end
% end
% plot(aa(1:m+1,1)*dT,aa(1:m+1,dim)-mean(aa(:,dim))); 
% %axis off
% axis([0 255*dT -3 3])
% xlabel('Time (sec)')
% ylabel('Displacement (\mum)')
% %print('-djpeg', '-r67',strcat('L:\cell dispersion analysis\beads analysis\tracking\',num2str(aa(m+1)+qq),'.jpg')); %-djpeg
% print('-djpeg', '-r67',strcat('M:\Fluorescence Beads\experimental data\1642\',num2str(aa(m+1)+qq),'.jpg')); %-djpeg
