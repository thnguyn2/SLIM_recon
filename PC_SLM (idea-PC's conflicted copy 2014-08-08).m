clc;
clear all;
close all;
fname='G:\SLIM\test\35msbg\';
I_back=im2double(imread('1859.tif'));
[nrows,ncols]=size(I_back);
% I_back=0;
%I_back=mean(mean(I_back));
num=1855;
mm=num2str(num);
%xx1=250; xx2=255; %1-1040 - An area with background only
%yy1=100; yy2=105; %1-1388
xx1 = 1; xx2 = nrows;
yy1 = 1; yy2 = ncols;
fname_save=strcat(fname,mm,'.jpg'); %jpg
minPhase=-0.5; %min(min(uphi)); -0.5
maxPhase=4; %max(max(phase));  4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A=RawRead2(fname,num+0); A=A-I_back; %A->I0 0pi
% B=RawRead2(fname,num+3); B=B-I_back; %B->I1
% C=RawRead2(fname,num+2); C=C-I_back; %C->I2
% D=RawRead2(fname,num+1); D=D-I_back; %D->I3

A=im2double((imread('01563_p000001t00000001z077c04.tif')));
%A=RawRead2(fname,num+3);
A=A-I_back; %A->I0 0pi
%B=RawRead2(fname,num+0);
B=im2double((imread('01563_p000001t00000001z077c01.tif')));
B=B-I_back; %B->I1
%C=RawRead2(fname,num+1);
C=im2double((imread('01563_p000001t00000001z077c02.tif')));
C=C-I_back; %C->I2
%D=RawRead2(fname,num+2);
D=im2double((imread('01563_p000001t00000001z077c03.tif')));
D=D-I_back; %D->I3

Gs=D-B; Gc=A-C; del_phi=atan2(Gs,Gc);

%beta
L=(A-C+D-B)./(sin(del_phi)+cos(del_phi))/4; %E0*E1
g1=(A+C)/2; %E0^2+E1^2=s
g2=L.*L; %E0^2*E1^2=p
x1=g1/2-sqrt(g1.*g1-4*g2)/2; x2=g1/2+sqrt(g1.*g1-4*g2)/2; %solutions
beta1=sqrt(x1./x2); beta2=1./beta1;
figure
imagesc(del_phi);
%beta1=x1;beta2=x2;
%get constant from  average over pixels
cL=L(xx1:xx2,yy1:yy2); cbeta1=beta1(xx1:xx2,yy1:yy2);
L1=real(mean2(cbeta1))/mean2(cL)    %[Tan]: this is 1/<Uo>^2
LL=L1*L;
beta=LL*2.3; %real beta, 2.3 for 40X and 4 for 10X
%beta=x1*2.5;

%%For trial only... beta = real(sqrt(x1./x2))*2.3;
figure
imagesc(beta);

phi=atan2(beta.*sin(del_phi),1+beta.*cos(del_phi));
%unwrap2
uphi=unwrap2(phi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uphi=medfilt2(uphi);
%uphi=imsmooth(uphi,3);
%uphi=uphi-bg;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
height=uphi*0.552/2/pi/1*1000; %heigth=height+5;
imagesc(uphi,[-0.4,1.5]);colorbar
%imagesc(height,[-15,15]);colorbar
axis image
axis off
%colormap(gray)
%print('-djpeg','-r300', fname_save); %-djpeg
writeTIFF(uphi,'phase_rec.tif')
% %print -dtiff -r300 fname_save
% phase=uphi-minPhase;
% uniPhase=phase/maxPhase;
% mPhase=num2str(maxPhase);
% imwrite(uint16((2^16-1)*uniPhase),strcat(fname,mm,'rec_',mPhase,'_.tif'),'tiff','Compression','none');
% %imwrite(uint16((2^16-1)*height/10),strcat(fname,mm,'rec_height_',num2str(10),'_.tif'),'tiff','Compression','none');
% imwrite(uint16((2^16-1)*(del_phi+pi)/2/pi),strcat(fname,mm,'rec_beta_',num2str(6.28),'_.tif'),'tiff','Compression','none');
% T_cos=cos(uphi);
% T_cos=(-T_cos+1)/2;
% imwrite(uint16((2^16-1)*(T_cos)),strcat(fname,mm,'rec_cos_',mm,'.tif'),'tiff','Compression','none');