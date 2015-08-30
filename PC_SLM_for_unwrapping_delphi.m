%This program is for testing whether unwrapping delta phi first is better
%than doing it directly on phi
clc;
clear all;
close all;
imdir = 'Data_with_phase_wrapping_during_mitosis\'
im_prefix = '51_0_1_0_';
%fname='G:\SLIM\test\35msbg\';
%fname = 'D:\Data_from_catalin\74 134 189 246\';
fname = '';
addpath(fname);
% A=im2double((imread('0.tif')));
% B=im2double((imread('1.tif')));
% C=im2double((imread('2.tif')));
% D=im2double((imread('3.tif')));


A=im2double((imread(strcat(imdir,im_prefix,'4.tif'))));
B=im2double((imread(strcat(imdir,im_prefix,'1.tif'))));
C=im2double((imread(strcat(imdir,im_prefix,'2.tif'))));
D=im2double((imread(strcat(imdir,im_prefix,'3.tif'))));

bgsubtraction =0;
if (bgsubtraction)
    I_back=im2double(imread('1859.tif'));
    A=A-I_back; %A->I0 0pi
    %B=RawRead2(fname,num+0);
    B=B-I_back; %B->I1
    %C=RawRead2(fname,num+1);
    C=C-I_back; %C->I2
    %D=RawRead2(fname,num+2);
    D=D-I_back; %D->I3
end

[nrows,ncols]=size(A);
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
beta=LL*2.56; %real beta, 2.3 for 40X and 4 for 10X
%beta=x1*2.5;

%%For trial only... beta = real(sqrt(x1./x2))*2.3;
figure
imagesc(beta);

phi=atan2(beta.*sin(del_phi),1+beta.*cos(del_phi));
%unwrap2
%uphi=unwrap2(phi);
uphi = phi;
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
writeTIFF(uphi,strcat(fname,'phase_rec.tif'));

%%Added section for recovery....
[nrows,ncols]=size(uphi);
gamma_amp = ones(nrows,ncols);
%Create a background mask
bg_thresh = 0.03;
%Filter the image to avoid very large value of the noise
h_denoise = fspecial('gaussian',[9 9],2);
uphi_denoised = imfilter(uphi,h_denoise,'same');
idx = find(abs(uphi_denoised)>bg_thresh);
mask = zeros(size(uphi));
mask(idx)=1;
[xx,yy]=meshgrid(1:ncols,1:nrows);
libim = uphi_denoised;
libim(idx)=NaN;
uphibg = inpaint_nans(libim,2);
uphi = uphi - uphibg;
gamma = gamma_amp.*exp(i*uphi);
figure
imagesc(uphi);


%Do interpolation for the remaining images of the object region to get a
%nice background


x0=1;
y0=1;
half_cs_length=50;
inverse = 1;
if (inverse)
    bw = 30;%This is the bandwidth parameter of the correlation ufnciton
    h=fspecial('gaussian',[round(6*bw)+1 round(6*bw)+1],bw); %Transfer function of the low-pass filter...

    h1 = zeros(nrows,ncols);
    h1(1:size(h,1),1:size(h,2))=h;
    kernel_size=size(h,1);
    h1 = circshift(h1,[-round((kernel_size-1)/2) -round((kernel_size-1)/2)]);
    gpu_compute_en =0; %1-Enable GPU computing
    %First, initialize tk and lk. Here, gk = t v h;
    lambda_weight =1;
    beta_weight=0;
    tol = 1e-4; %We don't need to find the best in each step since we will tweak 2 variables t and g at the same time
    niter=1800;
    if (gpuDeviceCount()==0)
          gpu_compute_en = 0; %if there is no gpu, compute the result on cpu instead
    end

    maxbg_phase = 0.5;
     method = 'relax';
     if (gpu_compute_en==0)
        [gk,tk] = estimate_gt(gamma,h,niter,lambda_weight,beta_weight,tol,method);
     else %Compute gk and tk on gpu
        d = gpuDevice();
        reset(d); %Reset the device and clear its memmory
        [gk,tk] = estimate_gt_gpu(gamma,h,niter,lambda_weight,beta_weight,tol,x0,y0,half_cs_length,method,mask);
     end

     figure(5);
     subplot(211);imagesc(angle(tk));colorbar; axis off;
     title('Phase after solving');
     subplot(212);imagesc(abs(tk));colorbar; axis off;
     title('Absorbance after solving');
     phasemap =angle(tk)-mean(mean(angle(tk(1:200,1:200)))); 

     writeTIFF(unwrap2(cast(phasemap,'double')),strcat('phase_t_edge.tif'));
     writeTIFF(abs(tk),strcat('amp_t_edge.tif'));
     writeTIFF(angle(gamma),strcat('org_phase.tif'));
end



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