clc;
clear all;
close all;
bgpath = '.\Data\mikhail\cells_background\';
addpath(bgpath);
dirpath = '.\Data\mikhail\cells_3\'; %path the files
addpath(dirpath);
imname='acq3.tif';
bgname='bg3.tif';
gamma_phase = im2double(imread(strcat(dirpath,imname)));
bg_phase = im2double(imread(strcat(bgpath,bgname)));
%For the image of Mikhail's
gamma_phase = gamma_phase(250:1250,450:1500);
bg_phase = bg_phase(250:1250,450:1500);

%Generate a noisy map and see if it affects the stability of the system
noisestd = 0.05;
%gamma_phase = noisestd*randn(size(gamma_phase));
%gamma_phase = bg_phase;

%gamma_phase = gamma_phase-bg_phase;
[nrows,ncols]=size(gamma_phase);
gamma_amp = ones(nrows,ncols);
%Find all NAN pixels
nanim = isnan(gamma_phase);
nanidx = find(nanim==1);
if (~isempty(nanidx))
    new_gamma_phase = gamma_phase;
    new_gamma_phase(nanidx)=0;
    new_gamma_phase = medfilt2(new_gamma_phase);
    gamma_phase = new_gamma_phase;
end
figure(1)
subplot(121);
imagesc(gamma_phase);
colorbar;
title('Original phase image');

% %First, generate a mask for the background region\
maxbg_phase = 0.02; %This number is small is better. Too small will not have much effects
% %Generating a mask that denotes the location with object information
mask = zeros(nrows,ncols);
idx = find(abs(gamma_phase)>maxbg_phase);
mask(idx)=0;
% [xx,yy]=meshgrid(1:ncols,1:nrows);

% xfit = xx(idx); xfit = xfit(:);
% yfit = yy(idx); yfit = yfit(:);
% zfit = gamma_phase(idx);
% fitobj = fit([xfit yfit],zfit, 'poly22'); %Fit with polynomial 1 in both x and y
% fitbg = feval(fitobj,[xx(:) yy(:)]);
% fitbg=reshape(fitbg,[nrows ncols]);
%avg_filter = fspecial('gaussian',[40 40],5);
%fitbg = imfilter(bg_phase,avg_filter,'same');
%gamma_phase=fitbg;

%figure(1)
%subplot(122);
%imagesc(fitbg);title('Fitted background to correct for abberation');
gamma = gamma_amp.*exp(i*gamma_phase);



x0=1;
y0=1;
half_cs_length=50;
inverse = 1;

%load(strcat('.\Data\mikhail\pillar_20\','h_edge_915_to_957_855_to_969.mat'));
%h = h_mask;

bw = 20;%This is the bandwidth parameter of the correlation ufnciton
h=fspecial('gaussian',[round(4*bw)+1 round(4*bw)+1],bw); %Transfer function of the low-pass filter...

h1 = zeros(nrows,ncols);
h1(1:size(h,1),1:size(h,2))=h;
kernel_size=size(h,1);
h1 = circshift(h1,[-round((kernel_size-1)/2) -round((kernel_size-1)/2)]);
gpu_compute_en =1; %1-Enable GPU computing
 
%First, initialize tk and lk. Here, gk = t v h;
lambda_weight =1;
beta_weight=0;
tol = 1e-4; %We don't need to find the best in each step since we will tweak 2 variables t and g at the same time
niter=500;

%Downsample image and the filter by a factor of 3 to speed up
%h=h(1:3:end,1:3:end);

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
 writeTIFF(unwrap2(cast(phasemap,'double')),strcat(dirpath,'phase_t_edge.tif'));
 writeTIFF(abs(tk),strcat(dirpath,'amp_t_edge.tif'));
 writeTIFF(angle(gamma),strcat(dirpath,'org_phase.tif'));
