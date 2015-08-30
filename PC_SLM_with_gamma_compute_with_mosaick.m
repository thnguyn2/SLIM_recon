%This code is for obtainining the growth curve measurement
clc;
clear all;
close all;
x1 = 876;
x2 = 1932;
y1 = 1723;
y2 = 2779;
nrows = y2-y1+1;
ncols = x2-x1+1;
[xx,yy]=meshgrid(1:ncols,1:nrows);

%Location for drawing the cross section - Do not use for most of the time
x0=1;
y0=1;
half_cs_length=50;
%----------------------------------------


 
gpu_compute_en =1; %1-Enable GPU computing
  
fpath ='F:\SLIM_data_for_the_coherence_paper\Data for growth curve measurement\20x\';
addpath(fpath);
recpathname='F:\SLIM_data_for_the_coherence_paper\Data for growth curve measurement\rec_data\';
addpath(recpathname);
%Create a background mask
bg_thresh = 0.04;
%Filter the image to avoid very large value of the noise
h_denoise = fspecial('gaussian',[9 9],0.25);
 
lambda_weight =10;
beta_weight=0;
tol = 1e-4; %We don't need to find the best in each step since we will tweak 2 variables t and g at the same time
niter=100;
bw_array = linspace(23,23,1);
nbw = length(bw_array);
ntimestep = 70;
for bwidx=1:nbw %Go through each possible bandwidth
    for timestepidx = 1:ntimestep %So through each time stamp
     disp(['Working at time step ' num2str(timestepidx)]);
        bw =bw_array(bwidx);%This is the bandwidth parameter of the correlation ufnciton
        disp(['Working at scale' num2str(bw)]);
        if (~exist(strcat(recpathname,num2str(bw))))
            mkdir(recpathname,num2str(bw));  
        end
        if (~exist(strcat(recpathname,strcat(num2str(bw),'\org\'))))
            mkdir(recpathname,strcat(num2str(bw),'\org\'));  
        end
        h=fspecial('gaussian',[round(6*bw)+1 round(6*bw)+1],bw); %Transfer function of the low-pass filter...
        h1 = zeros(nrows,ncols);
        h1(1:size(h,1),1:size(h,2))=h;
        kernel_size=size(h,1);
        h1 = circshift(h1,[-round((kernel_size-1)/2) -round((kernel_size-1)/2)]);

         curfilename=strcat(fpath,num2str(timestepidx-1),'.tif');
        uphi = imread(curfilename,'PixelRegion',{[y1 y2],[x1 x2]}); %Read the current file into the memory
        uphi = cast(uphi,'double');
        figure(1);
        subplot(221);
        imagesc(uphi);colorbar;colormap jet
        title('SLIM phase image')
        uphi_denoised = imfilter(uphi,h_denoise,'same');
        idx = find(abs(uphi_denoised)>bg_thresh);
        mask = zeros(size(uphi));
        mask(idx)=1;
        mask = medfilt2(mask);
        mask = imclose(mask,strel('disk',5));
        subplot(222);
        imagesc(mask);
        libimphase = uphi_denoised;
        libimphase(idx)=NaN;
        phasebg = inpaint_nans(libimphase,2);
        subplot(223);
        imagesc(phasebg);colorbar;
        gamma = exp(i*(uphi_denoised-phasebg));
        if (gpuDeviceCount()==0)
                gpu_compute_en = 0; %if there is no gpu, compute the result on cpu instead
        end
        figure(2);
        imagesc(mask);
        method = 'relax';
        if (gpu_compute_en==0)
            [gk,tk] = estimate_gt(gamma,h,niter,lambda_weight,beta_weight,tol,method);
        else %Compute gk and tk on gpu
            d = gpuDevice();
            reset(d); %Reset the device and clear its memmory
            [gk,tk] = estimate_gt_gpu(gamma,h,niter,lambda_weight,beta_weight,tol,x0,y0,half_cs_length,method,mask);
        end 
        figure(5);
        imagesc(angle(tk));colorbar; axis off;
        writeTIFF(unwrap2(cast(angle(tk),'double')),strcat(recpathname,num2str(bw),'\x1',num2str(x1),'_x2',num2str(x2),...
            '_y1',num2str(y1),'_y2',num2str(y2),'_',num2str(timestepidx-1),'.tif'));
         writeTIFF(uphi_denoised,strcat(recpathname,num2str(bw),'\org\org_x1',num2str(x1),'_x2',num2str(x2),...
            '_y1',num2str(y1),'_y2',num2str(y2),'_',num2str(timestepidx-1),'.tif'));

    end
end
   
  
  
 