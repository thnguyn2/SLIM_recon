clc;
clear all;
close all;
for cellidx=1:1 %Go through all cells

    %fname=strcat('F:\SLIM_data_for_the_coherence_paper\BestShots\cells_',num2str(cellidx),'\');
    fname ='F:\SLIM_data_for_the_coherence_paper\BestShots\pillar_20\'
    bgpathname='F:\SLIM_data_for_the_coherence_paper\BestShots\cells_background\'
    addpath(bgpathname);

    %I_back=im2double(imread('1859.tif'));
    I_back=im2double(imread('bg3.tif'));

    [nrows,ncols]=size(I_back);
    % I_back=0;
    %I_back=mean(mean(I_back));
    num=1855;
    mm=num2str(num);
    xx1=250; xx2=255; %1-1040 - An area with background only
    yy1=100; yy2=105; %1-1388

    fname_save=strcat(fname,mm,'.jpg'); %jpg
    minPhase=-0.5; %min(min(uphi)); -0.5
    maxPhase=4; %max(max(phase));  
    usephaseinfoonly =1;%A flag that specifies whether we should use the phase only or both the amplitude and the phase
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % A=im2double((imread('01563_p000001t00000001z077c04.tif')));
    % B=im2double((imread('01563_p000001t00000001z077c01.tif')));
    % C=im2double((imread('01563_p000001t00000001z077c02.tif')));
    % D=im2double((imread('01563_p000001t00000001z077c03.tif')));
    A=im2double((imread(strcat(fname,'image_1_0_D.tif'))));
    B=im2double((imread(strcat(fname,'image_1_0_A.tif'))));
    C=im2double((imread(strcat(fname,'image_1_0_B.tif'))));
    D=im2double((imread(strcat(fname,'image_1_0_C.tif'))));

    [uphi,beta,del_phi,u0sqr,gamma]=combine_4_frame(A,B,C,D,I_back,2.3,...
        xx1,xx2,yy1,yy2);
   
    %%For trial only... beta = real(sqrt(x1./x2))*2.3;
    figure(1)
    imagesc(beta);title('Beta');

    %uphi=imsmooth(uphi,3);
    %uphi=uphi-bg;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(2)
    imagesc(uphi);colorbar
    axis image
    axis off
    writeTIFF(uphi,'phase_rec.tif')

    %%Added section for recovery....
    [nrows,ncols]=size(uphi);
    gamma_amp = ones(nrows,ncols);
    %Create a background mask
    bg_thresh = 0.11;
    %Filter the image to avoid very large value of the noise
    h_denoise = fspecial('gaussian',[9 9],2);
    uphi_denoised = imfilter(uphi,h_denoise,'same');
    idx = find(abs(uphi_denoised)>bg_thresh);
    mask = zeros(size(uphi));
    mask(idx)=1;
    mask = medfilt2(mask);
    mask = imclose(mask,strel('disk',5));
    %Fill out the region, assume to object is solid (true for pillar)
    mask = imfill(mask,'hole');
    premask = mask; %A mask with all sample defects available
%     mask = bwareaopen(mask,4000); %Filter out all the dirt and scratches on the surface
%     maskdiff = premask-mask;
%     defectidx = find(maskdiff==1);
%     uphi_denoised(defectidx)=0;
    [xx,yy]=meshgrid(1:ncols,1:nrows);
%    gamma = abs(gamma).*exp(i*uphi_denoised);
    libimphase = uphi_denoised;
    libimamp = abs(gamma); %Do not fix for the amplitude
    libimphase(idx)=NaN;
    libimamp(idx)=NaN;
    
    
    ampbg = inpaint_nans(libimamp,2);
    phasebg = inpaint_nans(libimphase,2);

    gamma = gamma./(ampbg.*exp(i*phasebg));

    
    if (usephaseinfoonly)%just use the phase instead of the amplitude
        gamma = gamma./abs(gamma);
    end

    x0=1;
    y0=1;
    half_cs_length=50;
    inverse = 1;
    bw_array = linspace(14,14,1);
    %bw_array = 35;
    nbw = length(bw_array);
        for bw_arrayidx = 1:nbw
            if (inverse)
                bw =round(bw_array(bw_arrayidx));%This is the bandwidth parameter of the correlation ufnciton
                disp(['Current bw: ' num2str(bw)]);
                %For old SLIM 40x
                % bw = 32 with amplitude
                % bw = 30 without amplitude
                h=fspecial('gaussian',[round(6*bw)+1 round(6*bw)+1],bw); %Transfer function of the low-pass filter...

                h1 = zeros(nrows,ncols);
                h1(1:size(h,1),1:size(h,2))=h;
                kernel_size=size(h,1);
                h1 = circshift(h1,[-round((kernel_size-1)/2) -round((kernel_size-1)/2)]);
                gpu_compute_en =1; %1-Enable GPU computing
                %First, initialize tk and lk. Here, gk = t v h;
                lambda_weight =10;
                beta_weight=0;
                tol = 1e-4; %We don't need to find the best in each step since we will tweak 2 variables t and g at the same time
                niter=100;
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
                 bgsubtract_str='_no_bgs_phase_only';
                 writeTIFF(unwrap2(cast(phasemap,'double')),strcat(fname,'phase_t_edge_gaussian_',num2str(bw),bgsubtract_str,'_.tif'));
                 writeTIFF(abs(tk),strcat(fname,'amp_t_edge_gaussian_',num2str(bw),bgsubtract_str,'.tif'));
                 writeTIFF(angle(gamma),strcat(fname,'org_phase_gaussian_',num2str(bw),bgsubtract_str,'.tif'));
            end
        end
end
