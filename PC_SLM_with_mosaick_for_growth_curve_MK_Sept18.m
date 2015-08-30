%This code is for obtainining the growth curve measurement
clc;
clear all;
close all;

nfovs=1;
ntimesteps = 1;


%Location for drawing the cross section - Do not use for most of the time
x0=1;
y0=1;
half_cs_length=50;
%----------------------------------------


fpath ='E:\Data for SLIM edge of the Halo removal paper\80 micron measurement\40xPh20.75\';
addpath(fpath);

gpu_compute_en =1; %1-Enable GPU computing
%Create a background mask
bg_thresh = 0.02;
%Filter the image to avoid very large value of the noise
h_denoise = fspecial('gaussian',[9 9],0.25);
lambda_weight =10;
beta_weight=0;
tol = 1e-4; %We don't need to find the best in each step since we will tweak 2 variables t and g at the same time
niter=100;
bw =23;%This is the bandwidth parameter of the correlation ufnciton
h=fspecial('gaussian',[round(6*bw)+1 round(6*bw)+1],bw); %Transfer function of the low-pass filter...
overwriteslim_data=0;
%Create the directory
nz = 467;
for fovidx=1:nfovs
    dirname = strcat('fov_',num2str(fovidx-1));
    foldername =strcat(fpath,dirname); 
    if (~exist(foldername));
        mkdir(strcat(fpath,dirname)); %Create the directory if the folder doesn't exist.
    else
        disp(['Working at FOV' num2str(fovidx)]);
        for timestepidx=1:ntimesteps
             for zidx = 1:nz
              slimimname = strcat('org_',num2str(fovidx-1),'_',num2str(timestepidx-1),'_1_',num2str(zidx-1),'.tif');
             recimname= strcat('rec_',num2str(fovidx-1),'_',num2str(timestepidx-1),'_1_',num2str(zidx-1),'.tif');
             nnimname = strcat('nn_',num2str(fovidx-1),'_',num2str(timestepidx-1),'_1_',num2str(zidx-1),'.tif');
             disp(['Working at timestep ' num2str(timestepidx),' z ', num2str(zidx)]);
          
             if (~exist(strcat(foldername,'\',slimimname)))
                Aname = strcat(num2str(fovidx-1),'_',num2str(timestepidx-1),'_1_',num2str(zidx-1),'_0.tif');
                Bname = strcat(num2str(fovidx-1),'_',num2str(timestepidx-1),'_1_',num2str(zidx-1),'_1.tif');
                Cname = strcat(num2str(fovidx-1),'_',num2str(timestepidx-1),'_1_',num2str(zidx-1),'_2.tif');
                Dname = strcat(num2str(fovidx-1),'_',num2str(timestepidx-1),'_1_',num2str(zidx-1),'_3.tif');
                A=im2double(imread(strcat(fpath,Aname)));
                B=im2double(imread(strcat(fpath,Bname)));
                C=im2double(imread(strcat(fpath,Cname)));
                D=im2double(imread(strcat(fpath,Dname)));
                nrows = size(A,1);
                ncols = size(A,2);
                h1 = zeros(nrows,ncols);
                h1(1:size(h,1),1:size(h,2))=h;
                kernel_size=size(h,1);
                h1 = circshift(h1,[-round((kernel_size-1)/2) -round((kernel_size-1)/2)]);        

                [uphi,beta,del_phi,u0sqr,gamma]=combine_4_frame(A,B,C,D,zeros(nrows,ncols),...
                    2.45,1,nrows,1,ncols);
                uphi = unwrap2(uphi);
                uphi = cast(uphi,'single');
                figure(1);
                imagesc(uphi);colormap gray;
                colorbar;
                writeTIFF(uphi,strcat(foldername,'\',slimimname));
                disp('Done...');
             else
                 if (~exist(strcat(foldername,'\',nnimname)))
                         uphi = imread(strcat(foldername,'\',slimimname)); %Read the current file into the memory
                         uphi = cast(uphi,'double');
                         figure(1);
                         imagesc(uphi);
                         nnuphi = (uphi>=0).*uphi;
                         figure(2);
                         imagesc(nnuphi);
                         colorbar;
                         writeTIFF(nnuphi,strcat(foldername,'\',nnimname));
         
                         
                  
                 end
                 
                 %Run the reconstruction if needed
%                  if ((~exist(strcat(foldername,'\',recimname)))|(overwriteslim_data))
%                          uphi = imread(strcat(foldername,'\',slimimname)); %Read the current file into the memory
%                          uphi = cast(uphi,'double');
%                          figure(1);
%                         subplot(221);imagesc(uphi);colorbar;colormap jet
%                         title('SLIM phase image')
%                         uphi_denoised = imfilter(uphi,h_denoise,'same');
%                         idx = find(abs(uphi_denoised)>bg_thresh);
%                         mask = zeros(size(uphi));
%                         mask(idx)=1;
%                         mask = medfilt2(mask);
%                         mask = imclose(mask,strel('disk',5));
%                       %  mask = imfill(mask,'hole');
%                         subplot(222);
%                         imagesc(mask);
%                         libimphase = uphi_denoised;
%                         libimphase(idx)=NaN;
%                         phasebg = inpaint_nans(libimphase,2);
%                         subplot(223);
%                         imagesc(phasebg);colorbar;
%                         gamma = exp(i*(uphi_denoised-phasebg));
%                         if (gpuDeviceCount()==0)
%                                 gpu_compute_en = 0; %if there is no gpu, compute the result on cpu instead
%                         end
%                         figure(2);
%                         imagesc(mask);
%                         method = 'relax';
%                         if (gpu_compute_en==0)
%                             [gk,tk] = estimate_gt(gamma,h,niter,lambda_weight,beta_weight,tol,method);
%                         else %Compute gk and tk on gpu
%                             d = gpuDevice();
%                             reset(d); %Reset the device and clear its memmory
%                             [gk,tk] = estimate_gt_gpu(gamma,h,niter,lambda_weight,beta_weight,tol,x0,y0,half_cs_length,method,mask);
%                         end 
%                         figure(5);
%                         imagesc(angle(tk));colorbar; axis off;colormap gray;
%                         writeTIFF(unwrap2(cast(angle(tk),'double')),strcat(foldername,'\',recimname));
                  end
             end
          end
    end
end

  
 