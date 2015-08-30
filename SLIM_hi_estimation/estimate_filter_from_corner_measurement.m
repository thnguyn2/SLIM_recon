function demo
    %This function estimate the filter kernel from the edgemeasurement
    clc;
    clear all;
    close all;
    imname='org_phase_t.tif';
    %First, load the image of the pillar
    im = im2double(imread(imname));
   
    %Extract the image of the corner
    im = im(130:335,430:600);
    
    figure(3);
    subplot(121);imagesc(im);title('Original corner image')
    npoints = size(im,2);
    scale = 15.875;%pixels/microns
    ncols = size(im,2);%Get the number of rows and number of columns in the current image
    nrows = size(im,1);
    distance_in_x_mu = ncols/scale;
    distance_in_y_mu = nrows/scale;
   
    
    xx = linspace(-distance_in_x_mu/2,distance_in_x_mu/2,ncols);%Measured distance in microns
    yy = linspace(-distance_in_y_mu/2,distance_in_y_mu/2,nrows);%Measured distance in microns
    
    h = fspecial('gaussian',12,2);
    %Next, filter the corner image to suppress the noise
    im = imfilter(im,h,'same');
    
    %Remove smoothed area due to the effects of the objective. This step
    %will generate some error. Works fine if the NAc << NAo
    r1 = 83;
    r2 = 108;
    c1 = 79;
    c2 = 92;
    %im = im([1:r1,r2:nrows],[1:c1,c2:ncols]);
    
    %Make sure that the min and the max value in the image are the same
    phi_max = max(im(:));
    phi_min = min(im(:));
    %Scale the phase to range of [-1,1]
    im = (im-phi_min)*2/(phi_max-phi_min)-1;
    subplot(122);imagesc(im);title('Noise suppressed image');colorbar
   
    gamma = exp(i*im);
    gammaconj = conj(gamma);
    
    %Next, generate the step in the corner profile
    phasestep = zeros(size(im));
    phasestep(r1+1:end,1:c1)=2;
    expiu = exp(i*phasestep);
    lhs = gammaconj.*expiu;
    
    

    %Another filtering step to make sure the transient is smooth
    h = fspecial('gaussian',12,0.5);
    %Next, filter the corner image to suppress the noise
    lhs = imfilter(lhs,h,'same');
   
    
    [gx,gy]=gradient(lhs);
    [gxx,gxy]=gradient(gx);
    
    figure(2);
    subplot(121);imagesc(angle(lhs));title('Angle of the lhs');colorbar
    subplot(122);imagesc(abs(gxy/(exp(i*2)-1)));
    
    
   
end

