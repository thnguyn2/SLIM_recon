function demo
    %This function estimate the filter kernel from the edgemeasurement
    clc;
    clear all;
    close all;
    %edgemeas = xlsread('Edgemeas_NA_0.09.xls'); %First column is the pixel idx, second colum is the measured phase
    
    dirpath = '.\Data_for_the_Halo_paper\mikhail\pillar_20\'; %path the files
    imname=strcat(dirpath,'acq3.tif');
    %First, load the image
    im = im2double(imread(imname));
    %im = im2double(imread('edge_im.tif'));
    %Rotate the image
    %im = imrotate(im,-2,'crop');
    figure(1);
    imagesc(im);colorbar
    im = im(915:957,855:969);
   % im = im(:,end:-1:1); %Flip the image
    figure(2);
    imagesc(im);title('Edge profile');
    im = im(:,1:end-1);
    profile = sum(im,1)/size(im,1); %Sum over rows
    npoints = size(im,2);
    hb_std = 2;
    hb = gauss1d(hb_std,npoints);
    hb = fftshift(hb);
    
    profilef = fft(profile);
    hbf = fft(hb);
    profilefrec = profilef.*conj(hbf)./(abs(hbf).^2+1e-2);
    profilerec = ifft(profilefrec);
    
    scale = 15.875;%pixels/microns
    ncols = size(im,2);
    distance_in_mu = ncols/scale;
    xx = linspace(-distance_in_mu/2,distance_in_mu/2,ncols);
    
    figure(1);
    plot(xx,profilerec,'g','linewidth',2);
    hold on
    figure(1);
    plot(xx,im(round(end/2+5),:),'b','linewidth',2);
    hold on;
    plot(xx,profile,'r','linewidth',2);
    axis([min(xx) max(xx) -0.3 0.3]);
    h_legend = legend('Current CS','Average CS');
    set(h_legend,'FontSize',14);
    h=ylabel('arg(\Gamma_{r,s}(x)) [rad]');
    set(h,'FontSize',14);
    h=xlabel('Distance ({\mu}m)');
    set(h,'FontSize',14);
    set(gca,'FontSize',14);
    grid on;

    profile = profilerec;
    phi_m=profile;

    %Find the offset signal.
    offsetval = (max(phi_m)+min(phi_m))/2;
    phi_m = phi_m-offsetval; %'m'=measurable
  
    %Find the zero-crossing around the central
    phi_abs = abs(phi_m);
    start_loc = round(npoints*0.45);
    end_loc = round(npoints*0.55);
    [min_val,idx]=min(phi_abs(start_loc:end_loc));
    zcr_coord = start_loc + idx; %Location of the zero-crossing
    disp(['Location of the zero-crossing' num2str(zcr_coord)]);
    phi_max = max(phi_m);
    phi_min = min(phi_m);

    step_signal = zeros(size(phi_m));
    step_signal(zcr_coord:end)=phi_max-phi_min;
    phi_filtered = -phi_m+step_signal;
    
    
    %Eliminate the effects of finite objective's NA that the zero-crossing is
    fitobj = fit(xx(:),phi_filtered(:),  'poly5');
    phi_filtered_2 = feval(fitobj,xx(:));
    
      
    phi_m =  phi_filtered_2;
    
    phi_o = max(phi_m);
    %Compute the integral array
    int_arr = (sin(phi_o-phi_m)./(sin(phi_o-phi_m)+sin(phi_m)));int_arr = flipud(int_arr(:));
    %int_arr = int_arr(22:end-22); %Avoid the two ends
    h_cs = int_arr(2:end)-int_arr(1:end-1);
    
    
    %Make sure our filter cross-section is symmetric by zero-out the imaginary
    %part of the fourier spectrum;
    hf = fft(h_cs);hf = real(hf);
    h_cs2 = ifft(hf);h_cs2=h_cs2(2:end);
    h_cs2 = h_cs2/max(h_cs2);
    figure(3); plot(xx(1:length(h_cs2)),h_cs2,'r','linewidth',2);
    h=ylabel('Amplitute');
    set(h,'FontSize',14);
    h=xlabel('Distance ({\mu}m)');
    set(h,'FontSize',14);
    set(gca,'FontSize',14);
    grid on;

  
    %Finally, create a 2D matrix that has the specified cross-section
    f_size = 2*floor((length(h_cs2)+1)/2)-1;
    half_f_size = (f_size-1)/2;
    h_mask = zeros(f_size,f_size);
    [max_val,max_idx]=max(h_cs2);
    coord_arr = -half_f_size:half_f_size;
    [xx,yy]=meshgrid(coord_arr,coord_arr);
    distancemap = round(sqrt(xx.^2+yy.^2));
    for x_coord=1:f_size
        for y_coord = 1:f_size
            if (distancemap(x_coord,y_coord)<half_f_size)
                h_mask(y_coord,x_coord)=h_cs2(max_idx+distancemap(x_coord,y_coord));
            end
        end
    end
    h_mask = h_mask/sum(sum(h_mask));
    [xx,yy]=meshgrid(linspace(-distance_in_mu/2,distance_in_mu/2,size(h_mask,1)),...
        linspace(-distance_in_mu/2,distance_in_mu/2,size(h_mask,1)));
    figure(4);surf(xx,yy,h_mask);colormap jet;
    axis([-distance_in_mu/2,distance_in_mu/2,-distance_in_mu/2,distance_in_mu/2,min(h_cs2),max(h_cs2)]);
    
    axis off;
%     set(gca,'FontSize',14);
%     h=xlabel('x ({\mu}m)');
%     set(h,'FontSize',14);
%     h=ylabel('y ({\mu}m)');
%     set(h,'FontSize',14);
   
    grid on;
    save(strcat(dirpath,'h_edge_915_to_957_855_to_969.mat'),'h_mask','h_cs2','im');
end

function h=gauss1d(std_x,npoints)
    %Geneate 1d gaussian filter
    x=linspace(-npoints/2,npoints/2,npoints);
    h = exp(-x.^2/2/std_x^2);
    h = h./sum(h);
end
function y=conv1d(x,h)
    %Produce 1D convolution with replicate
    m = length(x);
    n = length(h(:));
    x_pad = ones(m+n-1,1)*x(end);
    x_pad((n+1)/2:(n+1)/2+m-1)=x;
    y = conv(x_pad,h,'same');
    y = y((n+1)/2:(n+1)/2+m-1);
end