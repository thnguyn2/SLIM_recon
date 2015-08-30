function demo
    %This function estimate the filter kernel from the edgemeasurement
    clc;
    clear all;
    close all;
    imname='org_phase_t.tif';
    %First, load the image of the pillar
    im = im2double(imread(imname));
   
    %Extract the image of the edge
    im = im(329:410,372:654);
    im = im(:,end:-1:1);
    im = im(:,1:end-1);
    
    %%Take the average line cross-section
    profile = sum(im,1)/size(im,1); %Sum over rows
    
    npoints = size(im,2);
    scale = 15.875;%pixels/microns
    ncols = size(im,2);
    distance_in_mu = ncols/scale;
    
    xx = linspace(-distance_in_mu/2,distance_in_mu/2,ncols);%Measured distance in microns
    deltax = xx(2)-xx(1);
    %Next, filter the proflle to remove the additive noise
    profilerec = conv(profile,gauss1d(1,30),'same');%Avoid measurement noise
    
    %Baseline correction to correct for sample tilt
    leftmean = mean(profile(1:10)); 
    rightmean=mean(profile(end-10:end));
    baseline = (xx-xx(1))/(xx(end)-xx(1))*(rightmean-leftmean)+leftmean;
    profilerec = profilerec -baseline; %Now, perform baseline correction
    profile = profile - baseline; %Line profile after correction. The left and right end should have 0 values now.
    phi_m=profilerec; 
  
    %%Prepare the step function 
    %Find the zero-crossing around the central
    phi_abs = abs(phi_m);
    start_loc = round(npoints*0.45);
    end_loc = round(npoints*0.55);
    [min_val,idx]=min(phi_abs(start_loc:end_loc));
    zcr_coord = start_loc + idx; %Location of the zero-crossing
    disp(['Location of the zero-crossing' num2str(zcr_coord)]);
    xx = xx-xx(zcr_coord-1); %The location of the jump should be 0
    
    %Make sure the edge filter respons have the same max and min vlaue
    [phi_max,max_loc] = max(phi_m);
    [phi_min,min_loc] = min(phi_m);
    %Scale the right and the left part of the plot, make sure they are
    %symmetric across the x-axis
    lhs = phi_m(1:zcr_coord);
    rhs = phi_m(zcr_coord+1:end);
    %Normalize the max and the min for phase profile to -1 and 1
    lhs = lhs/abs(phi_min);
    rhs = rhs/abs(phi_max);
    phi_m = [lhs,rhs];
    %%-----Done with pre-processing here-----%%
    
    %Display the phase cross-section and raw cross-sections
    figure(1);
    subplot(121)
    plot(xx,phi_m,'g','linewidth',2);
    h=xlabel('Distance ({\mu}m)');
    set(h,'FontSize',14);
    set(gca,'FontSize',14);
    grid on;
    h=ylabel('Phi [rad]');
    set(h,'FontSize',14);
    title('Final phase cross-section')
    subplot(122)
    plot(xx,im(end/2+5,:)-baseline,'b','linewidth',2);
    hold on;
    plot(xx,profile,'r','linewidth',2);
    axis([min(xx) max(xx) -0.3 0.3]);
    h_legend = legend('Raw profile','Mean profile');    set(h_legend,'FontSize',14);
    h=ylabel('arg(\Gamma_{r,s}(x)) [rad]');    set(h,'FontSize',14);
    h=xlabel('Distance ({\mu}m)');     set(h,'FontSize',14); 
    set(gca,'FontSize',14);
    grid on;

    
    %First result, process with Tan's edge method
    [h_mask,h_cs2,xx1d] = edge_process_tan_method(phi_m,xx);
    [xx,yy]=meshgrid(linspace(-distance_in_mu/2,distance_in_mu/2,size(h_mask,1)),...
    linspace(-distance_in_mu/2,distance_in_mu/2,size(h_mask,1)));
    figure(2);
    subplot(121);
    mesh(xx,yy,h_mask);colormap jet;title('PSF (Orignal method)')
    axis([-distance_in_mu/2,distance_in_mu/2,-distance_in_mu/2,distance_in_mu/2,min(h_cs2),max(h_cs2)]);
    axis off;
    grid on;
    save('h_slim.mat','h_mask','h_cs2','im');
    
    %Least square estimation and fit a Guassian kernel p(xi)=b*exp(-a*xi^2)
    %to the distribution
    xi = xx1d;
    xi = xi(:);
    yi = real(h_cs2(:));
    [mu,sigma,a,pfit]=gaussfit(xi,yi,0.1);
    %Draw the profile of the filter
    figure(3); plot(xx1d,abs(h_cs2),'r','linewidth',2);
    h=ylabel('|h_i|');    set(h,'FontSize',14);
    h=xlabel('Distance ({\mu}m)');    set(h,'FontSize',14);    set(gca,'FontSize',14);
    grid on; hold on;
    plot(xx1d,pfit/max(pfit),'-.b','linewidth',2);
    legend('Measured','Fitted')
    disp(['std in pixels: ' num2str(sigma*scale)])
    
    figure(4);
    plot(xx1d,real(h_cs2),'g','linewidth',2);
    h=ylabel('Value');    set(h,'FontSize',14);
    grid on;
    h=xlabel('Distance ({\mu}m)');    set(h,'FontSize',14);    set(gca,'FontSize',14);
    hold on;
    plot(xx1d,imag(h_cs2),'b','linewidth',2);
    plot(xx1d,abs(h_cs2),'k','linewidth',2);
    plot(xx1d,pfit/max(pfit),'-.r','linewidth',2);
    legend('Re(h_{i1})','Im(h_{i1})','|h_{i1}|','Fitted')
    set(gca,'FontSize',14);
    axis([-9 9 -0.4 1.2])
   
    
end

function [h_mask,h_cs2,xx1d] = edge_process_tan_method(phi_m,xx)
%This is the solver using Tan's method based on separability assumption
    %Crop out the section withing min_loc and max_loc
    [phi_max,max_loc] = max(phi_m);
    [phi_min,min_loc] = min(phi_m);
    phi_m = [phi_m(1:min_loc) phi_m(max_loc:end)];
    step_signal = zeros(size(phi_m));
    step_signal(min_loc+1:end)=2;
    phi_filtered = -phi_m+step_signal;
    %Adjust the coordinate of the zero-crossing
    xx=[xx(1:min_loc) xx(max_loc:end)-xx(max_loc)+xx(min_loc+1)];
    xx = xx - xx(min_loc+1);%Set the begining of the jump to be 0 coordinate
    
    %Smooth_out the signal to prevent the noise
    phi_filtered = conv1d(phi_filtered,gauss1d(10,80));%phi_filter is the original phase profile + step jump. Therefore, it is a smooth signal
    lhs = exp(i*phi_filtered);
    dlhs = lhs(2:end)-lhs(1:end-1);
    hi1 = dlhs;
    hi1f = fft(hi1);hi1f=real(hi1f);
    h_cs2 = ifft(hi1f);
    h_cs2=h_cs2(2:end);
    h_cs2 = h_cs2/max(h_cs2);

    
%   %This is the old code 
%    phi_m = phi_filtered;
%    phi_o = 2;    
%     %Compute the integral array
%     int_arr = (sin(phi_o-phi_m)./(sin(phi_o-phi_m)+sin(phi_m)));int_arr = flipud(int_arr(:));
%     %int_arr = int_arr(22:end-22); %Avoid the two ends
%     h_cs = int_arr(2:end)-int_arr(1:end-1);%This is after taking the derivative
%    %Make sure our filter cross-section is symmetric by zero-out the imaginary
%    %part of the fourier spectrum;
%    hf = fft(h_cs);hf = real(hf);
%    h_cs2 = ifft(hf);

   
    
    xx = xx(1:length(h_cs2));
    xx1d=xx;
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
                h_mask(y_coord,x_coord)=real(h_cs2(max_idx+distancemap(x_coord,y_coord)));
            end
        end
    end

end

function [mu,sigma,a,pfit]=gaussfit(x,y,thresh)
    %thresh: a thresh hold value to compare. If y is less than
    %ymax*thresh=> just zero-it out.
    %The fitting distribution is P = a*exp(-(x-mu)^2/2/sigma^2)
    %% cutting here
    index = (y > max(y)*thresh);
    p=polyfit(x(index), log(y(index)), 2); 
    sigma=sqrt(-1/(2*p(1))); 
    mu=p(2)*sigma^2;
    a=exp(p(3)+mu^2/(2*sigma^2));
    pfit = a*exp(-(x-mu).^2/2/(sigma^2));
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
    x_pad = ones(m+n-1,1);
    x_pad((n+1)/2:(n+1)/2+m-1)=x;
    x_pad(1:(n+1)/2)=mean(x(1:10));
    x_pad((n+1)/2+m:m+n-1)=mean(x(end-10:end));
    
    y = conv(x_pad,h,'same');
    y = y((n+1)/2:(n+1)/2+m-1);
end