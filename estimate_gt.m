function [gk,tk] = estimate_gt(gamma,h,niter,lambda,beta_weight,tol,method)
    %This function compute the estimation for gk and tk given gamma
    
   % maxgamma = max(gamma(:));%This value should be very close to 1
   % gamma = gamma/maxgamma;

    obj_array=zeros(0,1);
    nrows = size(gamma,1);
    ncols = size(gamma,2);
    gk = ones(size(gamma)).*exp(i*angle(gamma)); %Initial estimation
    tk = ones(size(gamma));
     
    %Fourier transform the filter
    h1 = zeros(nrows,ncols);
    h1(1:size(h,1),1:size(h,2))=h;
    kernel_size=size(h,1);
    h1 = circshift(h1,[-(kernel_size-1)/2 -(kernel_size-1)/2]);
    hf = fft2(h1);
   
    %Create the filter for the derivative filters
    gy1 = fspecial('sobel');
    gx1 = gy1';
    gx = zeros(nrows,ncols);
    gy = zeros(nrows,ncols);
    kernel_size = size(gy1,1);
    gx(1:size(gx1,1),1:1:size(gx1,2))=gx1;
    gy(1:size(gx1,1),1:1:size(gx1,2))=gy1;
    gx = circshift(gx,[-(kernel_size-1)/2 -(kernel_size-1)/2]);
    gy = circshift(gy,[-(kernel_size-1)/2 -(kernel_size-1)/2]);
    gxf = fft2(gx);
    gyf = fft2(gy);
    
    
    %Next, solve with the iterative method
    [obj,val1,val2,val3] = objective_comp(gamma,hf,gxf,gyf,tk,gk,lambda,beta_weight,nrows,ncols);
    obj_array(end+1)=obj;
    cjgamma = conj(gamma);
    disp(['Iter ' num2str(0) ': current objective: ' num2str(obj) ', 1st: ' num2str(val1),...
       ', 2nd: ' num2str(val2), ', 3rd: ' num2str(val3)]);
    val3_array = zeros(0,1);
    for iter=1:niter
        tic
        %First, recover g from t
        tkf = fft2(tk);
        gk = (tk.*cjgamma+lambda*ifft2(tkf.*hf))./(conj(tk).*tk+lambda+1e-8);
        
        switch method
            case 'relax'
                beta = norm(gk,'fro');
                betasqr = beta^2;
                rhs = betasqr*gamma./conj(gk)+lambda*Hhg_comp(hf,gk);
                rhsf = fft2(rhs);
                tkf = rhsf./(betasqr+lambda*abs(hf).^2+....
                    beta_weight*abs(gxf).^2+beta_weight*abs(gyf).^2+1e-8); %Added factor for stability
                tk = ifft2(tkf);
            case 'cg'
                rhs = gk.*gamma + lambda*Hhg_comp(hf,gk);        
                tk = cgs(@(x)A_comp(x,hf,gxf,gyf,lambda,beta_weight,gk,nrows,ncols),rhs(:),tol,20);
                tk = reshape(tk,[nrows ncols]);
        end
        [obj,val1,val2,val3] = objective_comp(gamma,hf,gxf,gyf,tk,gk,lambda,beta_weight,nrows,ncols);
        %obj_array(end+1)=obj;
        te = toc;
        disp(['Iter ' num2str(iter) ': current objective: ' num2str(obj) ', 1st: ' num2str(val1),...
        ', 2nd: ' num2str(val2), ', 3rd: ' num2str(val3)]);
        %Make sure that our phase is not offsett
        gamma_phase = angle(tk);
        gamma_phase = gamma_phase - mean(mean(gamma_phase(1:100,1:100)));
        tk = abs(tk).*exp(i*gamma_phase);
        figure(4);
        imagesc(angle(tk));colorbar;title('Phase (flipped) of tk');
        colormap gray
        figure(3);
        plot(angle(tk(1200,:)));drawnow

    end
    



end

function y=A_comp(x,hf,gxf,gyf,lambda,beta_weight,gk,nrows,ncols)
    %This function computes the results of (diag(gk.^2)+lambda*H^H*H)*x
    x = reshape(x,[nrows ncols]);
    xf = fft2(x);
    HhHf = conj(hf).*hf;
    GxhGxf = conj(gxf).*gxf;
    GyhGyf = conj(gyf).*gyf;
    
    yf = lambda*HhHf.*xf+beta_weight*(GxhGxf+GyhGyf).*xf;
    y = ifft2(yf);
    y = y + x.*conj(gk).*gk; %This one is faster than abs(gk).^2  
    y = y(:);
end

function Hhg=Hhg_comp(hf,gk)
    %This function compute the product H^H*gk
    gkf = fft2(gk);
    Hhgf = conj(hf).*gkf;
    Hhg = ifft2(Hhgf);      
end