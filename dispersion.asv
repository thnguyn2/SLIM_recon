dir_home='K:\20101117\SampleI\Before\temp2\';
dir_fi=dir(strcat(dir_home,'*_rec*'));
n_set=size(dir_fi,1); % number of sets
dir_start=str2num(dir_fi(1).name(1:3)); % first directory

for curr_dir=1:n_set
% for curr_dir=4:n_set
    
    
raw_dir=strcat(dir_home,num2str(dir_start+curr_dir-1),'_rec\tif\');
fname=strcat(raw_dir,'rec_4','_p0');
fname_save=strcat(dir_home,num2str(dir_start+curr_dir-1),'_rec\tif\');

im_files=dir(strcat(raw_dir,'*z01.tif'));
%[z,t]=get_info(im_files);

FileName =strcat(dir_home,num2str(dir_start+curr_dir-1),'_rec\tif\tif-scale_512cut.tif');
ave_phi = im2double(imread(strcat(dir_home,num2str(dir_start+curr_dir-1),'_rec\tif\AVG_tif-scale_512cut.tif')));
%fname_save=strcat('E:\Process\127_rec\tif\')
ave_phi = ave_phi*4-0.5 ;
fft_ave_phi = fftshift(fft2(ave_phi));

NumIm = 601;
%NumIm = size(im_files);
time_step = 1;
RowSize = 512;
ColumnSize = 512;
numerator = zeros(RowSize,ColumnSize) ;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
denominator = zeros(RowSize,ColumnSize) ;
for iframe = 1:NumIm -1 
    
  B_1 = im2double(imread(FileName,iframe)) ;  
  B_1 = B_1*4-0.5;
  B_2 = im2double(imread(FileName,iframe+1)) ;                                                                                                                                                 
  B_2 = B_2*4-0.5;
  real_phi_1 = real(fftshift(fft2(B_1))) ;
  img_phi_1 = imag(fftshift(fft2(B_1))) ;
  
  real_phi_2 = real(fftshift(fft2(B_2))) ;
  img_phi_2 = imag(fftshift(fft2(B_2))) ;  
  modu_phi_dt = ((real_phi_2-real_phi_1)/time_step).^2 + ((img_phi_2-img_phi_1)/time_step).^2 ;
  numerator = numerator + modu_phi_dt*time_step ;
  modu_phi = abs(fftshift(fft2(B_1))-fft_ave_phi).^2 ;
  denominator = denominator + modu_phi*time_step ;  
end

omega_bar_squa_qx_qy = numerator ./denominator  ;
omega_bar_qx_qy = sqrt((omega_bar_squa_qx_qy)) ;
imagesc(omega_bar_qx_qy) ;

% min_omega_bar = min(min(omega_bar_qx_qy)) 
% max_omega_bar = max(max(omega_bar_qx_qy)) 
% omega_bar = (omega_bar_qx_qy - min_omega_bar )/(max_omega_bar-min_omega_bar);
% imwrite(uint16((2^16-1)*omega_bar),strcat('result','.tif'),'tiff','Compression','none');

% 
for i = 1:RowSize
    for j = 1:ColumnSize
       i,j
        omega_bar_squ = omega_bar_qx_qy(i,j)^2 ;
        omega_max = 2*pi/time_step ;
        f = @(delta_omega)omega_bar_squ *atan(omega_max/delta_omega)- delta_omega*omega_max+delta_omega^2*atan(omega_max/delta_omega);
        delta_omega = fzero(f,1) ;
        omega_bar_correct(i,j)= abs(delta_omega);
    end
end

% Intercept = -0.03355 ;
% B1 = 0.10398;
% B2 = 0.10368;
% B3 = 0.077;
% 
% for i = 1 : RowSize
%     for j = 1: ColumnSize
%         i, j
%         omega_bar_correct(i, j) = Intercept + B1 * omega_bar_qx_qy(i,j)+B2 * omega_bar_qx_qy(i,j)^2+B3* omega_bar_qx_qy(i,j)^3;
%     end
% end

min_omega_bar_correct = min(min(omega_bar_correct)) 
max_omega_bar_correct = max(max(omega_bar_correct)) 
omega_bar_img = (omega_bar_correct - min_omega_bar_correct)/(max_omega_bar_correct-min_omega_bar_correct);
imwrite(uint16((2^16-1)*omega_bar_img),strcat(fname_save,'result','.tif'),'tiff','Compression','none');
fid = fopen(strcat(fname_save,'result.txt'),'wt')
fprintf(fid, '0   %8.8f \n',min_omega_bar_correct)
fprintf(fid, '65535   %8.8f \n', max_omega_bar_correct)
fclose(fid)
%clear NumIm;
%clear im_files;
end