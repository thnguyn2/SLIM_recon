fname='G:\Slide1\2300.tif_Files\2300_p0';
% fname='G:\Data\Jun10\LiveNeuron\zsacn3\0018.tif_Files_rec'
I_back=im2double(imread('G:\Slide1\back\back.tif'));
%I_back=0;
fname_save='G:\Slide1\rec2300\'; 
minPhase=-0.5; %min(min(uphi));
maxPhase=4; %max(max(phase));
xx1=250; xx2=255; %1-1040
yy1=100; yy2=105; %1-1388
imRepTile=20000; %reptile, usually 64
imNum=1000; %time lapse, 
imZnum=1; % z stacks
tInterval=1;%1/2.6;
zInterval=4.08;%0.2 for 40X; 4.08 for 10X
xsize=1040;
ysize=1388;
Ncol=9;
Nrow=21;
fid = fopen(strcat(fname_save,'summary.txt'),'wt');
fprintf(fid, 'Time Interval %8.4f min\n',tInterval)
fprintf(fid, 'Z interval %8.4f um\n', zInterval)
fprintf(fid, 'xx1 = %4d, xx2 =%4d, yy1= %4d, yy2= %4d\n', xx1, xx2, yy1, yy2)
fprintf(fid, 'L1 = %8.4f\n', L1)
fprintf(fid, 'maximum phase %8.4f\n', maxPhase)
fprintf(fid, 'minimum phase %8.4f\n', minPhase)
fclose(fid);
for num=1:imNum;
for imRep=8729:imRepTile
for znum=1:imZnum;
A=RawReadREPZTimelapse(fname,imRep,num,znum,'c04'); A=A-I_back; %A=A-min(min(A));%A->I0 0pi  04
B=RawReadREPZTimelapse(fname,imRep,num,znum,'c01'); B=B-I_back; %B=B-min(min(B));%B->I1  01
C=RawReadREPZTimelapse(fname,imRep,num,znum,'c02'); C=C-I_back; %C=C-min(min(C));%C->I2  02
D=RawReadREPZTimelapse(fname,imRep,num,znum,'c03'); D=D-I_back; %D=D-min(min(D));%D->I3  03
%A=RawReadChannel(fname,num,'c01'); A=A-I_back; %A=A-min(min(A));%A->I0 0pi
%B=RawReadChannel(fname,num,'c04'); B=B-I_back; %B=B-min(min(B));%B->I1
%C=RawReadChannel(fname,num,'c03'); C=C-I_back; %C=C-min(min(C));%C->I2
%D=RawReadChannel(fname,num,'c02'); D=D-I_back; %D=D-min(min(D));%D->I3

Gs=D-B; Gc=A-C; del_phi=atan2(Gs,Gc);
%beta
L=(A-C+D-B)./(sin(del_phi)+cos(del_phi))/4; %E0*E1
g1=(A+C)/2; %E0^2+E1^2=s
g2=L.*L; %E0^2*E1^2=p
x1=g1/2-sqrt(g1.*g1-4*g2)/2; x2=g1/2+sqrt(g1.*g1-4*g2)/2; %solutions
beta1=sqrt(x1./x2); beta2=1./beta1;
%beta1=x1;beta2=x2;
%get constant from  average over pixels
cL=L(xx1:xx2,yy1:yy2); cbeta1=beta1(xx1:xx2,yy1:yy2);
%L1=real(mean2(cbeta1))/mean2(cL) %117.8434
kkk(imRep)=real(mean2(cbeta1))/mean2(cL);
real(mean2(cbeta1))/mean2(cL)
LL=L1*L;
beta=LL*2.3; %real beta, !!!!!!!!!!!!!!!!!!!!!!2.3 for 40X and 4 for 10X

phi=atan2(beta.*sin(del_phi),1+beta.*cos(del_phi));
%unwrap2
uphi=unwrap2(phi);

uphi=medfilt2(uphi);
%uphi=imsmooth(uphi,3);
uphi=uphi-bg; %%%%%background substraction

backdiff(num)=sum(sum(uphi));
%uphi=uphi-backdiff;

imagesc(uphi,[-0.4,0.8]);colorbar  %-0.4, 1.5 for dry mass
title(strcat('p=',num2str(imRep),' t=',num2str(num*tInterval, '%04.0f\n'),'min'));
%title(strcat('p=',num2str(imRep)));
%title(strcat('z=',num2str(znum*0.2, '%04.1f\n'), '\mum')); %*900/imNum
%title(strcat('t=',num2str(num*tInterval, '%04.0f\n'),'sec', ', z=',num2str(znum*zInterval, '%04.1f\n'), '\mum')); %*900/imNum
%title(strcat('t=',num2str(num*tInterval, '%04.0f\n'),'sec' )); %*900/imNum 6.1f
axis image
axis off

mm=num2str(num)
zmm=num2str(znum)
pmm=num2str(imRep)
if znum<10
    zmm=strcat('0',zmm)
end
if num<10
    mm=strcat('0',mm)
end
if imRep<10
    pmm=strcat('0',pmm)
end
print('-djpeg', strcat(fname_save,'p',pmm,'t',mm,'z',zmm,'.jpg')); %-djpeg

phase=uphi-minPhase; 
uniPhase=phase/maxPhase;
mPhase=num2str(maxPhase);
imwrite(uint16((2^16-1)*uniPhase),strcat(fname_save,'tif\rec_',mPhase,'_','p',pmm,'t',mm,'z',zmm,'.tif'),'tiff','Compression','none');
%imwrite(uint16((2^16-1)*(del_phi+pi)/2/pi),strcat(fname_save,'delphi\rec_','p',pmm,'t',mm,'z',zmm,'.tif'),'tiff','Compression','none');
T_cos=cos(uphi);
T_cos=(-T_cos+1)/2;
imwrite(uint16((2^16-1)*(T_cos)),strcat(fname_save,'cos\p',pmm,'t',mm,'z',zmm,'.tif'),'tiff','Compression','none');

uphi_nng=uphi; 
for k=1:size(uphi_nng,1)
    for m=1:size(uphi_nng,2)
        if uphi_nng(k,m)<0
            uphi_nng(k,m)=0;
        end
    end
end
phase1=uphi_nng-minPhase; 
uniPhase1=phase1/maxPhase;
mPhase=num2str(maxPhase);
imwrite(uint16((2^16-1)*uniPhase1),strcat(fname_save,'tif_nng\rec_',mPhase,'_','p',pmm,'t',mm,'z',zmm,'.tif'),'tiff','Compression','none');
%imwrite(uint16((2^16-1)*(del_phi+pi)/2/pi),strcat(fname_save,'delphi\rec_','p',pmm,'t',mm,'z',zmm,'.tif'),'tiff','Compression','none');
T_cos1=cos(uphi_nng);
T_cos1=(-T_cos1+1)/2;
imwrite(uint16((2^16-1)*(T_cos1)),strcat(fname_save,'cos_nng\p',pmm,'t',mm,'z',zmm,'.tif'),'tiff','Compression','none');
%if floor((floor((imRep-1)/Ncol)+1)/2)~=(floor((imRep-1)/Ncol)+1)/2
%uniPhase_st(xsize*(floor((imRep-1)/Ncol))+1:xsize*(floor((imRep-1)/Ncol)+1),ysize*(mod(imRep-1,Nrow))+1:ysize*(mod(imRep-1,Nrow)+1))=uniPhase;
%uniPhase1_st(xsize*(floor((imRep-1)/Ncol))+1:xsize*(floor((imRep-1)/Ncol)+1),ysize*(mod(imRep-1,Nrow))+1:ysize*(mod(imRep-1,Nrow)+1))=uniPhase1;
%else
  %uniPhase_st(xsize*(floor((imRep-1)/Ncol))+1:xsize*(floor((imRep-1)/Ncol)+1),ysize*(Nrow-mod(imRep-1,Nrow)-1)+1:ysize*(Nrow-mod(imRep-1,Nrow)))=uniPhase;
  %uniPhase1_st(xsize*(floor((imRep-1)/Ncol))+1:xsize*(floor((imRep-1)/Ncol)+1),ysize*(Nrow-mod(imRep-1,Nrow)-1)+1:ysize*(Nrow-mod(imRep-1,Nrow)))=uniPhase1;
  
%end
end
end
%imwrite(uint16((2^16-1)*uniPhase_st),strcat(fname_save,'tif_st\rec_',mPhase,'_','p',pmm,'t',mm,'z',zmm,'.tif'),'tiff','Compression','none');
%imwrite(uint16((2^16-1)*uniPhase1_st),strcat(fname_save,'tif_nng_st\rec_',mPhase,'_','p',pmm,'t',mm,'z',zmm,'.tif'),'tiff','Compression','none');
end

