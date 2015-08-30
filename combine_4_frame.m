function [uphi,beta,del_phi,u0sqr,gamma]=combine_4_frame(A,B,C,D,I_back,alpha,xx1,xx2,yy1,yy2)
    %Compute the SLIM image from 4 frames A, B, C, D are 0, pi/2, pi,
    %3*pi/2
    if (nargin>=5)
        A=A-I_back; %A->I0 0pi
        B=B-I_back; %B->I1
        C=C-I_back; %C->I2
        D=D-I_back; %D->I3
    end
    if (nargin==4)
        alpha = 2.45;
    end
    Gs=D-B; Gc=A-C; del_phi=atan2(Gs,Gc);

    %beta
    L=(A-C+D-B)./(sin(del_phi)+cos(del_phi))/4; %E0*E1
    g1=(A+C)/2; %E0^2+E1^2=s
    g2=L.*L; %E0^2*E1^2=p
    x1=g1/2-sqrt(g1.*g1-4*g2)/2; x2=g1/2+sqrt(g1.*g1-4*g2)/2; %solutions
    beta1=sqrt(x1./x2); beta2=1./beta1;
    if (nargin>=8)
        cL=L(xx1:xx2,yy1:yy2); cbeta1=beta1(xx1:xx2,yy1:yy2);
        L1=real(mean2(cbeta1))/mean2(cL)    %[Tan]: this is 1/<Uo>^2
        LL=L1*L;
        beta=LL*alpha; %real beta, 2.3 for 40X and 4 for 10X
    end
    phi=atan2(beta.*sin(del_phi),1+beta.*cos(del_phi));
    %unwrap2
    uphi=unwrap2(phi);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    uphi=medfilt2(uphi);
    u0sqr = A./(1/alpha+beta).^2;
    gamma = u0sqr.*(1+beta.*exp(i*del_phi));
end