% function IABfrom4PiPSFfit(PSF,startpar)
% startpar=[dx dy dz phase norm, frequency] all vectors length 4 except frequency;
% dx=par(1:3);
% dy=par(4:6);
% dz=par(7:9);
% normf=par(10:12);
% phaseshift=par(13);
% frequency=par(14);
%shift PSF by dx,dy,dz

s=size(PSF);
rx=ceil((9+1)/2);rz=40;mp=ceil((s(1)+1)/2);mpz=ceil((s(3)+1)/2);
PSFs=PSF(mp-rx:mp+rx,mp-rx:mp+rx,mpz-rz:mpz+rz,:);
% PSFs=PSFal(mp-rx:mp+rx,mp-rx:mp+rx,mpz-rz:mpz+rz,:);

% PSFs(:,:,:,1)=PSFs(:,:,:,1)*.8;
frequency=0.2247;
phaseshiftp=1.42;
if phaseshiftp>1
    phaseshiftp=phaseshiftp-2;
end
phaseshift=phaseshiftp*pi;

startpar=[0 0 0 0 0 0 0 0 0 1 1 1];
% startpar=[0 0 0 0 0 0 0 0 0 1 1 1 phaseshift frequency];
fixpar=[phaseshift frequency];
PSFstart=recoverPSF(startpar,PSFs,fixpar);

options=optimset('fminsearch');
options.Display='iter';
fitmn=fminsearch(@errPSF,startpar,options,PSF,fixpar);

options=optimset('lsqcurvefit');
% options.Algorithm='levenberg-marquardt';
% options.Diagnostics='on';
options.FinDiffRelStep=0.01;
PSFc=PSFs(2:end-1,2:end-1,2:end-1,:);
fitpar=lsqcurvefit(@recoverPSF,startpar,PSFs,PSFc,[],[],options,fixpar);


PSFrecovered=recoverPSF(fitpar,PSFs,fixpar);
PSFrecoveredms=recoverPSF(fitmn,PSFs,fixpar);

errfmin=errPSF(fitmn,PSF,fixpar)
errlsq=errPSF(fitmn,PSF,fitpar)
errfmin-errlsq

function err=errPSF(par,PSF,fixpar)
PSFrec=recoverPSF(par,PSF,fixpar);
err=(((PSFrec-PSF(2:end-1,2:end-1,2:end-1,:)).^2));
err=sum(err(:));
end

function PSFm=recoverPSF(par,PSF,fixpar)

dx=par(1:3);
dy=par(4:6);
dz=par(7:9);
normf=par(10:12);
% normf=ones(1,3);
% phaseshift=par(13);
% frequency=par(14);
phaseshift=fixpar(1);
frequency=fixpar(2);

phaseshifts=[0 phaseshift pi phaseshift+pi]-pi;


%shift PSF and rescale
xn=1:size(PSF,1);yn=1:size(PSF,2);zn=1:size(PSF,3);
[Xq,Yq,Zq]=meshgrid(yn,xn,zn);
PSFS=zeros(size(PSF));
PSFS(:,:,:,1)=PSF(:,:,:,1);
for k=2:4
    PSFS(:,:,:,k)=interp3(PSF(:,:,:,k),Xq-dx(k-1),Yq-dy(k-1),Zq-dz(k-1),'cubic',0)/normf(k-1);
end

% calculate IAB
[I,A,B]=make4Pimodel(PSFS,phaseshifts,frequency);

%make PSF again
PSFm=makePSF(I,A,B,frequency, phaseshifts, [1 normf]);
PSFm=PSFm(2:end-1,2:end-1,2:end-1,:);
end


function [I,A,B]=make4Pimodel(allPSFs,phaseshifts,frequency)
%re-weight every PSF by relative transmission?
I1=(allPSFs(:,:,:,1)+allPSFs(:,:,:,3))/2;
I2=(allPSFs(:,:,:,2)+allPSFs(:,:,:,4))/2;
Iall=(I1+I2)/2;

z=(1:size(allPSFs,3))'-round(size(allPSFs,3)/2);
[A12,B12]=makeAB(allPSFs(:,:,:,1),allPSFs(:,:,:,2),Iall,z,frequency,phaseshifts(1),phaseshifts(2));
[A41,B41]=makeAB(allPSFs(:,:,:,4),allPSFs(:,:,:,1),Iall,z,frequency,phaseshifts(4),phaseshifts(1));
[A23,B23]=makeAB(allPSFs(:,:,:,2),allPSFs(:,:,:,3),Iall,z,frequency,phaseshifts(2),phaseshifts(3));
[A34,B34]=makeAB(allPSFs(:,:,:,3),allPSFs(:,:,:,4),Iall,z,frequency,phaseshifts(3),phaseshifts(4));
A=(A12+A23+A34+A41)/4;
B=(B12+B23+B34+B41)/4;
I=Iall;
end

function [A,B]=makeAB(P1,P2,I,z,frequency,phase1,phase2)
    A=zeros(size(I));B=zeros(size(I));
    for k=1:length(z)
        a1=2*frequency*z(k)+phase1;
        a2=2*frequency*z(k)+phase2;
        A(:,:,k)=(sin(a1).*(P2(:,:,k)-I(:,:,k))-sin(a2).*(P1(:,:,k)-I(:,:,k)))./(cos(a2).*sin(a1)-cos(a1).*sin(a2));
        B(:,:,k)=(-cos(a1).*(P2(:,:,k)-I(:,:,k))+cos(a2).*(P1(:,:,k)-I(:,:,k)))./(cos(a2).*sin(a1)-cos(a1).*sin(a2));
    end
end

function PSFo=makePSF(I,A,B,frequency, phaseshifts, normf)
s=size(I);
PSF=zeros(s(1)*s(2),s(3),4);
z=(1:s(3))-round(s(3)/2);
Ir=reshape(I,s(1)*s(2),s(3));
Br=reshape(B,s(1)*s(2),s(3));
Ar=reshape(A,s(1)*s(2),s(3));
for k=1:4
    PSF(:,:,k)=normf(k)*(Ir+Ar.*cos(2*frequency*z+phaseshifts(k))+Br.*sin(2*frequency*z+phaseshifts(k)));
end

PSFo=reshape(PSF,s(1),s(2),s(3),4);
end