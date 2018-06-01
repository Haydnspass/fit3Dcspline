function calibrate_4pi(p)
%for now now segmentation of calibration. 
ph=p;
% f=figure('Name','Bead calibration');
f=figure(44);
tg=uitabgroup(f);
t1=uitab(tg,'Title','prefit');
tgprefit=uitabgroup(t1);
ph.tabgroup=   tgprefit;
ph.isglobalfit=false;
ph.outputfile={};

settings_3D=readstruct('settings_3D.txt'); %later to settings, specify path in gui    
ph.settings_3D=settings_3D;
[beads,ph]=images2beads_globalfit(ph); %later extend for transformN

sstack=size(beads(1).stack.image);

% ybeads=getFieldAsVectorInd(beads,'pos',2);
xbeads=getFieldAsVectorInd(beads,'pos',1);

% y4pi=settings_3D.y4pi; %later: load from here! add field to GUI for selecting file.
% x4pi=settings_3D.x4pi;
% height4pi=settings_3D.height4pi; 
% width4pi=settings_3D.width4pi; 
% mirror4pi=settings_3D.mirror4pi;
fw=round((p.zcorrframes-1)/2);
framerange=round(sstack(3)/2-fw:sstack(3)/2+fw);
numchannels=length(settings_3D.y4pi);
allPSFs=zeros(sstack(1),sstack(2),sstack(3),numchannels);
for k=1:numchannels
%     h4pi=ph.settings_3D.height4pi;
%     indbh=ybeads>=(k-1)*h4pi+1 & ybeads<= k*h4pi;
    w4pi=ph.settings_3D.width4pi;
    indbh=xbeads>=(k-1)*w4pi+1 & xbeads<= k*w4pi;
    
    xposh=getFieldAsVectorInd(beads(indbh),'pos',1)';
    yposh=getFieldAsVectorInd(beads(indbh),'pos',2)';
%     ph.tabgroup=uitabgroup(uitab(tgprefit,'Title',num2str(k)));
    [allstacks,filenumber]=bead2stack(beads(indbh));
    [allPSFs(:,:,:,k),shiftedstack,shift,indgood]=registerPSF3D_g(allstacks,[],struct('framerange',framerange));
   
    
    beadtrue{k}(:,1)=xposh(indgood)-shift(indgood,2)'; %experimentally: this works :)
    beadtrue{k}(:,2)=yposh(indgood)-shift(indgood,1)';
    beadtrue{k}(:,3)=shift(indgood,3)'; %this is not yet tested, could be minus
    
end




%furhter align PSFs. Maybe this is not needed. But this gives x,y shift
%between PSFs in different channels. Not necessary to take into account?
framerange=round(max(1,sstack(3)/2-2*fw):min(sstack(3)/2+2*fw,sstack(3)));
[~,PSFaligned,shift,indgood]=registerPSF3D_g(allPSFs,[],struct('framerange',framerange,'removeoutliers',false),{},filenumber);

%try to correct for that shift:
for k=1:numchannels
    
    beadtrue{k}(:,1)= beadtrue{k}(:,1)-shift(k,2);
    beadtrue{k}(:,2)=beadtrue{k}(:,2)-shift(k,1);
    beadtrue{k}(:,3)=beadtrue{k}(:,3)-shift(k,3);
end

%calculate transformN
transform=interfaces.LocTransformN;
pt.mirror=0; %mirror already taken care of when reading in images
% pt.mirror=settings_3D.mirror4pi(1)*2; %2 for y coordinate
% pt.xrange=[settings_3D.x4pi(1) settings_3D.x4pi(1)+settings_3D.width4pi];
% pt.yrange=[settings_3D.y4pi(1) settings_3D.y4pi(1)+settings_3D.height4pi];
settings_3D=ph.settings_3D;
pt.xrange=[1 settings_3D.width4pi];
pt.yrange=[1 settings_3D.height4pi];
pt.unit='pixel';
pt.type='projective';
transform.setTransform(1,pt)
iAaa=1:size(beadtrue{1},1);
for k=2:length(beadtrue)
%     pt.mirror=settings_3D.mirror4pi(k)*2;
%     pt.yrange=[(k-1)*settings_3D.height4pi+1 k*settings_3D.height4pi];
    pt.xrange=[(k-1)*settings_3D.width4pi+1 k*settings_3D.width4pi];
    transform.setTransform(k,pt)
    tab=(uitab(tgprefit,'Title',['T' num2str(k)]));ph.ax=axes(tab);
    [~ ,iAa,iBa]=transform_locs_simpleN(transform,1, beadtrue{1},k,beadtrue{k},ph);
    %if use PSFaligned: add to beadtrue{k} etc the shift. sign?
    iAaa=intersect(iAa,iAaa);
end

out.transformation=transform;

tab=(uitab(tgprefit,'Title','frequency'));ph.ax=axes(tab);

%extend: also fix phi3-phi1=phi4-phi3=pi?
[phaseh,frequency]=getphaseshifts(allPSFs,ph.ax);
phaseshifts=[phaseh(1) phaseh(2) phaseh(1)+pi phaseh(2)+pi]; 
phaseshifts=phaseshifts-phaseshifts(1)-pi;
[I,A,B]=make4Pimodel(allPSFs,phaseshifts,frequency);

plotI(:,:,:,1)=I;plotI(:,:,:,2)=A;plotI(:,:,:,3)=B;plotI(:,:,:,4)=allPSFs(:,:,:,1);
tab=(uitab(tgprefit,'Title','IAB'));imageslicer(plotI,'Parent',tab)

%fit with this intermediate PSF model -> x,y,z, phi for all beads
fit4pidir=strrep(pwd,'SMAP',['ries-private' filesep 'PSF4Pi']);
if exist(fit4pidir,'file')
    addpath(fit4pidir);
end
mp=ceil((size(A,1)-1)/2);
dd=floor((p.ROIxy-1)/2);
PSF.Aspline=single(getsmoothspline(A(mp-dd:mp+dd,mp-dd:mp+dd,:),p));
PSF.Bspline=single(getsmoothspline(B(mp-dd:mp+dd,mp-dd:mp+dd,:),p));
PSF.Ispline=single(getsmoothspline(I(mp-dd:mp+dd,mp-dd:mp+dd,:),p));
PSF.frequency=frequency;
PSF.phaseshifts=phaseshifts;

%do fitting
fitroi=13;
sim=size(allPSFs);
mp=floor((sim(1)-1)/2)+1;
droi=floor((fitroi-1)/2);
rangeh=mp-droi:mp+droi;
phi0=phaseshifts;
%fit calibrations stack
shared=[0,0,1,1,1,1];
% [Pc,CRLB1, LL,update, error,residualc] =  kernel_MLEfit_Spline_LM_multichannel_4pi(allPSFs(rangeh, rangeh, :, :)*10000,PSF, shared,zeros(6,4,size(allPSFs,3)),50,phi0);
% [Pc,CRLB1, LL,update, error,residualc] =  CPUmleFit_LM_MultiChannel(single(allPSFs(rangeh, rangeh, :, :)*10000),uint32(shared),zeros(6,4,size(allPSFs,3)),50,phi0);

z0=round(size(PSF.Aspline,3)/2);
dTAll=zeros(6,4,size(allPSFs,3),'single');
iterations=50;
imstack=allPSFs(rangeh, rangeh, :, :)*10000;
[Pc,CRLB1 LL] = CPUmleFit_LM_MultiChannel(single(imstack(:, :, :, :)),uint32(shared),int32(iterations),single(PSF.Ispline), single(PSF.Aspline),single(PSF.Bspline),single(dTAll),single(phi0),single(z0));


mean(Pc(:,1:8),1)-droi+1

%this is not all the same -> PSFs in channels not well aligned. 

%previous transform: corresponding beads: %this should be part of a new
%images2beads function
    %round pos1: beadtrue{1}
    ph.isglobalfit=true;
    ph.transformation=transform;
    [beads,ph]=images2beads_globalfitN(ph); 
    

[imstack,fn,dxy]=bead2stack(beads);
sim=size(imstack);
imsqueeze=reshape(imstack,sim(1),sim(2),[],sim(end));
dTAll=reshape(dxy,size(dxy,1),sim(end),[]);
% phi0 = [0, pi, pi / 2, 1.5 * pi];

shared=[1,1,1,1,1,1];
% shared=[0,0,0,0,0,0,0];
%     dTAll=permute(dTAll,[2 3 1]);
    
   % dT: coord, channel,loc
% [P_EL, update_EL, error_EL, model_EL] =  kernel_MLEfit_Spline_LM_newFit(imstack(:, :, :, :), PSF, RoiPixelsize, 50, phi0, z0);




global P
% [P,CRLB1, LL,update, error,residual] =  kernel_MLEfit_Spline_LM_multichannel_4pi(imsqueeze(rangeh, rangeh, :, :),PSF, shared,dTAll,50,phi0);
imstack=imsqueeze(rangeh, rangeh, :, :);
[P,CRLB1 LL] = CPUmleFit_LM_MultiChannel(single(imstack(:, :, :, :)),uint32(shared),iterations,single(PSF.Ispline), single(PSF.Aspline),single(PSF.Bspline),single(dTAll),single(phi0),z0);
% imageslicer(residual)

shared(1:2)=0;
[Pu,CRLB1 LL] = CPUmleFit_LM_MultiChannel(single(imstack(:, :, :, :)),uint32(shared),iterations,single(PSF.Ispline), single(PSF.Aspline),single(PSF.Bspline),single(dTAll),single(phi0),z0);
dx21=Pu(:,2)-Pu(:,1);
% imageslicer(residual)
 
phase=mod(reshape(P(:,6),[],sim(4)),2*pi);
zphase=phase/2/frequency*p.dz;
zastig=reshape(P(:,5),[],sim(4))*p.dz;

xfit=reshape(P(:,1),[],sim(4));
yfit=reshape(P(:,2),[],sim(4));

z_phi = reshape(z_from_phi_JR(P(:, 5), phase(:), frequency, ceil(sim(3)/2)-.7),[],sim(4))*p.dz;
figure(68)
subplot(2,3,1)
plot(zastig)
xlabel('frame')
ylabel('z_astig')
subplot(2,3,2)
plot(phase)
xlabel('frame')
ylabel('phase')
subplot(2,3,3)
plot(zastig,zphase)
xlabel('z_astig')
ylabel('z_phase')
subplot(2,3,4)
plot(z_phi)
xlabel('frame')
ylabel('z_phi')
subplot(2,3,5)
plot(xfit,yfit,'+')
xlabel('x')
ylabel('y')
subplot(2,3,6)
hold off
plot(zastig,xfit)
hold on
xlabel('z_astig')
ylabel('x')
% plot(zastig,yfit)
    %cut out rois
    %fit unlinked
    %determine true x,y,z,phi for each bead


% find bead pairs
% cut out beads and move according to fit to perfect overlap
% calculate I, A, B for every bead
% perform registerPSF on I,A,B and all 4 channels together
% tilted coverslip: do we need to adjust phase for every bead to have same
% A, B, I? Use fitted phase phi for htis?

%alternatively: 
%as before for 2 channels average a 4 channel PSF (with fringes). Then use
%this for A,B,I


out.cal4pi.coeff=PSF;
out.cal4pi.dz=ph.dz;
out.cal4pi.x0=ceil((ph.ROIxy+1)/2);
out.cal4pi.z0=ceil((size(PSF.Aspline,3)+2)/2);
out.cal4pi.transformation=transform;
out.Xrange=ph.xrange;
out.Yrange=ph.yrange;
out.EMon=ph.emgain;


parameters=rmfield(ph,{'tabgroup','ax','status','fileax','smappos'});
out.parameters=parameters;
  
p.status.String='save calibration';drawnow
if ~isempty(p.outputfile)
        save(p.outputfile,'-struct','out');
end
end


function [stack,filenumber,dT]=bead2stack(beads)
ss=size(beads(1).stack.image);
if length(ss)==3
    stack=zeros(ss(1),ss(2),ss(3),length(beads));
    for k=length(beads):-1:1
        stack(:,:,:,k)=beads(k).stack.image;
        filenumber(k)=beads.filenumber;
    end
elseif length(ss)==4
    stack=zeros(ss(1),ss(2),ss(3),length(beads),ss(4));
    numpar=6;
    dT=zeros(numpar,ss(4),ss(3),length(beads));
    for k=length(beads):-1:1
        stack(:,:,:,k,:)=beads(k).stack.image;
        filenumber(k)=beads.filenumber;
        for zz=1:ss(3)
            dT(1:2,:,zz,k)=squeeze(beads(k).shiftxy(1,[2 1],:));
        end
    end    
end
end

function [phaseshiftso,frequencyo]=getphaseshifts(allPSFs,ax)
ss=size(allPSFs);
range=(ss(1)+1)/2+[-1 1];
fw=20;
frange=round(ss(3)/2+(-fw:fw)');

f=(1:ss(3))'-ss(3)/2;
intall=[];
 kapprox=1*pi*4/max(f);
for k=1:ss(4)
    inth=squeeze(sum(sum(allPSFs(range,range,:,k),1),2));
    intall(:,k)=inth;
end
normn=sum(intall,2)/2;   % normalize by i1+i2+i3+i4
intnf=intall./normn;
intn=intnf(frange,:);
%k phi1 phi2 As Bs...
    lb1=[(max(intn(:))-min(intn(:)))/4 -inf -inf  min(intn(:)) -inf -inf ];
    ub1=[max(intn(:)) inf 0  max(intn(:)) inf 0 ];
    s1=[(max(intn(:))-min(intn(:)))/2 0 0 0.5 0 0];
    [~,indmax]=max(intn(:,1));
     phasestart1=pi/2-indmax*kapprox+pi; if phasestart1<0,phasestart1=phasestart1+2*pi;end;

lba=horzcat(0,-pi,-pi,lb1,lb1,lb1,lb1);
uba=horzcat(inf,3*pi,3*pi,ub1,ub1,ub1,ub1);

startpa=[kapprox,phasestart1, phasestart1+pi/2,s1,s1,s1,s1];
fitpg=lsqcurvefit(@zintpg,startpa,f(frange),intn,lba,uba);
fitted=zintpg(fitpg,f(frange));
fitted=reshape(fitted,length(frange),4);

hold (ax,'off')
plot(ax,f,intnf,'-+')
% plot(ax,f(frange),intn)
hold(ax,'on')
% fst=zintpg(startpa,f(frange));
% plot(ax,f(frange),fst(:,:))
plot(ax,f(frange)',fitted','k');
phaseshiftso=fitpg([2 3]);
frequencyo=fitpg(1)/2;
title(ax,['frequency: ' num2str(frequencyo,2) ', phaseshift/pi: ' num2str((phaseshiftso(2)-phaseshiftso(1))/pi,2)])
%   fnc=@(k,phi1,phi2,A11,A21,A31,B11,B21,B31,A12,A22,A32,B12,B22,B32,A13,A23,A33,B13,B23,B33,A14,A24,A34,B14,B24,B34,x) zintg(k,phi1,A11,A21,A31,B11,B21,B31,phi2,A12,A22,A32,B12,B22,B32,phi3,A13,A23,A33,B13,B23,B33,phi4,A14,A24,A34,B14,B24,B34,x);

end

function into=zintpg(p,xdat)

into=horzcat(zintp([p(1) p(2) p(4:9)],xdat),zintp([p(1) p(3) p(10:15)],xdat),zintp([p(1) p(2)+pi p(16:21)],xdat),zintp([p(1) p(3)+pi p(22:27)],xdat));
% into=vertcat(zint(x,k,phi1,A11,A21,A31,B11,B21,B31),zint(x,k,phi2,A12,A22,A32,B12,B22,B32),zint(x,k,phi3,A13,A23,A33,B13,B23,B33),zint(x,k,phi4,A14,A24,A34,B14,B24,B34));
end
function into=zint(f,k,phi,A1,A2,A3,B1,B2,B3)
Bg=B1+B2*f+B3*f.^2;
Am=A1+A2*f+A3*f.^2;
os=sin(k*f+phi);
into=Bg+Am.*os;
end

function into=zintp(p,f)
k=p(1);phi=p(2);A1=p(3);A2=p(4);A3=p(5);B1=p(6);B2=p(7);B3=p(8);
into=zint(f,k,phi,A1,A2,A3,B1,B2,B3);
% Bg=B1+B2*f+B3*f.^2;
% Am=A1+A2*f+A3*f.^2;
% os=sin(k*f+phi);
% into=Bg+Am.*os;
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
% figure(188)
% ax=gca;
% Aall(:,:,:,1)=A12;Aall(:,:,:,2)=A23;Aall(:,:,:,3)=A34;Aall(:,:,:,4)=A41;
% imageslicer(Aall,'Parent',ax)
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

function cspline=getsmoothspline(V,p)

if ~isfield(p,'pixelsize')
    pixelsizeh=100;
else
    pixelsizeh=p.pixelsize{1}(1);
end
lambdax=p.smoothxy/pixelsizeh/100000;
lambdaz=p.smoothz/p.dz*100;
lambda=[lambdax lambdax lambdaz];
b3_0t=bsarray(double(V),'lambda',lambda);
zhd=1:1:b3_0t.dataSize(3);
dxxhd=1;
[XX,YY,ZZ]=meshgrid(1:dxxhd:b3_0t.dataSize(1),1:dxxhd:b3_0t.dataSize(2),zhd);
corrPSFhdt = interp3_0(b3_0t,XX,YY,ZZ,0);
cspline = Spline3D_interp(corrPSFhdt);
end
