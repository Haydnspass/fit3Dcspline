function calibrate_4pi(p)
%for now now segmentation of calibration. 
ph=p;
f=figure('Name','Bead calibration');
tg=uitabgroup(f);
t1=uitab(tg,'Title','prefit');
tgprefit=uitabgroup(t1);
ph.tabgroup=   tgprefit;
ph.isglobalfit=false;
ph.outputfile={};
    
    
[beads,ph]=images2beads_globalfit(ph);


sstack=size(beads(1).stack.image);


ybeads=getFieldAsVectorInd(beads,'pos',2);
y4pi=p.smappos.settings.y4pi;
x4pi=p.smappos.settings.x4pi;
height4pi=p.smappos.settings.height4pi; 
width4pi=p.smappos.settings.width4pi; 
mirror4pi=p.smappos.settings.mirror4pi;
fw=round((p.zcorrframes-1)/2);
framerange=round(sstack(3)/2-fw:sstack(3)/2+fw);
numchannels=length(y4pi);
allPSFs=zeros(sstack(1),sstack(2),sstack(3),numchannels);
for k=1:numchannels
    indbh=ybeads>=y4pi(k)&ybeads<=y4pi(k)+height4pi;
    
    xposh=getFieldAsVectorInd(beads(indbh),'pos',1)';
    yposh=getFieldAsVectorInd(beads(indbh),'pos',2)';
%     ph.tabgroup=uitabgroup(uitab(tgprefit,'Title',num2str(k)));
    [allstacks,filenumber]=bead2stack(beads(indbh));
    [allPSFs(:,:,:,k),shiftedstack,shift,indgood]=registerPSF3D_g(allstacks,[],struct('framerange',framerange));
    
%     [splinefit,indgood,posbeads,shift]=getstackcal_g(beads(indbh),ph);
    
    beadtrue{k}(:,1)=xposh(indgood)-shift(indgood,2); %experimentally: this works :)
    beadtrue{k}(:,2)=yposh(indgood)-shift(indgood,1);
    beadtrue{k}(:,3)=shift(indgood,3); %this is not yet tested, could be minus
    
end
%calculate transformN
transform=interfaces.LocTransformN;
pt.mirror=mirror4pi(1)*2; %2 for y coordinate
pt.xrange=[x4pi(1) x4pi(1)+width4pi];
pt.yrange=[y4pi(1) y4pi(1)+height4pi];
pt.units='pixel';
pt.type='projective';
transform.setTransform(1,pt)
for k=2:length(beadtrue)
    pt.mirror=mirror4pi(k)*2;
    pt.xrange=[x4pi(k) x4pi(k)+width4pi];
    pt.yrange=[y4pi(k) y4pi(k)+height4pi];
    transform.setTransform(k,pt)
    tab=(uitab(tgprefit,'Title',['T' num2str(k)]));ph.ax=axes(tab);
    transform_locs_simpleN(transform,1, beadtrue{1},k,beadtrue{k},ph)
end



%furhter align PSFs. Maybe this is not needed.
framerange=round(sstack(3)/2-2*fw:sstack(3)/2+2*fw);
[~,PSFaligned,shift,indgood]=registerPSF3D_g(allPSFs,[],struct('framerange',framerange,'removeoutliers',false),{},filenumber);

tab=(uitab(tgprefit,'Title','frequency'));ph.ax=axes(tab);
[phaseshifts,frequency]=getphaseshifts(allPSFs,ph.ax);
phaseshifts=phaseshifts-phaseshifts(1);
[I,A,B]=make4Pimodel(allPSFs,phaseshifts,frequency);




if ~p.isglobalfit %only normal local calibration, just call old proram
    calibrate3D_g(p);
    return
end

global S beadpos1 beadpos2
ph=p;
ph.isglobalfit=false;
%set spatial calibration
ph.outputfile={};
splitpos=256;% later to GUI?

switch p.Tmode
    case {'up-down','up-down mirror'}
        if max(p.yrange)<splitpos %defined only in upper part
            yrange1=p.yrange;
            yrange2=p.yrange+splitpos;yrange2(yrange2<splitpos)=splitpos;
        else
            yrange1=([p.yrange splitpos]);yrange1(yrange1>splitpos)=splitpos;yrange1=unique(yrange1);
            yrange2=([p.yrange+ splitpos]);yrange2(yrange2>splitpos)=splitpos;yrange2=unique(yrange2);
        end
            
%         yrange=unique([p.yrange splitpos p.yrange+splitpos]);
%         yrange1=yrange(yrange<=splitpos);yrange2=yrange(yrange>=splitpos);
        xrange1=p.xrange;xrange2=p.xrange;
    case {'right-left','right-left mirror'}
        xrange=unique([p.xrange splitpos p.xrange+splitpos]);  
        xrange1=xrange(xrange<=splitpos);xrange2=xrange(xrange>=splitpos);
        yrange1=p.xrange;yrange2=p.xrange;
end


f=figure('Name','Bead calibration');
tg=uitabgroup(f);
t1=uitab(tg,'Title','first channel');
ph.tabgroup=  uitabgroup(t1);  
    
ph.yrange=yrange1;ph.xrange=xrange1;
[S1,beadpos1,parameters1]=calibrate3D_g(ph);


t2=uitab(tg,'Title','second channel');
ph.tabgroup=  uitabgroup(t2);  

ph.yrange=yrange2;ph.xrange=xrange2;
[S2,beadpos2,parameters2]=calibrate3D_g(ph);


yrangeall=[yrange1(1:end-1) splitpos yrange2(2:end)];
S1.Yrangeall=yrangeall;
S2.Yrangeall=yrangeall;
S2.posind=[1,2];
% [S,beadpos]=calibrate3D_g(ph);

% Later: also do test-fitting with corresponding spline coefficients
p.tabgroup=uitab(tg,'Title','transformation');
% find transform
if p.makeT || isempty(p.Tfile)
    transform=transform_locs_simple(beadpos1{1},beadpos2{1},p);
else
    transform=load(p.Tfile);
end

ph=p;
ph.outputfile=[];
ph.isglobalfit=true;
ph.Tfile=transform;
ph.outputfile=p.outputfile;
t4=uitab(tg,'Title','global cal');
ph.tabgroup=  uitabgroup(t4);  
ph.yrange=yrange1;ph.xrange=xrange1;
[S,beadpos,parameters_g]=calibrate3D_g(ph);

SXY(1:length(S1))=S1;
SXY(end+1:end+length(S2))=S2;
SXY_g=S;
transformation=parameters_g.transformation;

if ~isempty(p.outputfile)
    if p.smap
        parameters1.smappos.P=[]; parameters2.smappos.P=[]; parameters_g.smappos.P=[];
        save(p.outputfile,'SXY','SXY_g','parameters_g','parameters1','parameters2','transformation');
    else
        save(p.outputfile,'gausscal','cspline_all','gauss_sx2_sy2','gauss_zfit','cspline','parameters');
    end
end
end


function [stack,filenumber]=bead2stack(beads)
ss=size(beads(1).stack.image);
stack=zeros(ss(1),ss(2),ss(3),length(beads));
for k=length(beads):-1:1
    stack(:,:,:,k)=beads(k).stack.image;
    filenumber(k)=beads.filenumber;
end
end

function [phaseshiftso,frequencyo]=getphaseshifts(allPSFs,ax)
ss=size(allPSFs);
range=(ss(1)+1)/2+[-1 1];
fw=20;
frange=round(ss(3)/2+(-fw:fw)');
% figure(88);
f=(1:ss(3))'-ss(3)/2;
intall=[];frall=[];
startpa=0;
   kapprox=1*pi*4/max(f);
   col=lines(ss(4));
for k=1:ss(4)
    inth=squeeze(sum(sum(allPSFs(range,range,:,k),1),2));

    lb=[0 -pi (max(inth)-min(inth))/4 -inf -inf  min(inth) -inf -inf ];
    ub=[inf 3*pi max(inth) inf 0  max(inth) inf 0 ];
%     fnc=@(k,phi,A1,A2,A3,B1,B2,B3,x) zint(x,k,phi,A1,A2,A3,B1,B2,B3);


    %approximate phase from max
 
    [~,indmax]=max(inth(f>0));
    phasestart=pi/2-indmax*kapprox; if phasestart<0,phasestart=phasestart+2*pi;end;
    
    startpar=[1*pi*4/max(f),phasestart,max(inth)/2,0,-max(inth)/max(f)^2/10,mean(inth),max(inth)/max(f)*0,-max(inth)/max(f)^2/10];
    fitp=lsqcurvefit(@zintp,startpar,f(frange),inth(frange),lb,ub);
  
    intall=vertcat(intall,inth(frange));
    frall=vertcat(frall,f(frange));
    startpa=horzcat(startpa,fitp(2:end));
%     fr=fit(f(frange),inth(frange),fnc,'StartPoint',startpar);
   
    
    fh=f(frange);
    plot(ax,f,inth,'-x','Color',col(k,:))
    hold(ax,'on')
     plot(ax,fh,zintp(fitp,fh),'r')
%      plot(fh,fitp(3)+fitp(4)*fh+fitp(5)*fh.^2,'g');
%      plot(fh,fitp(6)+fitp(7)*fh+fitp(8)*fh.^2,'c');
      plot(fh,fitp(3)+fitp(4)*fh+fitp(5)*fh.^2 + fitp(6)+fitp(7)*fh+fitp(8)*fh.^2,'Color',col(k,:));
      plot(fh,-(fitp(3)+fitp(4)*fh+fitp(5)*fh.^2) + fitp(6)+fitp(7)*fh+fitp(8)*fh.^2,'Color',col(k,:));
%     plot(ax,fh,fr.A1+fr.A2*fh+fr.A3*fh.^2)
   
    frequency(k)=fitp(1)/2;
    phaseshifts(k)=fitp(2);
    kapprox=fitp(1);
    
end
startpa(1)=fitp(1);

%global fit with same k
%XXXXX later:consider linking also phi1 and phi3 and phi2,phi4 to sum up to
%pi.

% startparall=horzcat(startpar,startpar(2:end),startpar(2:end),startpar(2:end));
lba=horzcat(lb,lb(2:end),lb(2:end),lb(2:end));
uba=horzcat(ub,ub(2:end),ub(2:end),ub(2:end));

fitpg=lsqcurvefit(@zintpg,startpa,f(frange),intall,lba,uba);
fitted=zintpg(fitpg,f(frange));
fitted=reshape(fitted,length(frange),4);
plot(ax,f(frange)',fitted','k');
phaseshiftso=fitpg([2 9 16 23]);
frequencyo=fitpg(1)/2;
%   fnc=@(k,phi1,A11,A21,A31,B11,B21,B31,phi2,A12,A22,A32,B12,B22,B32,phi3,A13,A23,A33,B13,B23,B33,phi4,A14,A24,A34,B14,B24,B34,x) zintg(k,phi1,A11,A21,A31,B11,B21,B31,phi2,A12,A22,A32,B12,B22,B32,phi3,A13,A23,A33,B13,B23,B33,phi4,A14,A24,A34,B14,B24,B34,x);
% 
%  fr=fit(frall,intall,fnc,'StartPoint',startparall);
end

% function into=zintg(k,phi1,A11,A21,A31,B11,B21,B31,phi2,A12,A22,A32,B12,B22,B32,phi3,A13,A23,A33,B13,B23,B33,phi4,A14,A24,A34,B14,B24,B34,x)
% into=vertcat(zint(x,k,phi1,A11,A21,A31,B11,B21,B31),zint(x,k,phi2,A12,A22,A32,B12,B22,B32),zint(x,k,phi3,A13,A23,A33,B13,B23,B33),zint(x,k,phi4,A14,A24,A34,B14,B24,B34));
% end
function into=zintpg(p,xdat)

into=vertcat(zintp([p(1) p(2:8)],xdat),zintp([p(1) p(9:15)],xdat),zintp([p(1) p(16:22)],xdat),zintp([p(1) p(23:29)],xdat));
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
