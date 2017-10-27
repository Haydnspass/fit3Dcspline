function [splinefit,indgood]=getstackcal_g(beads,p)
global stackcal_testfit
isastig=contains(p.modality,'astig')||contains(p.modality,'2D');
alignzastig=isastig&contains(p.zcorr,'astig');
zcorr=contains(p.zcorr,'corr');
sstack=size(beads(1).stack.image);
    halfstoreframes=round((size(beads(1).stack.image,3)-1)/2);
    if isastig    
        for B=length(beads):-1:1
            if  halfstoreframes<length(beads(B).stack.framerange)
                dframe(B)=beads(B).stack.framerange(halfstoreframes+1)-beads(B).f0;
            else
                dframe(B)=NaN;
            end
        end
        
    %remove outliers:
        badind=abs(dframe-nanmedian(dframe))>10|isnan(dframe);
        beads(badind)=[];
    

        psfx=[beads(:).psfx0];psfy=[beads(:).psfy0];
        dpsfx=(psfx-median(psfx(~isnan(psfx))))*10;
        dpsfy=(psfy-median(psfy(~isnan(psfy))))*10;
    else
        dframe=0;
        dpsfx=0;dpsfy=0;
    end
    
    allstacks=zeros(sstack(1),sstack(2),sstack(3),length(beads))+NaN;
    allstackst=zeros(sstack(1),sstack(2),sstack(3),length(beads))+NaN;
    goodvs=[];
    for B=length(beads):-1:1
        allstacks(:,:,:,B)=beads(B).stack.image;
        allstackst(:,:,:,B)=beads(B).stack.imagetar;
        stackh=allstacks(:,:,:,B);
        goodvs(B)=sum(~isnan(stackh(:)))/numel(stackh);
        shiftxy(B,1:2)=beads(B).shiftxy;
    end
    
    mstack=nanmean(allstacks,4);
    mstack=mstack-nanmin(mstack(:));
    mstack=mstack/nansum(mstack(:));
    for k=length(beads):-1:1
    	stackh=(allstacks(:,:,:,k));
        stackh=stackh-nanmin(stackh(:));
        stackh=stackh/nansum(stackh(:));
        dstack(k)=sum((stackh(:)-mstack(:)).^2);
    end
    dstack=dstack/mean(dstack);  
    
    mstack=nanmean(allstackst,4);
    mstack=mstack-nanmin(mstack(:));
    mstack=mstack/nansum(mstack(:));
    for k=length(beads):-1:1
    	stackh=(allstackst(:,:,:,k));
        stackh=stackh-nanmin(stackh(:));
        stackh=stackh/nansum(stackh(:));
        dstackt(k)=sum((stackh(:)-mstack(:)).^2);
    end
    dstackt=dstackt/mean(dstackt);  
    
    devs=(dpsfx.^2+dpsfy.^2+dstack+dstackt)./goodvs;

    if zcorr
        
        fw2=round((p.zcorrframes-1)/2);
        
    else
        fw2=2;
    end
   
%     ax=axes('Parent',uitab(p.tabgroup,'Title','scatter'));

    [~,sortinddev]=sort(devs);
    allrois=allstacks(:,:,:,sortinddev);
    allroist=allstackst(:,:,:,sortinddev);
    shiftxys=shiftxy(sortinddev,:);
    if alignzastig
        zshift=dframe(sortinddev)-round(median(dframe));
    else
        zshift=[];
    end
    
%     focusreference=round(median(dframe));
    midrange=halfstoreframes+1-round(median(dframe));
     framerange=max(1,midrange-fw2):min(midrange+fw2,size(stackh,3));
    p.status.String='calculate shift of individual PSFs';drawnow
    filenumber=[beads(:).filenumber];
    [corrPSF,shiftedstack,shift,beadgood]=registerPSF3D_g(allrois,allroist,struct('shiftxy',shiftxys,'framerange',framerange,'alignz',zcorr,'zshiftf0',zshift,'beadfilterf0',false),{},filenumber(sortinddev));
    
    corrPSFr=corrPSF(1:size(allrois,1),:,:);
    corrPSFt=corrPSF(size(allrois,1)+1:end,:,:);
    %undo sorting by deviation to associate beads again to their
    %bead number
    [~,sortback]=sort(sortinddev);
    shiftedstack=shiftedstack(:,:,:,sortback);
    beadgood=beadgood(sortback);
%     shiftxys=shiftxys(sortback,:);
    indgood=beadgood;
    allrois=allstacks;
  

        %cut out the central part of the PSF correspoinding to the set
        %Roisize in x,y and z

        scorrPSF=size(corrPSFr);
        x=round((scorrPSF(1)+1)/2);y=round((scorrPSF(2)+1)/2);

        dRx=round((p.ROIxy-1)/2);
        if isnan(p.ROIz)
            p.ROIz=size(corrPSFr,3);
        end
            dzroi=round((p.ROIz-1)/2);
        
        rangex=x-dRx:x+dRx;
        rangey=y-dRx:y+dRx;

        z=midrange;%always same reference: z=f0
        rangez=max(1,z-dzroi):min(size(corrPSFr,3),z+dzroi);
        z0reference=find(rangez>=z,1,'first');
        
        %careful: like this each PSF is normalized individually. This might
        %not be the right approahc. Then normalize by one only
        %normalize PSF
        centpsfr=corrPSFr(rangex,rangey,z-1:z+1); %cut out rim from shift
         centpsft=corrPSFt(rangex,rangey,z-1:z+1);
%         centpsf=corrPSF(2:end-1,2:end-1,2:end-1); %cut out rim from shift
        minPSFr=nanmin(centpsfr(:));minPSFt=nanmin(centpsft(:));
        corrPSFnr=corrPSFr-minPSFr;corrPSFnt=corrPSFt-minPSFt;
%         corrPSFn=corrPSF;
        intglobalr=nanmean(nansum(nansum(corrPSFnr(rangex,rangey,z-1:z+1),1),2));
        intglobalt=nanmean(nansum(nansum(corrPSFnt(rangex,rangey,z-1:z+1),1),2));
        corrPSFnr=corrPSFnr/intglobalr;corrPSFnt=corrPSFnt/intglobalt;

        shiftedstack(1:size(allrois,1),:,:,:)=shiftedstack(1:size(allrois,1),:,:,:)/intglobalr;
        corrPSFnr(isnan(corrPSFnr))=0;
        corrPSFnr(corrPSFnr<0)=0;
        corrPSFsr=corrPSFnr(rangex,rangey,rangez);
        shiftedstack(size(allrois,1)+1:end,:,:,:)=shiftedstack(size(allrois,1)+1:end,:,:,:)/intglobalt;        
        corrPSFnt(isnan(corrPSFnt))=0;
        corrPSFnt(corrPSFnt<0)=0;
        corrPSFst=corrPSFnt(rangex,rangey,rangez);
        
        PSFgood=true;

        %calculate effective smoothing factor. For dz=10 nm, pixelsize= 130
        %nm, a value around 1 produces visible but not too much smoothing.
        lambdax=p.smoothxy/p.cam_pixelsize_um(1)/100000;
        lambdaz=p.smoothz/p.dz*100;
        lambda=[lambdax lambdax lambdaz];
        %calculate smoothed bsplines
        b3_0r=bsarray(double(corrPSFsr),'lambda',lambda);
        b3_0t=bsarray(double(corrPSFst),'lambda',lambda);

        %calculate smoothed volume
        zhd=1:1:b3_0r.dataSize(3);
        dxxhd=1;
        [XX,YY,ZZ]=meshgrid(1:dxxhd:b3_0r.dataSize(1),1:dxxhd:b3_0r.dataSize(2),zhd);
        p.status.String='calculating cspline coefficients in progress';drawnow
        corrPSFhdr = interp3_0(b3_0r,XX,YY,ZZ,0);
        corrPSFhdt = interp3_0(b3_0t,XX,YY,ZZ,0);
        
        %calculate cspline coefficients
%         spline = Spline3D_v2(corrPSFhd);
%         coeff = spline.coeff;
        coeffr = Spline3D_interp(corrPSFhdr);
        coefft = Spline3D_interp(corrPSFhdt);
       
        %assemble output structure for saving
        bspline.bslpine={b3_0r,b3_0t};
        cspline.coeff={coeffr, coefft};
        cspline.z0=z0reference;%round((b3_0.dataSize(3)+1)/2);
        cspline.dz=p.dz;
        cspline.x0=dRx+1;
        bspline.z0=round((b3_0r.dataSize(3)+1)/2);
        bspline.dz=p.dz;            
        splinefit.bspline=bspline;
        p.z0=cspline.z0;
        
        splinefit.PSF={corrPSFr,corrPSFt};
        
        splinefit.PSFsmooth={corrPSFhdr,corrPSFhdt};
        splinefit.cspline=cspline;

        %plot graphs
        if PSFgood      
            dL=size(shiftedstack,1)/2;
            ax=axes(uitab(p.tabgroup,'Title','PSFz'));
             framerange0=max(p.fminmax(1)):min(p.fminmax(2));
             halfroisizebig=(size(shiftedstack,2)-1)/2;         
            ftest=z;
            xt=x;
            yt=y;
            zpallr=squeeze(shiftedstack(xt,yt,:,beadgood));
            zpallt=squeeze(shiftedstack(xt+dL,yt,:,beadgood));
%             zpall2=squeeze(allrois(xt,yt,:,beadgood));
            xpall=squeeze(shiftedstack(:,yt,ftest,beadgood));
%             xpall2=squeeze(allrois(:,yt,ftest,beadgood));
%             for k=1:size(zpall,2)
%                 zpall2(:,k)=zpall2(:,k)/nanmax(zpall2(:,k));
%                 xpall2(:,k)=xpall2(:,k)/nanmax(xpall2(:,k));                
%             end           
            zprofiler=squeeze(corrPSFnr(xt,yt,:));
            zprofilet=squeeze(corrPSFnt(xt,yt,:));
%             mphd=round((size(corrPSFhd,1)+1)/2);
                 
            xprofiler=squeeze(corrPSFnr(:,yt,ftest));
             xprofilet=squeeze(corrPSFnt(:,yt,ftest));
            mpzhd=round((size(corrPSFhdr,3)+1)/2+1);
            dzzz=round((size(corrPSFnr,3)-1)/2+1)-mpzhd;
            dxxx=0.1;
            xxx=1:dxxx:b3_0r.dataSize(1);
%             zzzt=0*xxx+mpzhd+dzzz-1;
            zzzt=0*xxx+ftest;
            xbsr= interp3_0(b3_0r,0*xxx+b3_0r.dataSize(1)/2+.5,xxx,zzzt);
            xbst= interp3_0(b3_0t,0*xxx+b3_0t.dataSize(1)/2+.5,xxx,zzzt);
            zzz=1:dxxx:b3_0r.dataSize(3);xxxt=0*zzz+b3_0r.dataSize(1)/2+.5;
            zbsr= interp3_0(b3_0r,xxxt,xxxt,zzz); 
            zbst= interp3_0(b3_0t,xxxt,xxxt,zzz); 
            hold(ax,'off')
             plot(ax,framerange0,zpallr(1:length(framerange0),:),'c')
             hold(ax,'on')
             dF=max(framerange0);
             plot(ax,framerange0+dF,zpallt(1:length(framerange0),:),'c')
            plot(ax,framerange0',zprofiler(1:length(framerange0)),'k*')
            plot(ax,framerange0'+dF,zprofilet(1:length(framerange0)),'k*')
            plot(ax,zzz+rangez(1)+framerange0(1)-2,zbsr,'b','LineWidth',2)
             plot(ax,zzz+rangez(1)+framerange0(1)-2+dF,zbst,'b','LineWidth',2)
            xlabel(ax,'frames')
            ylabel(ax,'normalized intensity')
            ax.XLim(2)=max(framerange0+dF);ax.XLim(1)=min(framerange0);
            title(ax,'Profile along z for x=0, y=0');
            
            legend('individual PSFs','average PSF','smoothed spline')
            
            xrange=-halfroisizebig:halfroisizebig;
             ax=axes(uitab(p.tabgroup,'Title','PSFx'));
            hold(ax,'off')
            plot(ax,[xrange xrange+dL],xpall,'c')
            hold(ax,'on')
            plot(ax,[xrange xrange+dL],vertcat(xprofiler,xprofilet),'k*-')
            plot(ax,(xxx-(b3_0r.dataSize(1)+1)/2),xbsr,'b','LineWidth',2)
            plot(ax,(xxx-(b3_0r.dataSize(1)+1)/2)+dL,xbst,'b','LineWidth',2)
            xlabel(ax,'x (pixel)')
            ylabel(ax,'normalized intensity')
            title(ax,'Profile along x for y=0, z=0');
             legend('individual PSFs','average PSF','smoothed spline')
            
            drawnow
            
            %quality control: refit all beads
            if isempty(stackcal_testfit)||stackcal_testfit  %not implemented yet in fitter. Fix later
                ax=axes(uitab(p.tabgroup,'Title','validate'));
                testallrois(:,:,:,:,1)=allrois(:,:,:,beadgood);
                testallrois(:,:,:,:,2)=allroist(:,:,:,beadgood);
                
                corrPSFfit=corrPSF/max(corrPSF(:))*max(testallrois(:)); %bring back to some reasonable photon numbers;
                corrPSFfitf(:,:,:,1,1)=corrPSFfit(1:size(corrPSFfit,2),:,:);
                corrPSFfitf(:,:,:,1,2)=corrPSFfit(size(corrPSFfit,2)+1:end,:,:);
                zref=testfit(corrPSFfitf,cspline.coeff,[0 0],p,{'k','LineWidth',2},ax);
                

                
                testallrois(isnan(testallrois))=0;
                zall=testfit(testallrois,cspline.coeff,-shiftxy,p,{},ax);
                drawnow
            end
        end 
end

function zs=testfit(teststack,coeff,shiftxy, p,linepar,ax)
if nargin<4
    linepar={};
elseif ~iscell(linepar)
    linepar={linepar};
end
d=round((size(teststack,1)-p.ROIxy)/2);
            range=d+1:d+p.ROIxy;

numstack=size(teststack,4);
t=tic;
% f=figure(989);ax2=gca;hold off
    for k=1:size(teststack,4)
        if toc(t)>1
            p.status.String=['fitting test stacks: ' num2str(k/numstack,'%1.2f')];drawnow
            t=tic;
        end
        if contains(p.modality,'2D')
            fitmode=6;
        else
            fitmode=5;
        end
        fitstack=single(squeeze(teststack(range,range,:,k,:)));
        coeffh(:,:,:,:,1)=single(coeff{1});
        coeffh(:,:,:,:,2)=single(coeff{2});
        shared=[1 1 1 0 0];
        nfits=size(fitstack,3);
        npar=5;
        dT=zeros(npar,2,nfits);
        dT(1,2,:)=shiftxy(k,1);
        dT(2,2,:)=shiftxy(k,2);
%         [PM,CRLBM, LLM,update, error] =  kernel_MLEfit_Spline_LM_multichannel_finalized(fitstack,coeffh, shared,dT,50);
        sharedA = repmat(shared,[1 size(fitstack,3)]);
        [P,CRLB, LL] =CPUmleFit_LM_MultiChannel(fitstack,int32(sharedA),50,coeffh,single(dT));
        
    
      
%         [P] =  mleFit_LM(single(squeeze(teststack(range,range,:,k))),fitmode,100,single(coeff),0,1);
        
        z=(1:size(P,1))'-1;

        znm=(P(:,3)-p.z0)*p.dz;
        plot(ax,z,znm,linepar{:})
        hold(ax,'on')
        xlabel(ax,'frame')
        ylabel(ax,'zfit (nm)')
        zs(:,k)=P(:,5);
        
        
        if 0% imageslicer to test
            coord=P(:,1:5);
            coord2=coord;
            coord2(:,1)=coord2(:,1)+squeeze(dT(1,2,:));
            coord2(:,2)=coord2(:,2)+squeeze(dT(2,2,:));
            img1=renderPSF(coeff{1},coord,size(fitstack,1));
            img2=renderPSF(coeff{2},coord2,size(fitstack,1));
            imall=[fitstack(:,:,:,1),img1; fitstack(:,:,:,2),img2];
            res=[fitstack(:,:,:,1)-img1, 0*img1;fitstack(:,:,:,2)-img2,0*img2];
            ims(:,:,:,1)=imall;ims(:,:,:,2)=res;
            imageslicer(ims)
        end
        
% test for the returned photons and photons in the raw image        
%         phot=P(:,3); bg=P(:,4);
%         totsum=squeeze(nansum( nansum(teststack(range,range,:,k),1),2));
%         totsum=totsum-squeeze(min(min(teststack(range,range,:,k),[],1),[],2))*length(range)^2;
%         photsum=phot+0*bg*length(range)^2;
%         plot(ax2,z,(photsum-totsum)./totsum,'.')
%         hold(ax2,'on')
    end
    
end

function teststripes(coeff,p,ax)
%not used, can be called to test for stripe artifacts.
tt=tic;

zr=0:0.2:p.ROIz;
xr=0:0.05:p.ROIxy;
hz=zeros(1,length(zr)-1);
hx=zeros(1,length(xr)-1);
hy=hx;
while toc(tt)<30
    nn=rand(11,11,10000,'single');
%     P=callYimingFitter(nn,single(coeff),50,5,0,1);
    [P] =  mleFit_LM(nn,5,50,single(coeff),0,1);
    hz=histcounts(P(:,5),zr)+hz;
    hx=histcounts(P(:,1),xr)+hx;
    hy=histcounts(P(:,2),xr)+hy;
    
end

hz(1)=[];hz(end)=[];
hz(1)=0;hz(end)=0;

indx=(hx==0);
hx(indx)=[];
indy=(hy==0);
hy(indy)=[];
hx(1)=[];hx(end)=[];
hy(1)=[];hy(end)=[];
hzx=myxcorr(hz-mean(hz),hz-mean(hz));
hxx=myxcorr(hx-mean(hx),hx-mean(hx));
hyx=myxcorr(hy-mean(hy),hy-mean(hy));
hzx(1)=0;hzx(end)=0;
ax2=axes(ax.Parent);
subplot(1,2,1,ax);
subplot(1,2,2,ax2);
findx=find(~indx);findy=find(~indy);
plot(ax,zr(2:end-2),hz,zr(2:end-2),hzx/max(hzx)*(quantile(hz,.99)));
ax.YLim(2)=(quantile(hz,.99))*1.1;
ax.YLim(1)=min(quantile(hz,.01),quantile(hzx/max(hzx)*(quantile(hz,.99)),.01));
plot(ax2,xr(findx(2:end-1)),hx,xr(findx(2:end-1)),hxx/max(hxx)*max(hx),xr(findy(2:end-1)),hy,xr(findy(2:end-1)),hyx/max(hyx)*max(hy));
end
