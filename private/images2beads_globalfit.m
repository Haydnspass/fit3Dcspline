function [b,p]=images2beads_globalfit(p)
% addpath('bfmatlab')
fs=p.filtersize;
h=fspecial('gaussian',2*round(fs*3/2)+1,fs);
fmax=0;
roisize=p.ROIxy;
roisizeh=round(1.5*(p.ROIxy-1)/2); %create extra space if we need to shift;
rsr=-roisizeh:roisizeh;
filelist=p.filelist;
b=[];
ht=uitab(p.tabgroup,'Title','Files');
tg=uitabgroup(ht);
l=load(p.Tfile);
if p.isglobalfit
    transform=l.transformation;
    p.transformation=transform;
    p.mirror=contains(transform.tinfo.mirror.targetmirror,'up-down');
else
    p.mirror=false;
end
for k=1:length(filelist)
    ax=axes(uitab(tg,'Title',num2str(k)));
    p.fileax(k)=ax;
    if isfield(p,'smap') && p.smap
        
        imstack=readfile_ome(filelist{k});
        if isempty(imstack)
            disp('using simple reader')
            imstack=readfile_tif(filelist{k});
        end
          
    else
        imstack=readfile_tif(filelist{k});
    end
    
    if p.emgain
        imstack=imstack(:,end:-1:1);
    end
       
    if isfield(p,'framerangeuse')
        imstack=imstack(:,:,p.framerangeuse(1):p.framerangeuse(end));
    end
    
    imstack=imstack-min(imstack(:)); %fast fix for offset;
    
%     imageslicer(imstack)%%%%%%%%XXXXXXX
    
    mim=max(imstack,[],3);
    mim=filter2(h,mim);
    imagesc(ax,mim);
    axis(ax,'equal');
    axis(ax,'off')
    title(ax,'Maximum intensity projection')
    if isfield(p,'beadpos') %passed on externally
        maxima=round(p.beadpos{k});
    else
        maxima=maximumfindcall(mim);
        int=maxima(:,3);
        try
        mimc=mim(roisize:end-roisize,roisize:end-roisize);
        mmed=quantile(mimc(:),0.3);
        imt=mimc(mimc<mmed);
            sm=sort(int);
        mv=mean(sm(end-5:end));
%         cutoff=mean(imt(:))+max(2.5*std(imt(:)),(mv-mean(imt(:)))/10);
        cutoff=mean(imt(:))+max(2.5*std(imt(:)),(mv-mean(imt(:)))/15);
        catch
            cutoff=quantile(mimc(:),.95);
        end
        if any(int>cutoff)
            maxima=maxima(int>cutoff,:);
        else
            [~,indm]=max(int);
            maxima=maxima(indm,:);
        end
    end
    
    %remove beads that are closer together than mindistance
    if isfield(p,'mindistance')&&~isempty(p.mindistance)
        indgoodb=true(size(maxima,1),1);
        for bk=1:size(maxima,1)
            for bl=bk+1:size(maxima,1)
                if  sum((maxima(bk,1:2)-maxima(bl,1:2)).^2)<p.mindistance^2
                    indgoodb(bk)=false;
                    indgoodb(bl)=false;
                end
            end
        end 
        maxima=maxima(indgoodb,:);
    end 
    
    hold (ax,'on')
    if p.isglobalfit
        %calculate in nm on chip (reference for transformation)
        maximanm=(maxima(:,1:2)+p.smappos.roi{k}([1 2]));
        maximanm(:,1)=maximanm(:,1)*p.smappos.pixelsize{k}(1)*1000;
        maximanm(:,2)=maximanm(:,2)*p.smappos.pixelsize{k}(end)*1000;

        %transform reference to target

        indref=transform.getRef(maximanm(:,1),maximanm(:,2));
        maximaref=maxima(indref,:);
        [x,y]=transform.transformCoordinatesFwd(maximanm(indref,1),maximanm(indref,2));
    %     [x,y]=transform.transformCoordinatesFwd(maximanm(indref,2),maximanm(indref,1));

        maximatargetf=[];
        maximatargetf(:,1)=x/p.smappos.pixelsize{k}(1)/1000-p.smappos.roi{k}(1);
        maximatargetf(:,2)=y/p.smappos.pixelsize{k}(end)/1000-p.smappos.roi{k}(2);
    %     maximatargetf(:,2)=x/p.smappos.pixelsize{k}(1)/1000-p.smappos.roi{k}(1);
    %     maximatargetf(:,1)=y/p.smappos.pixelsize{k}(end)/1000-p.smappos.roi{k}(2);


        if 0 %for testing
            maximatargetf(:,1)=maximatargetf(:,1)+1;
            maximatargetf(:,2)=maximatargetf(:,2)+0.5;
        maximatargetfm(:,1)=maximatargetf(:,1)-0.1+2;
        maximatargetfm(:,2)=maximatargetf(:,2)+0.1;
        maximatar=round(maximatargetfm);
        else 
        maximatar=round(maximatargetf);
        end
        dxy=maximatargetf-maximatar;       
        plot(ax,maximaref(:,1),maximaref(:,2),'ko',maximatar(:,1),maximatar(:,2),'kd')   
    else
        plot(ax,maxima(:,1),maxima(:,2),'ko')
        maximaref=maxima;
        maximatar=maxima;
        dxy=zeros(size(maximatar));
    end
    hold (ax,'off')
    drawnow
    numframes=size(imstack,3);
    bind=length(b)+size(maximaref,1);

    for l=1:size(maximaref,1)
        b(bind).loc.frames=(1:numframes)';
        b(bind).loc.filenumber=zeros(numframes,1)+k;
        b(bind).filenumber=k;
        b(bind).pos=maximaref(l,1:2);
        b(bind).postar=maximatar(l,1:2);
        b(bind).shiftxy=dxy(l,:);
        try
            b(bind).stack.image=imstack(b(bind).pos(2)+rsr,b(bind).pos(1)+rsr,:);
            b(bind).stack.imagetar=imstack(b(bind).postar(2)+rsr,b(bind).postar(1)+rsr,:);
            b(bind).stack.framerange=1:numframes;
            b(bind).isstack=true;
            
        catch err
            b(bind).isstack=false;
%             err
        end
        
            b(bind).roi=p.smappos.roi{k};
        bind=bind-1;
    end
    fmax=max(fmax,numframes);
end
b=b([b(:).isstack]);

p.fminmax=[1 fmax];

        if isfield(p,'files')&&~isempty(p.files)
            p.cam_pixelsize_um=p.files(k).info.cam_pixelsize_um;
        else
            p.cam_pixelsize_um=[1 1]/10; %?????
        end      

p.pathhere=fileparts(filelist{1});
end


