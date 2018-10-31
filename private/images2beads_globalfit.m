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
is4pi=contains(p.modality,'4Pi');

if p.isglobalfit 
   if ischar(p.Tfile)%extend later to 4Pi
        l=load(p.Tfile);
        transform=l.transformation;
   else
       transform=p.Tfile;
   end
    p.transformation=transform;
    p.mirror=contains(transform.tinfo.mirror.targetmirror,'up-down')|contains(transform.tinfo.mirror.targetmirror,'left-right');
else
    p.mirror=false;
end


for k=1:length(filelist)
    ax=axes(uitab(tg,'Title',num2str(k)));
%     axis(ax,'image')
    p.fileax(k)=ax;
    [imstack, p.roi{k}, p.pixelsize{k},settings3D]=readbeadimages(filelist{k},p);
    
    if is4pi
        if ~isempty(settings3D)
            p.settings_3D=settings3D;
            p.settings_3D.file='multifile';
        end
        if isfield(p,'settings_3D') && ~isempty(p.settings_3D) %calibration file: cut out and mirror already here!
            imstack=cutoutchannels(imstack,p.settings_3D);
            p.roi{k}=[0 0 size(imstack,1) size(imstack,2)]; %check x,y
        else
            disp('no  settings_3D found. Use default. Specify settings file?')
            wx=size(imstack,2)/4;wy=size(imstack,1);
            p.settings_3D=struct('y4pi',[0 0 0 0],'x4pi',[0 wx 2*wx 3*wx], 'width4pi',wx,'height4pi',wy,'mirror4pi',[0 0 0 0],'pixelsize_nm',100,'offset',100,'conversion',0.5);
        end
    end
    
    if p.emgain
        imstack=imstack(:,end:-1:1,:);
    end
       
    if isfield(p,'framerangeuse')
        imstack=imstack(:,:,p.framerangeuse(1):p.framerangeuse(end));
    end
    
    imstack=imstack-min(imstack(:)); %fast fix for offset;
    
%     imageslicer(imstack)%%%%%%%%XXXXXXX
    
    mim=max(imstack,[],3);

%     mim=mean(imstack,3);
    mim=filter2(h,mim);
    imagesc(ax,mim);
    axis(ax,'image');
    axis(ax,'off')
    title(ax,'Maximum intensity projection')
    if isfield(p,'beadpos') %passed on externally
        maxima=round(p.beadpos{k});
    else
        if isfield(p,'roimask')&&~isempty(p.roimask)
            mim=double(mim).*double(p.roimask);
        end
        maxima=maximumfindcall(mim);
        %XXXXXXXXXXX
%         maxima(:,1)=(maxima(:,1)+2*round(rand(size(maxima,1),1)))-1; % for testing if positions match
        int=maxima(:,3);
        try
        mimc=mim(roisize:end-roisize,roisize:end-roisize);
        mmed=myquantile(mimc(:),0.3);
        imt=mimc(mimc<=mmed);
            sm=sort(int);
        mv=mean(sm(end-5:end));
%         cutoff=mean(imt(:))+max(2.5*std(imt(:)),(mv-mean(imt(:)))/10);
        cutoff=mean(imt(:))+max(2.5*std(imt(:)),(mv-mean(imt(:)))/15);
        catch
            cutoff=myquantile(mimc(:),.95);
        end
        cutoff=cutoff*p.cutoffrel;
        if any(int>cutoff)
            maxima=maxima(int>cutoff,:);
        else
            [~,indm]=max(int);
            maxima=maxima(indm,:);
        end
    end
    indgoodb=true(size(maxima,1),1);
    %remove beads that are closer together than mindistance
    if isfield(p,'mindistance')&&~isempty(p.mindistance)
        
        for bk=1:size(maxima,1)
            for bl=bk+1:size(maxima,1)
                if  sum((maxima(bk,1:2)-maxima(bl,1:2)).^2)<p.mindistance^2
                    indgoodb(bk)=false;
                    indgoodb(bl)=false;
                end
            end
        end 
       if isfield(p,'settings_3D')
           w=p.settings_3D.width4pi;
           xm=mod(maxima(:,1),w);
           xm(xm>w/2)=xm(xm>w/2)-w;
           indgoodb(abs(xm)<p.mindistance/2)=false;
       end
       h=size(imstack,1);
       ym=maxima(:,2);
       ym(ym>h/2)=ym(ym>h/2)-h;
       indgoodb(abs(ym)<p.mindistance/2)=false;
%         maxima=maxima(indgoodb,:);
         maxima=maxima(indgoodb,:);
    end 
    
  
    if p.isglobalfit
        %calculate in nm on chip (reference for transformation)
        maximanm=(maxima(:,1:2)+p.roi{k}([1 2]));
        
        if isfield(transform.tinfo,'units') &&strcmp(transform.tinfo.units,'pixels')
            pixfac=[1 1];
        else
            pixfac=[p.pixelsize{k}(1)*1000 p.pixelsize{k}(end)*1000];
        end
            maximanm(:,1)=maximanm(:,1)*pixfac(1);
            maximanm(:,2)=maximanm(:,2)*pixfac(end);

        %transform reference to target

            indref=transform.getRef(maximanm(:,1),maximanm(:,2));
            maximaref=maxima(indref,:);
            [x,y]=transform.transformCoordinatesFwd(maximanm(indref,1),maximanm(indref,2));
        %     [x,y]=transform.transformCoordinatesFwd(maximanm(indref,2),maximanm(indref,1));

            maximatargetf=[];
            maximatargetf(:,1)=x/pixfac(1)-p.roi{k}(1);
            maximatargetf(:,2)=y/pixfac(end)-p.roi{k}(2);
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
          
    else
        indgoodr=maxima(:,1)>p.xrange(1)&maxima(:,1)<p.xrange(end)&maxima(:,2)>p.yrange(1)&maxima(:,2)<p.yrange(end);
        maxima=maxima(indgoodr,:);
        maximaref=maxima;
        maximatar=maxima;
        dxy=zeros(size(maximatar));
    end
   
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
        
            b(bind).roi=p.roi{k};
        bind=bind-1;
    end
    fmax=max(fmax,numframes);
     hold (ax,'on')
    if p.isglobalfit
        plot(ax,maximaref(:,1),maximaref(:,2),'ko',maximatar(:,1),maximatar(:,2),'kd') 
    else
        plot(ax,maxima(:,1),maxima(:,2),'ko')
    end
    hold (ax,'off')
    drawnow
end
indgoodbead=[b(:).isstack];
b=b(indgoodbead);

% plot

    
p.fminmax=[1 fmax];

%         if isfield(p,'files')&&~isempty(p.files)
%             p.cam_pixelsize_um=p.files(k).info.cam_pixelsize_um;
%         else
%             p.cam_pixelsize_um=[1 1]/10; %?????
%         end      

p.pathhere=fileparts(filelist{1});
end


