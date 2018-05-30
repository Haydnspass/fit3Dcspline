function [imstack, roi, pixelsize,settings3D]=readbeadimages(file,p)
[path,f,ext]=fileparts(file);
indq=strfind(f,'_q');
settings3D=[];
multichannel=false;
if ~isempty(indq)
    allfiles=dir([path filesep f(1:indq+1)  '*' ext]);
    for k=1:length(allfiles)
        files{k}=[path filesep allfiles(k).name];
        disp(allfiles(k).name(indq:end))
    end
    file=files;
    multichannel=true;
end
pixelsize=100;
    
    if isfield(p,'smap') && p.smap
        
        try
             r=imageloaderAll(file,[],p.smappos.P);

             imstack=r.getmanyimages(1:r.metadata.numberOfFrames,'mat');
             roi=r.metadata.roi;
             pixelsize=r.metadata.cam_pixelsize_um;
             r.close;
        catch err
            err
            imstack=readfile_tif(file);
            roi=[0 0 size(imstack,1) size(imstack,2)]; %check x,y
        end
        if isempty(imstack)
            disp('using simple reader')
            warndlg('using simple reader, this might create problems if only part of the camera chip is used.','replace');
            if multichannel
                imstack=[];
                for k=1:length(file)
                    imstack=vertcat(imstack,readfile_tif(file{k}));
                end
            else
                 imstack=readfile_tif(file);
            end
            roi=[0 0 size(imstack,1) size(imstack,2)];
        end
          
    else
        imstack=readfile_tif(file);
        roi=[0 0 size(imstack,1) size(imstack,2)]; %check x,y
    end
    if multichannel
        wx=size(imstack,1)/4;wy=size(imstack,2);
        settings3D=struct('x4pi',[0 0 0 0],'y4pi',[0 wy 2*wy 3*wy], 'width4pi',wx,'height4pi',wy,'mirror4pi',[0 0 0 0],'pixelsize_nm',100,'offset',100,'conversion',0.5);
    end
end