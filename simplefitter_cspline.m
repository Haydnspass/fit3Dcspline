function simplefitter_cspline(p)
%p.filename
fittime=0;
fitsperblock=50000;
imstack=zeros(p.roifit,p.roifit,fitsperblock,'single');
peakcoordinates=zeros(fitsperblock,3);
indstack=0;
resultsind=1;
%results
% frame, x,y,z,phot,bg, errx,erry, errz,errphot, errbg,logLikelihood
%load calibration
if exist(p.calfile,'file')
    cal=load(p.calfile);
else
    errordlg('please select 3D calibration file')
end
p.dz=cal.cspline.dz;  %coordinate system of spline PSF is corner based and in units pixels / planes
p.z0=cal.cspline.z0;
p.coeff=cal.cspline.coeff;
p.dx=floor(p.roifit/2);
% readerome=bfGetReader(p.imagefile);
 p.status.String=['Open tiff file' ]; drawnow

 switch p.loader
     case 1
        reader=mytiffreader(p.imagefile);
        numframes=reader.info.numberOfFrames;
     case 2
        reader=bfGetReader(p.imagefile);
        numframes=reader.getImageCount;
     case 3 %fiji
          ij=p.mij.imagej;
          ijframes=ij.getFrames;
          for k=1:length(ijframes)
            if strcmp(ijframes(k).class,'ij.gui.StackWindow')&&~isempty(ijframes(k).getImagePlus)    
                reader=ijframes(k).getImagePlus.getStack;
                break
            end
          end
          numframes=reader.size;
 end


if p.preview
    frames=min(p.previewframe,numframes);
else
    frames=1:numframes;
end

%loop over frames, do filtering/peakfinding
h=fspecial('gaussian',ceil(3*p.peakfilter+1),p.peakfilter);
tshow=tic;
for F=frames
    image=getimage(F,reader,p);
%     image=reader.read(F);
%      image2=bfGetPlane(readerome,F);
    sim=size(image);
    imphot=(single(image)-p.offset)*p.conversion;
    bg=mywaveletfilter(imphot,3,false,true);
    impf=filter2(h,imphot-bg);
    maxima=maximumfindcall(impf);
    maxgood=maxima(maxima(:,3)>p.peakcutoff,:);
    
    if p.preview && size(maxgood,1)>2000
        p.status.String=('increase cutoff');
        break
    end
    
    %cut out images
    for k=1:size(maxgood,1)
        if maxgood(k,1)>p.dx && maxgood(k,2)>p.dx && maxgood(k,1)<= sim(2)-p.dx && maxgood(k,2)<=sim(1)-p.dx 
            indstack=indstack+1;
            if p.mirror
                imstack(:,:,indstack)=imphot(maxgood(k,2)-p.dx:maxgood(k,2)+p.dx,maxgood(k,1)+p.dx:-1:maxgood(k,1)-p.dx);
            else
                imstack(:,:,indstack)=imphot(maxgood(k,2)-p.dx:maxgood(k,2)+p.dx,maxgood(k,1)-p.dx:maxgood(k,1)+p.dx);
            end
            peakcoordinates(indstack,1:2)=maxgood(k,1:2);
            peakcoordinates(indstack,3)=F;
   
            if indstack==fitsperblock
                p.status.String=['Fitting...' ]; drawnow
                t=tic;
                resultsh=fitspline(imstack,peakcoordinates,p);
                fittime=fittime+toc(t);
                
                results(resultsind:resultsind+fitsperblock-1,:)=resultsh;
                resultsind=resultsind+fitsperblock;
                
                indstack=0;
            end
        end
       
    end
    if toc(tshow)>1
        tshow=tic;
        p.status.String=['Loading frame ' num2str(F) ' of ' num2str(numframes)]; drawnow
    end
              
end
closereader(reader,p);
p.status.String=['Fitting last stack...' ]; drawnow
if indstack<1
    p.status.String=['No localizations found. Increase cutoff?' ]; drawnow
else

t=tic;
resultsh=fitspline(imstack(:,:,1:indstack),peakcoordinates(1:indstack,:),p); %fit all the rest
fittime=fittime+toc(t);

results(resultsind:resultsind+indstack-1,:)=resultsh;
end
if p.preview
    figure(201)
    imagesc(impf);
     colorbar
    hold on
    plot(maxgood(:,1),maxgood(:,2),'wo')
    plot(results(:,2),results(:,3),'k+')
    hold off
   
    p.status.String=['Preview done. ' num2str(size(results,1)/fittime,'%3.0f') ' fits/s. ' num2str(size(results,1),'%3.0f') ' localizations.']; drawnow
else
    p.status.String=['Fitting done. ' num2str(size(results,1)/fittime,'%3.0f') ' fits/s. ' num2str(size(results,1),'%3.0f') ' localizations. Saving now.']; drawnow
    results(:,[13,15])=results(:,[2, 7])*p.pixelsize(1);
    results(:,[14,16])=results(:,[3, 8])*p.pixelsize(end);
    
    resultstable=array2table(results,'VariableNames',{'frame','x_pix','y_pix','z_nm','photons','background',' crlb_x','crlb_y','crlb_z','crlb_photons','crlb_background','logLikelyhood','x_nm','y_nm','crlb_xnm','crlb_ynm'});
    % 
    writenames=true;
    if contains(p.outputformat,'pointcloud')
        resultstable=resultstable(:,[1 13 14 4 15 9]); %for pointcloud-loader
        del='\t';
        disp('Load in http://www.cake23.de/pointcloud-loader/')
    elseif contains(p.outputformat,'ViSP')
        writenames=false;
        del='\t';
         resultstable=resultstable(:,[13 14 4 15 16 9 5 1]);
         [path,file]=fileparts(p.outputfile);
         p.outputfile=fullfile(path, [file '.3dlp']);
         disp('Load in Visp: https://science.institut-curie.org/research/multiscale-physics-biology-chemistry/umr168-physical-chemistry/team-dahan/softwares/visp-software-2/')
    else
        del=',';
    end

    writetable(resultstable,p.outputfile,'Delimiter',del,'FileType','text','WriteVariableNames',writenames);
    p.status.String=['Fitting done. ' num2str(size(results,1)/fittime,'%3.0f') ' fits/s. ' num2str(size(results,1),'%3.0f') ' localizations. Saved.']; drawnow
end
end

function results=fitspline(imstack,peakcoordinates,p)
if p.bidirectional
    fitmode=6;
else
    fitmode=5;
end

[Pcspline,CRLB,LL]=mleFit_LM(imstack,fitmode,50,single(p.coeff),0,1);
results=zeros(size(imstack,3),12);
results(:,1)=peakcoordinates(:,3);
if  p.mirror
    results(:,2)=p.dx-Pcspline(:,2)+peakcoordinates(:,1);
else
    results(:,2)=Pcspline(:,2)-p.dx+peakcoordinates(:,1);
    
end
% frame, x,y,z,phot,bg, errx,erry, errz,errphot, errbg,logLikelihood
results(:,3)=Pcspline(:,1)-p.dx+peakcoordinates(:,2); %x,y in pixels 
results(:,4)=(Pcspline(:,5)-p.z0)*p.dz;
results(:,5:6)=Pcspline(:,3:4);
results(:,7:8)=sqrt(CRLB(:,[2 1]));
results(:,9)=sqrt(CRLB(:,5)*p.dz);
results(:,10:11)=sqrt(CRLB(:,3:4));
results(:,12)=LL;

end

function img=getimage(F,reader,p)
switch p.loader
     case 1
         img=reader.read(F);
    case 2
        img=bfGetPlane(reader,F);
    case 3        
        ss=[reader.getWidth reader.getHeight reader.getSize];
        if F>0&&F<=ss(3)
            pixel=reader.getPixels(F);
            img=reshape(pixel,ss(1),ss(2))';
        else
            img=[];
        end
end

end

function closereader(reader,p)
switch p.loader
     case 1
         reader.close;
    case 2
    case 3        

end

end