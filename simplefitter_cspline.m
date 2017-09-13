function simplefitter_cspline(p)

%parameters:
% p.imagefile: fiename of data (char);
% p.calfile: filename of calibration data (char);
% p.offset=ADU offset of data;
% p.conversion=conversion e-/ADU;
% p.preview: true if preview mode (fit only current image and display
% results).
% p.previewframe=frame to preview;
% p.peakfilter=filtersize (sigma, Gaussian filter) for peak finding;
% p.peakcutoff=cutoff for peak finding
% p.roifit=size of the ROI in pixels
% p.bidirectional= use bi-directional fitting for 2D data
% p.mirror=mirror images if bead calibration was taken without EM gain
% p.status=handle to a GUI object to display the status;
% p.outputfile=file to write the localization table to;
% p.outputformat=Format of file;
% p.pixelsize=pixel size in nm;


fittime=0;
fitsperblock=50000;
imstack=zeros(p.roifit,p.roifit,fitsperblock,'single');
peakcoordinates=zeros(fitsperblock,3);
indstack=0;
resultsind=1;
%results
% frame, x,y,z,phot,bg, errx,erry, errz,errphot, errbg,logLikelihood
%load calibration
cal=load(p.calfile);
p.dz=cal.cspline.dz;  %coordinate system of spline PSF is corner based and in units pixels / planes
p.z0=cal.cspline.z0;
p.coeff=cal.cspline.coeff;
p.dx=floor(p.roifit/2);
% readerome=bfGetReader(p.imagefile);
reader=mytiffreader(p.imagefile);
numframes=reader.info.numberOfFrames;
if p.preview
    frames=min(p.previewframe,numframes);
else
    frames=1:numframes;
end

%loop over frames, do filtering/peakfinding
h=fspecial('gaussian',ceil(3*p.peakfilter+1),p.peakfilter);
tshow=tic;
for F=frames
    image=reader.read(F);
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
p.status.String=['Fitting...' ]; drawnow
t=tic;
resultsh=fitspline(imstack(:,:,1:indstack),peakcoordinates(1:indstack,:),p); %fit all the rest
fittime=fittime+toc(t);
results(resultsind:resultsind+indstack-1,:)=resultsh;

if p.preview
    figure(201)
    imagesc(impf);
    hold on
    plot(maxgood(:,1),maxgood(:,2),'wo')
    plot(results(:,2),results(:,3),'k+')
    hold off
    colorbar
    p.status.String=['Preview done. ' num2str(size(results,1)/fittime,'%3.0f') ' fits/s. ' num2str(size(results,1),'%3.0f') ' localizations. Saved.']; drawnow
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
        del=',\t';
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