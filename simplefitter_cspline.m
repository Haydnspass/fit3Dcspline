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
cal=load(p.calfile);
p.dz=cal.cspline.dz;  %coordinate system of spline PSF is corner based and in units pixels / planes
p.z0=cal.cspline.z0;
p.coeff=cal.cspline.coeff;
p.dx=floor(p.roifit/2);

reader=bfGetReader(p.imagefile);
numframes=reader.getImageCount();
if p.preview
    frames=min(p.previewframe,numframes);
else
    frames=1:numframes;
end

%loop over frames, do filtering/peakfinding
h=fspecial('gaussian',ceil(3*p.peakfilter+1),p.peakfilter);
tshow=tic;
for F=frames
    image=bfGetPlane(reader,F);
    sim=size(image);
    imphot=(single(image)-p.offset)*p.conversion;
    bg=mywaveletfilter(imphot,3,false,true);
    impf=filter2(h,imphot-bg);
    maxima=maximumfindcall(impf);
    maxgood=maxima(maxima(:,3)>p.peakcutoff,:);
    
    %cut out images
    for k=1:size(maxgood,1)
        if maxgood(k,1)>p.dx && maxgood(k,2)>p.dx && maxgood(k,1)<= sim(2)-p.dx && maxgood(k,2)<=sim(1)-p.dx 
            indstack=indstack+1;
            if p.mirror
                imstack(:,:,indstack)=image(maxgood(k,2)-p.dx:maxgood(k,2)+p.dx,maxgood(k,1)+p.dx:-1:maxgood(k,1)-p.dx);
            else
                imstack(:,:,indstack)=image(maxgood(k,2)-p.dx:maxgood(k,2)+p.dx,maxgood(k,1)-p.dx:maxgood(k,1)+p.dx);
            end
            peakcoordinates(indstack,1:2)=maxgood(k,1:2);
            peakcoordinates(indstack,3)=F;
   
            if indstack==fitsperblock
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
        p.status.String=['Fitting frame ' num2str(F) ' of ' num2str(numframes)]; drawnow
    end
              
end

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
end
p.status.String=['Fitting done. ' num2str(size(results,1)/fittime,'%3.0f') ' fits/s. ' num2str(size(results,1),'%3.0f') ' localizations.'];

resultstable=array2table(results,'VariableNames',{'frame','x_pix','y_pix','z_nm','photons','background',' crlb_x','crlb_y','crlb_z','crlb_photons','crlb_background','logLikelyhood'});

writetable(resultstable,p.outputfile);
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
    results(:,2)=p.dx-Pcspline(:,1)+peakcoordinates(:,1);
else
    results(:,2)=Pcspline(:,1)-p.dx+peakcoordinates(:,1);
    
end
% frame, x,y,z,phot,bg, errx,erry, errz,errphot, errbg,logLikelihood
results(:,3)=Pcspline(:,2)-p.dx+peakcoordinates(:,2); %x,y in pixels 
results(:,4)=(Pcspline(:,5)-p.z0)*p.dz;
results(:,5:6)=Pcspline(:,3:4);
results(:,7:8)=CRLB(:,1:2);
results(:,9)=CRLB(:,5)*p.dz;
results(:,10:11)=CRLB(:,3:4);
results(:,12)=LL;

end