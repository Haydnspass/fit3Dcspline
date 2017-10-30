function img=simulate2c
% transformfile='/Users/jonas/Documents/Data/ROI2_20per639_50msexp_1/ROI2_20per639_50msexp_1_MMStack2_T.mat';
cal3dfile='/Users/jonas/Documents/Data/ROI2_20per639_50msexp_1_3dcal.mat'; 
transformfile='/Users/jonas/Documents/Data/simulstack_affine_2_T.mat';

nbeads=25;
roisize=21;
dr=(roisize-1)/2;
zpos=1:200; %in slices
pixelsize=100;
Nphot=10000;
Nbg=100;
img=zeros(512,512,length(zpos),'single');

x=rand(nbeads,1)*(size(img,1)-roisize*2)+roisize;
y=rand(nbeads,1)*(size(img,2)-roisize*2)/2+roisize;

xnm=x*pixelsize; ynm=y*pixelsize;

t=load(transformfile); transformation=t.transformation;
t=load(cal3dfile); coeff=t.SXY.cspline.coeff;



[xtnm,ytnm]=transformation.transformCoordinatesFwd(xnm,ynm);
xt=xtnm/pixelsize;yt=ytnm/pixelsize;
figure(88);plot(x,y,'+',xt,yt,'x')
coord=horzcat(y,x);

channel=1;

makeimg;
coord=horzcat(yt,xt);
channel=2;
makeimg
img=img+Nbg;
imageslicer(img);

imgnoise=poissrnd(img);
outfile=[fileparts(cal3dfile) filesep 'simulstack.tif'];
saveastiff(uint16(imgnoise),outfile);

function makeimg
roipos=round(coord);
incoord=coord-roipos+dr;
incoord(:,4)=Nphot;
incoord(:,5)=0;
for zk=1:length(zpos)
    incoord(:,3)=zpos(zk);
    impsf=renderPSF(coeff{channel},incoord,roisize);
    for k=1:nbeads
        if all(roipos(k,:)>dr)&&all(roipos(k,:)<size(img,1)-dr+1)
       img(roipos(k,2)-dr:roipos(k,2)+dr,roipos(k,1)-dr:roipos(k,1)+dr,zk)=img(roipos(k,2)-dr:roipos(k,2)+dr,roipos(k,1)-dr:roipos(k,1)+dr,zk)+impsf(:,:,k);
    
        end
    end
   
end
end
end