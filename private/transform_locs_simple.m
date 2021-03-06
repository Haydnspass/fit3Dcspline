function transform=transform_locs_simple(locrefi,loctargeti,p)
df=10;
locref=reducepos(locrefi,df);
loctarget=reducepos(loctargeti,df);


tfile=p.Tfile;

if exist(tfile,'file')
    Tinitial=load(tfile);
    [loctT.x,loctT.y]=Tinitial.transformCoordinatesInv(loctarget.x(:),loctarget.y(:));
    mirrorinfo=Tinitial.tinfo.mirror;
    dx=0;dy=0;
%     if contains(mirrorinfo.targetmirror,'no')
%     cutout=true;
%     end
    %     pos=Tinitial.pos;
%     size=Tinitial.size;
else %all initial estimation:
%     approximate shift from size and position
    separator=256;
    sepscale=1; %maximum separation measure
    loctT=loctarget;
    separators=[512 512];
    switch p.Tmode
        case 'up-down'
            loctT.y=loctarget.y(:)-separator;
            targetmirror='no mirror';
            targetpos='bottom';
            separators(2)=separator;
        case 'up-down mirror'
            loctT.y=-loctarget.y(:)+2*separator;
            targetpos='bottom';
            targetmirror='up-down';
            separators(2)=separator;
        case 'right-left'
            loctT.x=loctarget.x(:)-separator;
            targetmirror='no mirror';
            targetpos='right';
            separators(1)=separator;
        case 'right-left mirror'
            loctT.x=-loctarget.x(:)+2*separator;   
            targetmirror='left-right';
            targetpos='right';
            separators(1)=separator;
    end
    mirrorinfo.targetmirror=targetmirror;
    %determine approximate shift
    xr=1:1:2*separator;yr=xr;
    ht=histcounts2(locref.x,locref.y,xr,yr);
    hr=histcounts2(loctT.x,loctT.y,xr,yr);
    G=fftshift(ifft2(conj(fft2(ht)).*fft2(hr)));
    h=fspecial('gaussian',13,sepscale);
    Gf=filter2(h,G);
    [~ ,indmax]=max(Gf(:));
    [x0,y0]=ind2sub(size(Gf),indmax);
    dx0=x0-ceil(size(Gf,1)/2);
    dy0=y0-ceil(size(Gf,2)/2);
%     loctT.x=loctT.x-dx;
%     loctT.y=loctT.y-dy;
  

end


[iAa,iBa,na,nb,nseen]=matchlocsall(locref,loctT,-dx0,-dy0,2*sepscale,1e5);

transform=interfaces.LocTransform;
t.type='polynomial';
% t.type='affine';
t.parameter=3;
transform.findTransform(locref.x(iAa),locref.y(iAa),loctarget.x(iBa),loctarget.y(iBa),t)
transform.findTransformZ(locref.x(iAa),locref.y(iAa),locref.z(iAa),loctarget.x(iBa),loctarget.y(iBa),loctarget.z(iBa),t)

 [xa, ya, za]=transform.transformCoordinatesInv((loctarget.x(iBa)),(loctarget.y(iBa)),(loctarget.z(iBa)));
    dx=xa-locref.x(iAa);
   dy=ya-locref.y(iAa);
   
  figure(88);plot(locref.x,locref.y,'b.',loctT.x,loctT.y,'r+',loctT.x-dx0,loctT.y-dy0,'g.',loctargeti.x,loctargeti.y,'rx',xa,ya,'cx') 
   
if isfield(p,'tabgroup')
    axh=axes(p.tabgroup);
else
    figure(99);
    axh=gca;
end
axes(axh);
dscatter(dx,dy)
title([num2str(std(dx)) ', ' num2str(std(dy))]);

transform.tinfo.targetpos=targetpos;
transform.tinfo.separator=separators;
transform.tinfo.mirror=mirrorinfo;
transform.tinfo.cam_pixelsize_nm=1;
transform.tinfo.units='pixels';

end


function pos=reducepos(posin,df)
    z0=ceil(size(posin.x,1)/2);
    for l=size(posin.x,2):-1:1
        framerange=abs(posin.frame(:,l)-z0)<=df;
        x(:,l)=posin.x(framerange,l);y(:,l)=posin.y(framerange,l);z(:,l)=posin.z(framerange,l);frame(:,l)=posin.frame(framerange,l);
    end
    [~,indsort]=sort(frame(:));
    pos.x=x(indsort);pos.y=y(indsort);pos.z=z(indsort);pos.frame=frame(indsort);
    
end