function calibrate_4pi(p)
%for now now segmentation of calibration. 
ph=p;
f=figure('Name','Bead calibration');
tg=uitabgroup(f);
t1=uitab(tg,'Title','prefit');
tgprefit=uitabgroup(t1);
ph.tabgroup=   tgprefit;
ph.isglobal=false;
ph.outputfile={};
    
    
[beads,ph]=images2beads_globalfit(ph);


s=size(beads(1).stack.image);
range=(s(1)+1)/2+[-1 1];

ybeads=getFieldAsVectorInd(beads,'pos',2);
y4pi=p.smappos.settings.y4pi;
x4pi=p.smappos.settings.x4pi;
height4pi=p.smappos.settings.height4pi; 
width4pi=p.smappos.settings.width4pi; 
mirror4pi=p.smappos.settings.mirror4pi;

for k=1:length(y4pi)
    indbh=ybeads>=y4pi(k)&ybeads<=y4pi(k)+height4pi;
    
    xposh=getFieldAsVectorInd(beads(indbh),'pos',1)';
    yposh=getFieldAsVectorInd(beads(indbh),'pos',2)';
    ph.tabgroup=uitabgroup(uitab(tgprefit,'Title',num2str(k)));
    [splinefit,indgood,posbeads,shift]=getstackcal_g(beads(indbh),ph);
    
    beadtrue{k}(:,1)=xposh(indgood)+shift(indgood,1); %could also be -shift. Test later!
    beadtrue{k}(:,2)=yposh(indgood)+shift(indgood,2);
    beadtrue{k}(:,3)=shift(indgood,3);
    
end
%calculate transformN
%transform_simple for transformN: 
figure(88);hold off
transform=interfaces.LocTransformN;
pt.mirror=mirror4pi(1)*2;
pt.xrange=[x4pi(1) x4pi(1)+width4pi];
pt.yrange=[y4pi(1) y4pi(1)+height4pi];
pt.units='pixel';
transform.setTransform(1,pt)
for k=2:length(beadtrue)
    plot(beadtrue{k}(:,1),beadtrue{k}(:,2),'.');hold on
    pt.mirror=(k)*2;
    pt.xrange=[x4pi(k) x4pi(k)+width4pi];
    pt.yrange=[y4pi(k) y4pi(k)+height4pi];
    transform.setTransform(k,pt)
    transform_locs_simpleN(transform, beadtrue{1},beadtrue{k},p)
end


figure(99);hold off;
for k=1:length(beads)
    intz=squeeze(sum(sum(beads(k).stack.image(range,range,:),1),2));
    if beads(k).pos(2)>p.smappos.settings.y4pi(1) && beads(k).pos(2)<p.smappos.settings.y4pi(1)+p.smappos.settings.height4pi 

    plot(intz)  
    hold on
    end
end


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