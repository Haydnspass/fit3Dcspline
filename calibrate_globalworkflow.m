function calibrate_globalworkflow(p)

if ~p.isglobalfit %only normal local calibratin, just call old proram
    calibrate3D_g(p);
    return
end

global S beadpos1 beadpos2
ph=p;
ph.isglobalfit=false;
%set spatial calibration
ph.outputfile={};
splitpos=p.Tsplitpos;
% splitpos=256;% later to GUI?
% splitpos=150 %challenge

if ~isfield(p,'yrange')
    p.yrange=[-inf inf];
end
if  ~isfield(p,'xrange')
    p.xrange=[-inf inf];
end

switch p.Tmode
    case {'up-down','up-down mirror'}
        if max(p.yrange)<splitpos %defined only in upper part
            yrange1=p.yrange;
            if contains(p.Tmode,'mirror')
                yrange2=sort(-p.yrange+2*splitpos);yrange2(yrange2<splitpos)=splitpos;
            else
                yrange2=sort(p.yrange+splitpos);yrange2(yrange2<splitpos)=splitpos;
            end
        else
            yrange1=([p.yrange splitpos]);yrange1(yrange1>splitpos)=splitpos;yrange1=unique(yrange1);
            yrange2=([p.yrange+ splitpos]);yrange2(yrange2<splitpos)=splitpos;yrange2=unique(yrange2);
        end
            
%         yrange=unique([p.yrange splitpos p.yrange+splitpos]);
%         yrange1=yrange(yrange<=splitpos);yrange2=yrange(yrange>=splitpos);
        xrange1=p.xrange;xrange2=p.xrange;
        yrangeall=[yrange1(1:end-1) splitpos yrange2(2:end)];
        yrange1(end)=yrange1(end)-p.mindistance; %do not take into account locs too close to separator
        yrange2(1)=yrange2(1)+p.mindistance;
        xrangeall=p.xrange;
        XYpos=[1,2];
        
        split='ud';
    case {'right-left','right-left mirror'}
        if max(p.xrange)<splitpos %defined only in upper part
            xrange1=p.xrange;
            xrange2=p.xrange+splitpos;xrange2(xrange2<splitpos)=splitpos;
        else
            xrange1=([p.xrange splitpos]);xrange1(xrange1>splitpos)=splitpos;xrange1=unique(xrange1);
            xrange2=([p.xrange+ splitpos]);xrange2(xrange2<splitpos)=splitpos;xrange2=unique(xrange2);
        end
        yrange1=p.yrange;yrange2=p.yrange;
        yrangeall=p.yrange;
        xrangeall=[xrange1(1:end-1) splitpos xrange2(2:end)];
        XYpos=[2,1];
        split='rl';
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



S1.Yrangeall=yrangeall;S1.Xrangeall=xrangeall;
S2.Yrangeall=yrangeall;S2.Xrangeall=xrangeall;
S2.posind=XYpos;
% [S,beadpos]=calibrate3D_g(ph);

% Later: also do test-fitting with corresponding spline coefficients
p.tabgroup=uitab(tg,'Title','transformation');
p.separator=splitpos;
% find transform
if p.makeT || isempty(p.Tfile)
    transform=transform_locs_simple(beadpos1{1},beadpos2{1},p);
else
    l=load(p.Tfile);
    transform=l.transformation;
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

if strcmp(split,'rl')
    SXY(1:length(S1),1)=S1;
    SXY(end+1:end+length(S2),1)=S2;
else
    SXY(1,1:length(S1))=S1;
    SXY(1,end+1:end+length(S2))=S2;
end
SXY_g=S;
transformation=parameters_g.transformation;
calibrationfigure=f;
if ~isempty(p.outputfile)
    if p.smap
        parameters1.smappos.P=[]; parameters2.smappos.P=[]; parameters_g.smappos.P=[];
        save(p.outputfile,'SXY','SXY_g','parameters_g','parameters1','parameters2','transformation');
        
    else
        save(p.outputfile,'gausscal','cspline_all','gauss_sx2_sy2','gauss_zfit','cspline','parameters');
    end
    filefig=strrep(p.outputfile,'.mat','.fig');
    savefig(calibrationfigure,filefig,'compact');
end