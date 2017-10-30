function zs=testfit(teststack,coeff,shiftxy, p,linepar,ax)
if nargin<4
    linepar={};
elseif ~iscell(linepar)
    linepar={linepar};
end
d=round((size(teststack,1)-p.ROIxy)/2);
            range=d+1:d+p.ROIxy;

numstack=size(teststack,4);
t=tic;
% f=figure(989);ax2=gca;hold off
    for k=1:size(teststack,4)
        if toc(t)>1
            p.status.String=['fitting test stacks: ' num2str(k/numstack,'%1.2f')];drawnow
            t=tic;
        end
        if contains(p.modality,'2D')
            fitmode=6;
        else
            fitmode=5;
        end
        fitstack=single(squeeze(teststack(range,range,:,k,:)));
        coeffh(:,:,:,:,1)=single(coeff{1});
        coeffh(:,:,:,:,2)=single(coeff{2});
        shared=[1 1 1 0 0]';
        nfits=size(fitstack,3);
        npar=5;
        dT=zeros(npar,2,nfits);
        dT(1,2,:)=shiftxy(k,2);
        dT(2,2,:)=shiftxy(k,1);
        iterations=50;
%         [PM,CRLBM, LLM,update, error] =  kernel_MLEfit_Spline_LM_multichannel_finalized(fitstack,coeffh, shared,dT,50);
        sharedA = repmat(shared,[1 size(fitstack,3)]);
        [P,CRLB, LL] =CPUmleFit_LM_MultiChannel(fitstack,int32(sharedA),iterations,coeffh,single(dT));
        %compare with sinlge fits
        [P1,CRLB1, LL1] =mleFit_LM(fitstack(:,:,:,1),5,iterations,single(coeff{1}),0,1);
        [P2,CRLB2, LL2] =mleFit_LM(fitstack(:,:,:,2),5,iterations,single(coeff{2}),0,1);
        %compare with all free fits
        [Pf,CRLBr, LLf] =CPUmleFit_LM_MultiChannel(fitstack,int32(sharedA*0),iterations,coeffh,single(dT));
        
        
        %define one as reference and plot differences 
        figure(78);
        rs=02;
%         subplot(2,2,1)
hold off
        plot(P1(:,1),P2(:,1)+0*shiftxy(k,1),'+')
        xlim([-1 1]*rs+(p.ROIxy-1)/2);
        ylim([-1 1]*rs+(p.ROIxy-1)/2);
        xlabel('x1 single fit')
        ylabel('x')
%         subplot(2,2,2)
hold on
        plot(Pf(:,1),Pf(:,2),'o')
        xlim([-1 1]*rs+(p.ROIxy-1)/2);
        ylim([-1 1]*rs+(p.ROIxy-1)/2);
        xlabel('x1 not linked')
        ylabel('x')
%         subplot(2,2,3)

        plot(P1(:,1),P(:,1),'x')
        
        plot([-1 1]*rs+(p.ROIxy-1)/2,[-1 1]*rs+(p.ROIxy-1)/2)
        xlim([-1 1]*rs+(p.ROIxy-1)/2);
        ylim([-1 1]*rs+(p.ROIxy-1)/2);
        
        
        title(shiftxy(k,1))
        xlabel('x1 single fit')
        ylabel('x ')
        legend('x2 individual fit','x2 not linked','x linked')
        
%         [P] =  mleFit_LM(single(squeeze(teststack(range,range,:,k))),fitmode,100,single(coeff),0,1);
        
        z=(1:size(P,1))'-1;

        znm=(P(:,3)-p.z0)*p.dz;
        plot(ax,z,znm,linepar{:})
        hold(ax,'on')
        xlabel(ax,'frame')
        ylabel(ax,'zfit (nm)')
        zs(:,k)=P(:,5);
        
        
        if 1% imageslicer to test
%             coord=P1(:,[1 2 5 3 4]);
            coord=P(:,[1 2 3 4 6]);
            coord2=coord;
            coord2(:,1)=coord2(:,1)+squeeze(dT(1,2,:));
            coord2(:,2)=coord2(:,2)+squeeze(dT(2,2,:));
            img1=renderPSF(coeff{1},coord,size(fitstack,1));
            img2=renderPSF(coeff{2},coord2,size(fitstack,1));
            imall=[fitstack(:,:,:,1),img1; fitstack(:,:,:,2),img2];
            res=[fitstack(:,:,:,1)-img1, 0*img1;fitstack(:,:,:,2)-img2,0*img2];
            ims(:,:,:,1)=imall;ims(:,:,:,2)=res;
            imageslicer(ims)
        end
        
% test for the returned photons and photons in the raw image        
%         phot=P(:,3); bg=P(:,4);
%         totsum=squeeze(nansum( nansum(teststack(range,range,:,k),1),2));
%         totsum=totsum-squeeze(min(min(teststack(range,range,:,k),[],1),[],2))*length(range)^2;
%         photsum=phot+0*bg*length(range)^2;
%         plot(ax2,z,(photsum-totsum)./totsum,'.')
%         hold(ax2,'on')
    end
    
end
