%  Copyright (c)2017 Ries Lab, European Molecular Biology Laboratory,
%  Heidelberg.
%  
%  This file is part of GPUmleFit_LM Fitter.
%  
%  GPUmleFit_LM Fitter is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%  
%  GPUmleFit_LM Fitter is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%  
%  You should have received a copy of the GNU General Public License
%  along with GPUmleFit_LM Fitter.  If not, see <http://www.gnu.org/licenses/>.
%  
%  
%  Additional permission under GNU GPL version 3 section 7
%  
%  If you modify this Program, or any covered work, by
%  linking or combining it with libraries required for interaction
%  with analysis programs such as Igor Pro or Matlab,
%  the licensors of this Program grant you additional permission
%  to convey the resulting work.

%%
classdef simplefitter_GUI<handle
    properties
        guihandles
        smappos
        fijipath
        mij
    end
    methods
        function obj=simplefitter_GUI(varargin)  
            %constructur: make GUI   
            addpath('shared')
            addpath('bfmatlab')
            
            figureheight=630;
            h=figure('Name','3D fittter','MenuBar','none','ToolBar','none');
            initPosition = h.Position;
            h.Position=[initPosition(1), initPosition(2)- figureheight+initPosition(4),450, figureheight];
            top=h.Position(4)-10;
            vsep=24;
            if ispc
                fontsize=12;
            else 
                fontsize=14;
            end
            xpos1=10;
            xw=100;
            hatitle='left';
            obj.guihandles.handle=h;
            obj.guihandles.title=uicontrol('style','text','String','Fit image stack with experimental PSF model. (c) 2017 Ries lab','Position',[xpos1,top-vsep+10,xw*4.5,vsep],'FontSize',10,'HorizontalAlignment',hatitle,'FontWeight','bold');
            
            obj.guihandles.loadert=uicontrol('style','text','String','Loader','Position',[xpos1,top-2*vsep,xw*1.5,vsep],'FontSize',fontsize);
            obj.guihandles.loader=uicontrol('style','popupmenu','String',{'simple tif','ome loader','ImageJ'},'Position',[xpos1+1.5*xw,top-2*vsep,xw*2.5,vsep],'FontSize',fontsize,'Callback',@obj.changeloader_callback,'Value',3);
                      
            obj.guihandles.selectfiles=uicontrol('style','pushbutton','String','Load raw images','Position',[xpos1,top-3*vsep,xw*1.5,vsep],'FontSize',fontsize,'Callback',{@obj.selectfiles_callback,1});
%             obj.guihandles.selectfiles.TooltipString='Select image files with bead stacks. You can select several files from different locations with the file select dialog box opend';
            obj.guihandles.imagefile=uicontrol('style','edit','String','','Position',[xpos1+1.5*xw,top-3*vsep,xw*2.5,vsep],'FontSize',fontsize);
%             obj.guihandles.imagefile.TooltipString='List of image files used for calibration. To change this list, use select camera files';
            
            obj.guihandles.selectcoeff=uicontrol('style','pushbutton','String','Load calibration','Position',[xpos1,top-4*vsep,xw*1.5,vsep],'FontSize',fontsize,'Callback',{@obj.selectfiles_callback,2});
%             obj.guihandles.selectcoeff.TooltipString='Select image files with bead stacks. You can select several files from different locations with the file select dialog box opend';
            obj.guihandles.calfile=uicontrol('style','edit','String','','Position',[xpos1+1.5*xw,top-4*vsep,xw*2.5,vsep],'FontSize',fontsize);
%             obj.guihandles.calfile.TooltipString='List of image files used for calibration. To change this list, use select camera files';
            
            obj.guihandles.camtext=uicontrol('style','text','String','Acquisition parameters:','Position',[xpos1,top-6*vsep,xw*4,vsep],'FontSize',fontsize,'HorizontalAlignment',hatitle,'FontWeight','bold');
            
            ha='right';
            
            obj.guihandles.conversiont=uicontrol('style','text','String','conversion (e-/ADU)','Position',[xpos1,top-7*vsep,xw*1.5,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.conversion=uicontrol('style','edit','String','0.1','Position',[xpos1+xw*1.5,top-7*vsep,xw*.5,vsep],'FontSize',fontsize);
            obj.guihandles.conversion.TooltipString=sprintf('conversion factor = conv/EMgain. conv is the gain stated in the camera spec sheet (e-/ADU)');
%             obj.guihandles.dzt.TooltipString=obj.guihandles.dz.TooltipString;
            obj.guihandles.offset=uicontrol('style','text','String','offset (ADU)','Position',[xpos1+2*xw,top-7*vsep,xw*1.5,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.offset=uicontrol('style','edit','String','100','Position',[xpos1+xw*3.5,top-7*vsep,xw*.5,vsep],'FontSize',fontsize);
            
            obj.guihandles.peaktext=uicontrol('style','text','String','Peak candidate finding:','Position',[xpos1,top-9*vsep,xw*4,vsep],'FontSize',fontsize,'HorizontalAlignment',hatitle,'FontWeight','bold');

            obj.guihandles.peakfiltert=uicontrol('style','text','String','Filter size (pixel)','Position',[xpos1,top-10*vsep,xw*1.5,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.peakfilter=uicontrol('style','edit','String','1.2','Position',[xpos1+xw*1.5,top-10*vsep,xw*.5,vsep],'FontSize',fontsize);
%             obj.guihandles.dz.TooltipString=sprintf('Distance in nm between frames. By convention, these are objective positions (not corrected for refractive index mismatch). \n A spacing between 10 nm and 50 nm works well ');
%             obj.guihandles.dzt.TooltipString=obj.guihandles.dz.TooltipString;
            obj.guihandles.peakcutofft=uicontrol('style','text','String','cutoff (photons)','Position',[xpos1+2*xw,top-10*vsep,xw*1.5,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.peakcutoff=uicontrol('style','edit','String','1','Position',[xpos1+xw*3.5,top-10*vsep,xw*.5,vsep],'FontSize',fontsize);
            
             obj.guihandles.fittxt=uicontrol('style','text','String','Fitting paramters:','Position',[xpos1,top-12*vsep,xw*4,vsep],'FontSize',fontsize,'HorizontalAlignment',hatitle,'FontWeight','bold');

            obj.guihandles.roifitt=uicontrol('style','text','String','ROI size (pixel)','Position',[xpos1,top-13*vsep,xw*1.5,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.roifit=uicontrol('style','edit','String','13','Position',[xpos1+xw*1.5,top-13*vsep,xw*.5,vsep],'FontSize',fontsize);
%             obj.guihandles.dz.TooltipString=sprintf('Distance in nm between frames. By convention, these are objective positions (not corrected for refractive index mismatch). \n A spacing between 10 nm and 50 nm works well ');
%             obj.guihandles.dzt.TooltipString=obj.guihandles.dz.TooltipString;
            obj.guihandles.bidirectional=uicontrol('style','checkbox','String','2D','Position',[xpos1+2*xw,top-13*vsep,xw*1,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.mirror=uicontrol('style','checkbox','String','mirror','Position',[xpos1+xw*3,top-13*vsep,xw*1,vsep],'FontSize',fontsize);
            
            
            
            obj.guihandles.selectoutput=uicontrol('style','pushbutton','String','set output','Position',[xpos1,top-15*vsep,xw*1.5,vsep],'FontSize',fontsize,'Callback',@obj.selectoutput_callback);
            obj.guihandles.selectoutput.TooltipString='Select image files with bead stacks. You can select several files from different locations with the file select dialog box opend';
            obj.guihandles.outputfile=uicontrol('style','edit','String','','Position',[xpos1,top-16*vsep,xw*4,vsep],'FontSize',fontsize);
            obj.guihandles.outputfile.TooltipString='List of image files used for calibration. To change this list, use select camera files';
            
            obj.guihandles.outputformat=uicontrol('style','popupmenu','String',{'csv','pointcloud-loader','ViSP'},'Position',[xpos1+1.5*xw,top-15*vsep,xw*2.5,vsep],'FontSize',fontsize);
            obj.guihandles.outputformat.TooltipString='Choose output format. CSV saves all fit parameters';
            
            
            obj.guihandles.pixelsizet=uicontrol('style','text','String','pixel size (nm)','Position',[xpos1,top-17*vsep,xw*1.5,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.pixelsize=uicontrol('style','edit','String','100','Position',[xpos1+xw*1.5,top-17*vsep,xw*.5,vsep],'FontSize',fontsize);
            
            
            obj.guihandles.preview=uicontrol('style','pushbutton','String','Preview frame:','Position',[xpos1,top-19*vsep,xw*1.5,vsep],'FontSize',fontsize, 'Callback',@obj.preview_callback);
%             obj.guihandles.dz.TooltipString=sprintf('Distance in nm between frames. By convention, these are objective positions (not corrected for refractive index mismatch). \n A spacing between 10 nm and 50 nm works well ');
%             obj.guihandles.dzt.TooltipString=obj.guihandles.dz.TooltipString;
            obj.guihandles.previewframe=uicontrol('style','edit','String','1','Position',[xpos1+1.5*xw,top-19*vsep,xw*.5,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
            obj.guihandles.localize=uicontrol('style','pushbutton','String','Localize','Position',[xpos1+2.5*xw,top-19*vsep,xw*1.5,vsep],'FontSize',fontsize, 'Callback',@obj.localize_callback,'FontWeight','bold');
            
%             obj.guihandles.openexternal=uicontrol('style','edit','String','1','Position',[xpos1+1.5*xw,top-16*vsep,xw*.5,vsep],'FontSize',fontsize,'HorizontalAlignment',ha);
%             obj.guihandles.openexternal=uicontrol('style','pushbutton','String','Open in','Position',[xpos1,top-18*vsep,xw*1.5,vsep],'FontSize',fontsize, 'Callback',@obj.openexternal_callback,'FontWeight','bold');
%                        
%             
            obj.guihandles.status=uicontrol('style','text','String','Status','Position',[xpos1,top-21*vsep,xw*4,vsep],'FontSize',fontsize);
            
      
        end
        function selectfiles_callback(obj,a,b,which)
            switch which
                case 1
                    switch obj.guihandles.loader.Value
                        case {1,2}
                            file=obj.guihandles.imagefile.String;
                            if isempty(file)
                                file='*.tif';
                            end
                            handle='imagefile';
                        case 3 %Fiji
                            if exist('settings.mat','file')
                                l=load('settings.mat');
                                fijipath=l.fijipath;
                            end
                            if ~exist('fijipath','var')
                                [file,fijipath]=uigetfile('Image*.exe','Select the ImageJ executable in the Fiji.app directory')
                                obj.fijipath=fijipath;
                                save('settings.mat','fijipath')
                            end
                            if ~isempty(obj.mij)
                                try
                                obj.mij.exit;
                                catch err
                                    err
                                end
                            end
                            wf=msgbox('Open the raw image file in ImageJ. Then set parameters, load 3D calibratin and use preview/localize. Don''t close ImageJ.');
                            waitfor(wf);
%                             addpath([fijipath filesep 'scripts'])
%                             path=pwd;
                            myMiji(true,fijipath);
%                             cd(path)
                            obj.mij=MIJ;
                             obj.guihandles.imagefile.String='from_ImageJ';
                            return
                    end
                case 2
                    file=obj.guihandles.calfile.String;
                    if isempty(file)
                        file='*_3dcal.mat';
                    end
                    handle='calfile';
            end
            [f,p]=uigetfile(file);
            if f
                obj.guihandles.(handle).String=[p f];
                if which==1
                    if isempty(obj.guihandles.outputfile.String)
                        obj.guihandles.outputfile.String=strrep([p f],'.tif','.csv');
                    end
                end
            end
 
        end
        function selectoutput_callback(obj,a,b)
            of=obj.guihandles.outputfile.String;
            if isempty(of)
           
                of='out.csv';
            end
            [f,p]=uiputfile(of);
            if f
            obj.guihandles.outputfile.String=[p,f];
            end
        end
        
        
        function preview_callback(obj,a,b)
            p=obj.getguiparameters;
            p.preview=true;
            simplefitter_cspline(p)
        end
        
        
        function localize_callback(obj,a,b)
            p=obj.getguiparameters;
            p.preview=false;
            simplefitter_cspline(p)
        end        
        
        function p=getguiparameters(obj)
            p.imagefile=obj.guihandles.imagefile.String;
            p.calfile=obj.guihandles.calfile.String;
            p.offset=str2double(obj.guihandles.offset.String);
            p.conversion=str2double(obj.guihandles.conversion.String);
            p.previewframe=str2double(obj.guihandles.previewframe.String);
            p.peakfilter=str2double(obj.guihandles.peakfilter.String);
            p.peakcutoff=str2double(obj.guihandles.peakcutoff.String);
            p.roifit=str2double(obj.guihandles.roifit.String);
            p.bidirectional=(obj.guihandles.bidirectional.Value);
            p.mirror=(obj.guihandles.mirror.Value);
            p.status=obj.guihandles.status;
            p.outputfile=obj.guihandles.outputfile.String;
            p.outputformat=obj.guihandles.outputformat.String{obj.guihandles.outputformat.Value};
            p.pixelsize=str2num(obj.guihandles.pixelsize.String);
            p.loader=(obj.guihandles.loader.Value);
            p.mij=obj.mij;
        end
        
        function changeloader_callback(obj,a,b)
        end
       
    end
end
