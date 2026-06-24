% main script for optogenetic stimulation run analysis 
function runThis_functionSimpleNPCorr(baseName,frames,nbinning,sMatsize,lineForSetParamsAK,sFOV,breakPoint,framePeriod,find_pulses,doNP_correct)
%runThis_functionSimpleNPCorr('D:\Data\mouse393\aligned-mouse393-004\set4\mouse393-004_Cycle00001_Ch2_',[],1,512,'mouse',331,[],[],1)

% baseName - floder name and first part of tiff file names
% frames - can specify which frames to include if you don't want the whole run
% nbinning - spatial binning factor
% sMatsize - size of the FOV in pixels
% lineForSetParamsAK - setting for cell size for automated cell finding - not using this method anymore
% sFOV - size of the FOV in microns
% 'breakpoint' - a breakpoint to restrict analysis - see code
% 'framePeriod' - can add manually instead of finding with xml file
% 'find_pulses' - can be set to 1 and it will use the algorithm to find them - or it can be an array of pulse times which will be used instead - e.g. [1000 1400 1800]
% 'doNP_correct' - set to '1' to do the neuropil correction - even if '1', both corrected and uncorrected data are saved for later

if nargin<7 || isempty(breakPoint),    breakPoint=0;  end
if nargin<8, framePeriod=[];  end
if nargin<9, find_pulses=0;  end
if nargin<10, doNP_correct=1;  end

gridmasks=0;        %set to 1 to make a grid of cell masks automatically - define grid parameters in makemaskgrid.m
baselinefraction = 10 ; % the baseline is the mean of smallest 1/n of the points.
%pct_for_hist=80; % in pct units

s(1).MATsize=sMatsize;
s(1).expdir=lineForSetParamsAK;
s(1).FOV=sFOV;

p=0;

% Create directory for analysis result
%Input_fname='D:\KaraLab\Cat060604\T-04-Jun-2006-11-04-004\T-04-Jun-2006-11-04-004_Cycle001_Ch2_';
ss=strfind(baseName,'\');
mainDir=baseName(1:ss(end))

fName=baseName(ss(end)+1:end);
tmpDir=[mainDir 'analyzed'];    

analysisDir=replace(tmpDir,'D:\Data','E:\Analysis');%12/7/23 - changed since new folder organizing scheme
if ~exist(analysisDir, 'dir'),    mkdir(analysisDir);  end
maskDir=replace(analysisDir,'analyzed','masks');%12/7/23 - changed since new folder organizing scheme
if ~exist(maskDir, 'dir'),    mkdir(maskDir);  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%read tif files
matFile=[analysisDir '\stack.mat'];
if exist(matFile)
    disp(['reading file: ',matFile]);
    load(matFile,'stack');
else
    disp('Making stack_ave');
    [stack, xsize, ysize]= read_stack_no_avg (mainDir,baseName,frames,nbinning);         
    save('-v7.3',matFile,'stack');
end

if breakPoint==1
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FOV is avg of all images
matFile=[maskDir '\FOV.mat'];
if exist(matFile)
    disp(['reading file: ',matFile]);
    load(matFile,'FOV');
    h=figure;
    imshow(FOV./max(max(FOV)));
    title('FOV');
else
    disp('Calculating FOV');
    FOV=mean(stack,3);  %FOV is average of all frames
    save(matFile, 'FOV');
    h=figure;
    imshow(FOV./max(max(FOV)));
    title('FOV');
    imwrite(FOV./max(max(FOV)), [maskDir '\FOV.tif']);
    saveas(h,[maskDir '\FOV.fig']);
    imwrite(FOV./max(max(FOV)), [analysisDir '\FOV.tif']);
    saveas(h,[analysisDir '\FOV.fig']);
end

if breakPoint==2
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run find_cellsAK

avtmp1=uint8(round(255*FOV/(max(FOV(:)))));  %normalize to uint8
matFileL=[maskDir '\labelimg.mat'];
matFileB=[maskDir '\binarymask.mat'];

if gridmasks
    [labelimg,binarymask] = makemaskgrid(avtmp1);
    save(matFileL, 'labelimg');
    save(matFileB, 'binarymask');
elseif exist(matFileL)
    disp(['reading file: ',matFileL]);
    load(matFileL,'labelimg');
    if exist(matFileB)
        load(matFileB,'binarymask');
    else
        binarymask=zeros(size(labelimg));
        binarymask(labelimg>0)=1;
        save(matFileB, 'binarymask');
    end
elseif exist(matFileB)
    disp(['reading file: ',matFileB]);
    load(matFileB,'binarymask');
    labelimg=bwlabel(binarymask,4);
    save(matFileL, 'labelimg');
else
    nf=1;
    [labelimg,binarymask] = find_cellsAK(avtmp1,s,nf);
    save(matFileL, 'labelimg');
    save(matFileB, 'binarymask');
end


%***************************************************
% save shaded cell imgs
total_cells=max(max(labelimg));
shadeimg = shadecellsCR(FOV, binarymask,1);
h=figure;
imshow(shadeimg);
title('FOV Shaded');
warning off all
imwrite(shadeimg, [maskDir '\FOVShaded.tif']);
imwrite(shadeimg, [analysisDir '\FOVShaded.tif']);
warning on all
hold on;
for i=1:total_cells
    [y,x] = find(labelimg==i);
    text(x(1),y(1),num2str(i),'Color','r');
end
saveas(h,[maskDir '\FOVShadedNumbered.fig']);
saveas(h,[maskDir '\FOVShadedNumbered.tif']);
saveas(h,[analysisDir '\FOVShadedNumbered.fig']);
saveas(h,[analysisDir '\FOVShadedNumbered.tif']);

clear shadeimg
if breakPoint==5
    return
end

%%%%%%%%%% xml file %%%%%%%%%%
if nargin<8 || isempty(framePeriod)
    xmlfiles = dir([mainDir '*.xml']);
    if ~isempty(xmlfiles)
        useXML = [mainDir xmlfiles(1).name];
        xmldata=populatexmldatafromnvgui(useXML);
        if length(unique(xmldata.Cycle))>1
            if length(unique(xmldata.Cycle))>length(unique(xmldata.ImageNum))
                if mod(length(xmldata.AbsTime),length(unique(xmldata.ImageNum)))>0%if not each Image from the last Cycle is collected because run was stopped -%if mod(length(xmldata.AbsTime),2)>0 %in case more than two depth planes
                    missing_fr=(length(unique(xmldata.ImageNum))*length(unique(xmldata.Cycle)))-length(xmldata.AbsTime);
                    len_at=length(xmldata.AbsTime);
                    for ii=1:missing_fr
                        xmldata.AbsTime(len_at+ii)=xmldata.AbsTime(len_at+ii-1)+(xmldata.AbsTime(2)-xmldata.AbsTime(1));
                    end
                end                
                %if mod(length(xmldata.AbsTime),2)>0,    xmldata.AbsTime(end+1)=xmldata.AbsTime(end)+(xmldata.AbsTime(2)-xmldata.AbsTime(1));    end
                at_by_set=reshape(xmldata.AbsTime,[length(unique(xmldata.ImageNum)),length(unique(xmldata.Cycle))])';%re-arrange at so each column is a single Cycle's times
            end
        end
        if exist('at_by_set')
            framePeriod=at_by_set(2,1)-at_by_set(1,1);
        else %for single depth plane (only 1 Cycle for whole run) runs
            framePeriod=xmldata.AbsTime(2,1)-xmldata.AbsTime(1,1);
        end
    else
        disp(['Unable to open ' useXML ' ,just using frame numbers instead of time']);
    end
end
%************************************************************
%calculate t courses for cells
sz=size(stack);
nbaseline = round(sz(3)/baselinefraction); % see baselinefraction, above

v1=reshape(labelimg,[],1);   
v3=reshape(stack,[],sz(3));
for ii = 1:max(labelimg(:))               
    v2=find(labelimg==ii);        
    timetmp = mean(v3(v2,:));%timetmp = mean(v3(v1==ii,:));%timetmp = sum(v3(v1==ii,:));%timetmp = sum(stack(labelimg==ii,:),3);    
    sorttmp = sort(timetmp(:));
    xbaseline = mean(sorttmp(1:nbaseline));%xbaseline = sum(sorttmp(1:nbaseline)) /nbaseline;  % baseline, to get fractional change
    timecoursesUnCorr(:,ii) = timetmp(:)/xbaseline;%normalized but not NP corrected
    timecoursesraw(:,ii) = timetmp; %this is now the mean for each cell - so takes into account the # of pixels
    %tc_raw_sizescale(:,ii)=timetmp./length(v2);%scale raw values by number of pixels in cell mask - raw is just total luminance
end

%%%%%% baseline normalization for NPcorrected and for NP signal itself
if doNP_correct
    outstruct=SubtrNeuropil(baseName,maskDir,s,stack,baselinefraction);    
    timecourses=outstruct.timecoursesC;%cell and neuropil tcourses are separately normalized then the subtraction - a 10% change in both is more relevant than the absolute change since the NP is usually dimmer 
    tcCnorm=outstruct.tcCnorm;
else
    timecourses=timecoursesUnCorr;
    tcCnorm=timecoursesUnCorr;
end

%all pixels that aren't part of cell masks - might be useful for overall luminance issues
v4=find(labelimg==0);
noncell_timetmp = mean(v3(v4,:));%noncell_timetmp=sum(v3(v1==0,:));    
sorttmp = sort(noncell_timetmp(:));
xbaseline = sum(sorttmp(1:nbaseline)) /nbaseline;  % baseline, to get fractional change
timecourses_noncell(:,1) = noncell_timetmp(:)/xbaseline;
timecoursesraw_noncell(:,1) = noncell_timetmp;% tc_raw_sizescale_noncell(:,1)=noncell_timetmp./length(v4);%scale raw values by number of pixels in cell mask - raw is just total luminance

%all pixels that are part of cell masks
v5=find(labelimg>0);
allcell_timetmp = mean(v3(v5,:));%allcell_timetmp=sum(v3(v1>0,:));
sorttmp = sort(allcell_timetmp(:));
xbaseline = sum(sorttmp(1:nbaseline)) /nbaseline;  % baseline, to get fractional change
timecourses_allcell(:,1) = allcell_timetmp(:)/xbaseline;
timecoursesraw_allcell(:,1) = allcell_timetmp; %tc_raw_sizescale_allcell(:,1)=allcell_timetmp./length(v5);%scale raw values by number of pixels in cell mask - raw is just total luminance

if exist('framePeriod')
    times=[1:size(timecoursesUnCorr,1)]*framePeriod; 
    x_ax_lab='Time (s)';
else
    times=1:size(timecoursesUnCorr,1);
    x_ax_lab='Frames';
end
        
%%%%%%% plotting time 
f1=figure; hold on
%plot(times,timecourses);
plot(1:size(timecoursesUnCorr,1),timecoursesUnCorr);
xlabel(x_ax_lab); ylabel('fluorescence');
legend; title(analysisDir);
saveas(f1,[analysisDir '\tcourse.eps'], 'psc2');
saveas(f1,[analysisDir '\tcourse.fig']);
saveas(f1,[analysisDir '\tcourse.jpg'], 'jpeg');


f2=figure; hold on
%plot(times,timecourses);
plot(1:size(timecourses,1),timecourses);
xlabel(x_ax_lab); ylabel('fluorescence');
legend; title([analysisDir ' tc_Corrected']);
saveas(f2,[analysisDir '\tcourseC.eps'], 'psc2');
saveas(f2,[analysisDir '\tcourseC.fig']);
saveas(f2,[analysisDir '\tcourseC.jpg'], 'jpeg');

if find_pulses
    if find_pulses==1 %this is when the frames aren't given and the algorithm is supposed to find them
        if exist([analysisDir '\opto_pulses.mat'])
            load([analysisDir '\opto_pulses.mat']);
        else
            stimframes=find_opto_pulses(baseName,analysisDir);%,dataframes,optostimframes)
            save([analysisDir '\opto_pulses.mat'],'stimframes');
        end
    else %this is for when the actual pulse frames are given as the parameter
        stimframes=find_opto_pulses(baseName,analysisDir,[],find_pulses);
        save([analysisDir '\opto_pulses.mat'],'stimframes');
    end
    
    timecourses_interp=timecourses;%this is now the normalized to baseline and then NP corrected tcourses
    timecoursesUnCorr_interp=timecoursesUnCorr;
    %timecoursesSurr_interp=outstruct.timecoursesS;
    tc_raw_interp=timecoursesraw;
    timecourses_noncell_interp=timecourses_noncell;
    timecourses_allcell_interp=timecourses_allcell;
    tcCnorm_interp=tcCnorm;%NP corrected but normalized after subtraction - never liked this as much - seems you should subtract the same % response  - 10% of NP not absolute brightness of NP
    
    for ii=1:length(stimframes)
        if stimframes(ii)<size(timecourses,1)
            nonstimF1=0; nonstimF2=0; ind1=1; ind2=1;
            while ~nonstimF1 || ~nonstimF2
                if ~ismember(stimframes,stimframes(ii)-ind1),   nonstimF1=1;
                else,                                           ind1=ind1+1;
                end
                if ~ismember(stimframes,stimframes(ii)+ind2),   nonstimF2=1;
                else,                                           ind2=ind2+1;
                end
            end
            for jj=1:size(timecourses,2)
                timecourses_interp(stimframes(ii),jj)=mean([timecourses(stimframes(ii)-ind1,jj) timecourses(stimframes(ii)+ind2,jj)]);                
                tcCnorm_interp(stimframes(ii),jj)=mean([tcCnorm(stimframes(ii)-ind1,jj) tcCnorm(stimframes(ii)+ind2,jj)]);                
                timecoursesUnCorr_interp(stimframes(ii),jj)=mean([timecoursesUnCorr(stimframes(ii)-ind1,jj) timecoursesUnCorr(stimframes(ii)+ind2,jj)]);                                
                %timecoursesSurr_interp(stimframes(ii),jj)=mean([outstruct.timecoursesS(stimframes(ii)-ind1,jj) outstruct.timecoursesS(stimframes(ii)+ind2,jj)]);                                                
                tc_raw_interp(stimframes(ii),jj)=mean([timecoursesraw(stimframes(ii)-ind1,jj) timecoursesraw(stimframes(ii)+ind2,jj)]);
            end
            timecourses_noncell_interp(stimframes(ii),1)=mean([timecourses_noncell(stimframes(ii)-ind1,1) timecourses_noncell(stimframes(ii)+ind2,1)]);
            timecourses_allcell_interp(stimframes(ii),1)=mean([timecourses_allcell(stimframes(ii)-ind1,1) timecourses_allcell(stimframes(ii)+ind2,1)]);
        end
    end
    
    %tc_interp_Corr=timecoursesUnCorr_interp-outstruct.contamination*timecoursesSurr_interp;%corrected timcourse by raw fluorescence subtraction
    
    f3=figure; hold on %plot the uncorrected but interpolated tcourse - this is baseline normalized
    %plot(times,timecourses);
    plot(1:size(timecourses,1),timecoursesUnCorr_interp);
    xlabel(x_ax_lab); ylabel('fluorescence');
    legend; title([analysisDir ' tc_UnCorrected_interp']);
    saveas(f3,[analysisDir '\tcourseUnCorr_interp.eps'], 'psc2');
    saveas(f3,[analysisDir '\tcourseUnCorr_interp.fig']);
    saveas(f3,[analysisDir '\tcourseUnCorr_interp.jpg'], 'jpeg');            
    
    %%%%%% plot the NP corrected and interpolated tcourse
    f4=figure; hold on
    %plot(times,timecourses);
    plot(1:size(timecourses,1),timecourses_interp);
    xlabel(x_ax_lab); ylabel('fluorescence');
    legend; title([analysisDir ' tc_Corrected_interp']);
    saveas(f4,[analysisDir '\tcourseCinterp.eps'], 'psc2');
    saveas(f4,[analysisDir '\tcourseCinterp.fig']);
    saveas(f4,[analysisDir '\tcourseCinterp.jpg'], 'jpeg');            
    
    %%%%%% plot the interpolated first and then NP corrected tcourse
    %This is identical to timecourses_interp - first NP correctd and then interp
    %{
    f5=figure; hold on
    %plot(times,timecourses);
    plot(1:size(timecourses,1),tc_interp_Corr);
    xlabel(x_ax_lab); ylabel('fluorescence');
    legend; title([analysisDir ' tc_interp_Corrected']);
    saveas(f5,[analysisDir '\tcourse_interp_Corr.eps'], 'psc2');
    saveas(f5,[analysisDir '\tcourse_interp_Corr.fig']);
    saveas(f5,[analysisDir '\tcourse_interp_Corr.jpg'], 'jpeg');                
    %}
    save([analysisDir '\tc_data_interp.mat'],'timecourses_interp','tcCnorm_interp','timecoursesUnCorr_interp','tc_raw_interp','timecourses_noncell_interp','timecourses_allcell_interp');
end

%%%%%%%%% save data
tcourseFile=[analysisDir '\tc_data.mat'];% the 'timecourses' variable is now the NP corrected tcourse
if doNP_correct
    save(tcourseFile, 'times','timecourses','tcCnorm','timecoursesUnCorr','timecoursesraw','timecourses_noncell','timecoursesraw_noncell','timecourses_allcell','timecoursesraw_allcell','outstruct');
else
    save(tcourseFile, 'times','timecourses','tcCnorm','timecoursesUnCorr','timecoursesraw','timecourses_noncell','timecoursesraw_noncell','timecourses_allcell','timecoursesraw_allcell');
end

xlsfile=[analysisDir '\timecourses.xlsx'];
xlswrite(xlsfile, cellstr('time'), 'TcourseNPcorrNorm', 'A1');
xlswrite(xlsfile, times', 'TcourseNPcorrNorm', 'A2');
xlswrite(xlsfile, timecourses, 'TcourseNPcorrNorm', 'B2');
xlswrite(xlsfile, cellstr('time'), 'TcourseRaw', 'A1');
xlswrite(xlsfile, times', 'TcourseRaw', 'A2');
xlswrite(xlsfile, timecoursesraw, 'TcourseRaw', 'B2');
xlswrite(xlsfile, cellstr('time'), 'TcourseOrigNorm', 'A1');
xlswrite(xlsfile, times', 'TcourseOrigNorm', 'A2');
xlswrite(xlsfile, timecoursesUnCorr, 'TcourseOrigNorm', 'B2');
      
end