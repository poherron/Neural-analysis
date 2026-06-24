% main script for visual orientation analysis 
function p=runThis_function(baseName,nframes_per_stim,nstim_per_run,nbinning,run_inds,base_inds,stim_inds,close_fig_flag,sMatsize,lineForSetParamsAK,sFOV,breakPoint,skip_opto)
%runThis_function('D:\Data\mouse292\aligned-mouse292-016\set4\mouse292-016_Cycle00001_Ch2_',60,8,1,[0:5],[40:50],[51:60],1,512,'mouse',441,0,0 

% baseName - floder name and first part of tiff file names
% nframes_per_stim - how many frames for a single orientation stimulus presentation (blank + stim)
% nstim_per_run - number of stimuli (8 orientations)
% nbinning - spatial binning factor
% run_inds - run indexes to be analyzed - zero based regadless of file names base
% base_inds - baseline frames in one trial 
% stim_inds - stimulus frames in one trial
% close_fig_flag - close figures after rounds of plotting
% sMatsize - size of the FOV in pixels
% lineForSetParamsAK - setting for cell size for automated cell finding - not using this method anymore
% sFOV - size of the FOV in microns
% breakPoint - to break the analysis at different points
% skip_opto - set to 1 to interpolate data across frames where opto pulses happened 

if nargin<12 || isempty(breakPoint),    breakPoint=0;    end
if nargin<13 || isempty(skip_opto),     skip_opto=0;    end


runThis_function_internal(baseName,nframes_per_stim,nstim_per_run,nbinning,run_inds,base_inds,stim_inds,close_fig_flag,sMatsize,lineForSetParamsAK,sFOV,breakPoint,skip_opto);

end

function p=runThis_function_internal(baseName,nframes_per_stim,nstim_per_run,nbinning,run_inds,base_inds,stim_inds,close_fig_flag,sMatsize,lineForSetParamsAK,sFOV,breakPoint,skip_opto)


gridmasks=0;        %set to 1 to make a grid of cell masks automatically - define grid parameters in makemaskgrid.m
%print_flag controls printing of figures
master_print_flag=0;
print_flag_pix_maps=1;
print_flag_allTc=1;
print_flag_allTc_with_dots=1;
print_flag_Tc_pages=1;
print_flag_fit=0;
print_flag_allTc_with_dots_sliding_bkg=1;
print_flag_Tc_pages_sliding_bkg=1;
print_flag_cell_maps=1;
print_flag_scatter=0;
print_flag_hist=0;

print_flag_pix_maps=print_flag_pix_maps*master_print_flag;
print_flag_allTc=print_flag_allTc*master_print_flag;
print_flag_allTc_with_dots=print_flag_allTc_with_dots*master_print_flag;
print_flag_Tc_pages=print_flag_Tc_pages*master_print_flag;
print_flag_fit=print_flag_fit*master_print_flag;
print_flag_allTc_with_dots_sliding_bkg=print_flag_allTc_with_dots_sliding_bkg*master_print_flag;
print_flag_Tc_pages_sliding_bkg=print_flag_Tc_pages_sliding_bkg*master_print_flag;
print_flag_cell_maps=print_flag_cell_maps*master_print_flag;
print_flag_scatter=print_flag_scatter*master_print_flag;
print_flag_hist=print_flag_hist*master_print_flag;


% master_alpha(1) - alpha resp, master_alpha(2) - alpha sel
master_alpha=[0.05 0.05];
pct_for_hist=80; % in pct units
show_R2_flag=1;
nframes_per_run = nframes_per_stim .* nstim_per_run;

s(1).MATsize=sMatsize;
s(1).expdir=lineForSetParamsAK;
s(1).FOV=sFOV;

p=0;

% stack - images averaged across multiple runs
global stack;

% Create directory for analysis result
ss=strfind(baseName,'\');
mainDir=baseName(1:ss(end))
tmpDir=[mainDir 'analyzed'];
analysisDir=replace(tmpDir,'Data','Analysis'); % new folder organizing scheme
if ~exist(analysisDir, 'dir'),    mkdir(analysisDir);  end
maskDir=replace(analysisDir,'analyzed','masks');% new folder organizing scheme
if ~exist(maskDir, 'dir'),    mkdir(maskDir);  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%read tif files
matFile=[analysisDir '\stack.mat'];
if exist(matFile)
    disp(['reading file: ',matFile]);
    load(matFile,'stack');
else
    disp('Making stack_ave');
    [stack, xsize, ysize]= read_stack_ave (baseName, nframes_per_run, run_inds, nbinning);         
    save(matFile, 'stack');
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
    disp('Calculating VOF');
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

if close_fig_flag,    close (h);    end

if breakPoint==2
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate and save sliding base

%create spatial filter. Will be used to calculate smoothed (_sm) images
lowpass = 1; % (um) normally its around 1 or 2 microns.
lowpass_pix = lowpass * s(1).MATsize / s(1).FOV / nbinning;
kernel_size = ceil(lowpass_pix * 5);
sp_filter=fspecial('gaussian', kernel_size, lowpass_pix);    % spatial lowpass filter in pixel unit (after binning)

matFile=[analysisDir '\base.mat'];
if exist(matFile)
    disp(['reading file: ',matFile]);
    load(matFile,'base');
    matFile=[analysisDir '\base_sm.mat'];
    load(matFile,'base_sm');
else
    disp('Calculating base');
    [base, base_sm] = calc_base_sliding(base_inds, nframes_per_stim, sp_filter);
    save(matFile, 'base');
    matFile=[analysisDir '\base_sm.mat'];
    save(matFile, 'base_sm');
end

if breakPoint==3
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate and save dF
% if spatial filter is specified then dir_ratio and ori_ratio are smoothed
matFile=[analysisDir '\dir_dF.mat'];
if exist(matFile)
    disp(['reading file: ',matFile]);
    load(matFile,'dir_dF');
    matFile=[analysisDir '\dir_dF_sm.mat'];
    load(matFile,'dir_dF_sm');
    matFile=[analysisDir '\ori_dF.mat'];
    load(matFile,'ori_dF');
    matFile=[analysisDir '\ori_dF_sm.mat'];
    load(matFile,'ori_dF_sm');
    matFile=[analysisDir '\dir_ratio.mat'];
    load(matFile,'dir_ratio');
    matFile=[analysisDir '\ori_ratio.mat'];
    load(matFile,'ori_ratio');
else
    disp('Calculating dF');
    [dir_dF, dir_dF_sm, ori_dF, ori_dF_sm, dir_ratio, ori_ratio] = calc_dF_sliding_base(base, stim_inds, nframes_per_stim, nstim_per_run, sp_filter);
    save(matFile, 'dir_dF');
    matFile=[analysisDir '\dir_dF_sm.mat'];
    save(matFile, 'dir_dF_sm');
    matFile=[analysisDir '\ori_dF.mat'];
    save(matFile, 'ori_dF');
    matFile=[analysisDir '\ori_dF_sm.mat'];
    save(matFile, 'ori_dF_sm');
    matFile=[analysisDir '\dir_ratio.mat'];
    save(matFile, 'dir_ratio');
    matFile=[analysisDir '\ori_ratio.mat'];
    save(matFile, 'ori_ratio');
end
clear stack base base_sm;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save dF images
tifFile=[analysisDir '\dir_dF_1.tif'];
if exist(tifFile)
    disp('Files dir_dF_N.tif ... are already in analysis directory');
else
    [max_dir_dF, max_dir_ratio_change, max_ori_dF, max_ori_ratio_change] = write_dF_images_SY (analysisDir,dir_dF_sm, ori_dF_sm, dir_ratio, ori_ratio);

    disp(['Doing map_params']);
    dir_dF_params = calc_map_params (dir_dF_sm);
    dir_ratio_params = calc_map_params (dir_ratio);
    ori_dF_params = calc_map_params (ori_dF_sm);
    ori_ratio_params = calc_map_params (ori_ratio);
    save ([analysisDir '\dir_dF_params'], 'dir_dF_params');
    save ([analysisDir '\dir_ratio_params'], 'dir_ratio_params');
    save ([analysisDir '\ori_dF_params'], 'ori_dF_params');
    save ([analysisDir '\ori_ratio_params'], 'ori_ratio_params');
    write_intensity_maps (dir_dF_params, [analysisDir '\dir_dF_'], max_dir_dF);
    write_intensity_maps (dir_ratio_params, [analysisDir '\dir_ratio_'], max_dir_ratio_change);
    write_intensity_maps (ori_dF_params, [analysisDir '\ori_dF_'], max_ori_dF);
    write_intensity_maps (ori_ratio_params, [analysisDir '\ori_ratio_'], max_ori_ratio_change);

    disp(['Doing angle_map']);
    dir_angle = write_angle_map (dir_dF_params.th, [analysisDir '\dir_']);
    ori_angle = write_angle_map (ori_dF_params.th, [analysisDir '\ori_']);
    save ([analysisDir '\dir_angle'], 'dir_angle');
    save ([analysisDir '\ori_angle'], 'ori_angle');

    disp(['Doing polar_map']);
    dir_dF_polar = write_polar_map (dir_dF_params, [analysisDir '\dir_dF_'], max_dir_dF);
    dir_ratio_polar = write_polar_map (dir_ratio_params, [analysisDir '\dir_ratio_'], max_dir_ratio_change);
    ori_dF_polar = write_polar_map (ori_dF_params, [analysisDir '\ori_dF_'], max_ori_dF);
    ori_ratio_polar = write_polar_map (ori_ratio_params, [analysisDir '\ori_ratio_'], max_ori_ratio_change);
    save ([analysisDir '\dir_dF_polar'], 'dir_dF_polar');
    save ([analysisDir '\dir_ratio_polar'], 'dir_ratio_polar');
    save ([analysisDir '\ori_dF_polar'], 'ori_dF_polar');
    save ([analysisDir '\ori_ratio_polar'], 'ori_ratio_polar');

    tune_max = 0.4; % in HLS map, if tune >= tune_max, color will be saturated. if not, color will be unsaturated.

    disp(['Doing HLS_map']);
    dir_dF_HLS = write_HLS_map (dir_dF_params, [analysisDir '\dir_dF_'], max_dir_dF, tune_max);
    dir_ratio_HLS = write_HLS_map (dir_ratio_params, [analysisDir '\dir_ratio_'], max_dir_ratio_change, tune_max);
    ori_dF_HLS = write_HLS_map (ori_dF_params, [analysisDir '\ori_dF_'], max_ori_dF, tune_max);
    ori_ratio_HLS = write_HLS_map (ori_ratio_params, [analysisDir '\ori_ratio_'], max_ori_ratio_change, tune_max);
    save ([analysisDir '\dir_dF_HLS'], 'dir_dF_HLS');
    save ([analysisDir '\dir_ratio_HLS'], 'dir_ratio_HLS');
    save ([analysisDir '\ori_dF_HLS'], 'ori_dF_HLS');
    save ([analysisDir '\ori_ratio_HLS'], 'ori_ratio_HLS');

    write_hue_colormap(dir_dF_params, [analysisDir '\pix_maps_']);


    %log.input_file = Input_fname;
    log.baseName=baseName;
    log.nframes_per_stim = nframes_per_stim;
    log.nstim_per_run = nstim_per_run;
    log.nframes_per_trial = nframes_per_run;
    log.run_inds = run_inds;
    log.nbinning = nbinning;
    log.base_inds = base_inds;
    log.stim_inds = stim_inds;
    log.sMatsize=sMatsize;
    log.lineForSetParamsAK=lineForSetParamsAK;
    log.sFOV=sFOV;
    log.breakPoint=breakPoint;
    log.lowpass = lowpass;
    log.kernel_size = kernel_size;
    log.tune_max = tune_max;
    log.max_dir_dF = max_dir_dF;
    log.max_dir_ratio_change = max_dir_ratio_change;
    log.max_ori_dF = max_ori_dF;
    log.max_ori_ratio_change = max_ori_ratio_change;

    save ([analysisDir '\log'], 'log');
end

if print_flag_pix_maps==1
    printmapsPK(1,0,'',nstim_per_run,analysisDir);
end

clear dir_dF;
clear dir_dF_sm;
clear ori_dF;
clear ori_dF_sm;

clear dir_dF_params;
clear dir_ratio_params;
clear ori_dF_params;
clear ori_ratio_params;
clear dir_angle;
clear ori_angle;
clear dir_dF_polar;
clear dir_ratio_polar;
clear ori_dF_polar;
clear ori_ratio_polar;
clear dir_dF_HLS;
clear dir_ratio_HLS;
clear ori_dF_HLS;
clear ori_ratio_HLS;

if breakPoint==4
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

if close_fig_flag
    close all
    close all hidden
end
pause(1)

%************************************************************
%calculate t courses for cells
matFile=[analysisDir '\tcourse.mat'];
if exist(matFile)
    disp(['reading file: ',matFile]);
    load(matFile);
else
    disp('Calculating T courses');    
    tc = get_tcourses_from_files(labelimg, baseName, nframes_per_run, run_inds,nbinning);
    
    %%%% interpolate data points during opto stim when PMT shutter closes -
    %%%% can load it or try to automatically compute it
    if skip_opto
        if exist([analysisDir '\opto_pulses.mat'])
            load([analysisDir '\opto_pulses.mat']); 
        else
            stimframes=find_opto_pulses(baseName,analysisDir);%,dataframes,optostimframes)
            save([analysisDir '\opto_pulses.mat'],'stimframes');
        end    
        if breakPoint==6
            return
        end
        tc_orig=tc;
        for ii=1:length(stimframes)
            for jj=1:size(tc,2)
                if stimframes(ii)<size(tc,1)                
                tc(stimframes(ii),jj)=mean([tc_orig(stimframes(ii)-1,jj) tc_orig(stimframes(ii)+1,jj)]);    
                end
            end
        end
        [tcavg,tcavgnorm,tcnorm,cellavg,tcnormstd,tcavgnormstd]=tcourseprocCRSY(tc, [],nframes_per_run);
        save(matFile,'tc_orig','tc','tcavg','tcavgnorm','tcnorm','cellavg','tcnormstd','tcavgnormstd');
    else
        [tcavg,tcavgnorm,tcnorm,cellavg,tcnormstd,tcavgnormstd]=tcourseprocCRSY(tc, [],nframes_per_run);
        save(matFile,'tc','tcavg','tcavgnorm','tcnorm','cellavg','tcnormstd','tcavgnormstd');
    end
        
end


%save as colorcoded figure
h20=figure;
imagesc(tcavgnorm'-1);
axis xy;
xlabel('Frame #');
ylabel('Cell #');
title(strrep(analysisDir,'\','\\'));
colorbar;
if print_flag_allTc==1
    print;
end
saveas(h20,[analysisDir '\allTc_colored.fig']);

%tc's are already only for run_inds

figFile=[analysisDir '\allTc_with_dots.fig'];
if exist(figFile)
    disp('File allTc_with_dots.fig is already in directory');
else
    %plot all cells on one figure (normalized and avg over all trials)
    h121=figure;
    plot(tcavgnorm);
    saveas(h121,[analysisDir '\allTc.fig']);

    %plot avg Tc with dots
    alltc=mean(tcavgnorm ,2);
    h5=figure;
    plot(alltc);
    hold on;
    title(strrep(analysisDir,'\','\\'));
    for r=1:nstim_per_run
        plot((r-1)*nframes_per_stim+base_inds,alltc((r-1)*nframes_per_stim+base_inds)+0.01,'ok');
        plot((r-1)*nframes_per_stim+stim_inds,alltc((r-1)*nframes_per_stim+stim_inds)+0.01,'or');
    end
    grid;
    %##################### shading
    hold on
    ax=axis;
    for ngray=1:nstim_per_run
        rectangle('Position',[nframes_per_stim*(ngray-1)+stim_inds(1),ax(3),length(stim_inds),ax(4)-ax(3)],'FaceColor',[0.95 0.95 0.95])
    end
    %replot
    plot(alltc);
    for r=1:nstim_per_run
        plot((r-1)*nframes_per_stim+base_inds,alltc((r-1)*nframes_per_stim+base_inds)+0.01,'ok');
        plot((r-1)*nframes_per_stim+stim_inds,alltc((r-1)*nframes_per_stim+stim_inds)+0.01,'or');
    end
    axis(ax);
    %####################
    set(gca, 'Layer', 'top');


    if print_flag_allTc_with_dots==1
        print;
    end
    saveas(h5,[analysisDir '\allTc_with_dots.fig']);

    frrr=1:nframes_per_stim;
    for r=1:nstim_per_run/8
        h6=figure;
        set(h6,'PaperUnits','normalized');
        set(h6,'PaperPosition',[0.05 0.05 0.9 0.9]);

        for w=1:8
            subplot(4,2,w);
            plot(alltc(((r-1)*8+w-1)*nframes_per_stim+frrr));
            hold on;
            title([strrep(analysisDir,'\','\\') ' ori # ' num2str((r-1)*4+w)]);
            plot(base_inds,alltc(((r-1)*8+w-1)*nframes_per_stim+base_inds)+0.01,'ok');
            plot(stim_inds,alltc(((r-1)*8+w-1)*nframes_per_stim+stim_inds)+0.01,'or');
            grid;

            %##################### shading
            hold on
            ax=axis;
            for ngray=1:nstim_per_run
                rectangle('Position',[nframes_per_stim*(ngray-1)+stim_inds(1),ax(3),length(stim_inds),ax(4)-ax(3)],'FaceColor',[0.95 0.95 0.95])
            end
            %replot
            plot(alltc(((r-1)*8+w-1)*nframes_per_stim+frrr));
            plot(base_inds,alltc(((r-1)*8+w-1)*nframes_per_stim+base_inds)+0.01,'ok');
            plot(stim_inds,alltc(((r-1)*8+w-1)*nframes_per_stim+stim_inds)+0.01,'or');
            axis(ax);
            %####################
            set(gca, 'Layer', 'top');

        end
        if print_flag_allTc_with_dots==1
            print;
        end
        saveas(h6,[analysisDir '\allTc_with_dots_' num2str(r) '.fig']);
        hold off;
    end
end



%****************************************
%print t courses
disp('Printing T courses');
sz=size(tcnorm);
n_perpage=14;
npage=ceil(sz(2)/n_perpage);
n=1;

warning off all
for i=1:npage
    h=figure;
    set(h,'PaperUnits','normalized');
    set(h,'PaperPosition',[0.05 0.05 0.9 0.9]);
    %main comment

    subplot('Position',[0.05 0.95 0.9 0.05]);
    set(gca,'FontSize',6);
    axis off;
    text(0.1,0.8,analysisDir,'FontSize',6);

    text(0.1,0.6,[num2str(i) '/' num2str(npage)],'FontSize',6);
    % image
    subplot('Position',[0.05 0.8 0.2 0.14]);
    set(gca,'FontSize',6);
    %axis square;
    image(FOV*64/max(max(FOV)));
    colormap gray;

    ker=fspecial('gaussian',90,30);
    avg_f=filter2KO(ker,FOV);
    avg_img_2=FOV./avg_f;
    avg_img_2=avg_img_2-min(min(avg_img_2));

    msk=im2bw(labelimg,0.0001);
    im3c=avg_img_2./max(max(avg_img_2));
    im3c(:,:,2)=avg_img_2./max(max(avg_img_2));
    im3c(:,:,3)=avg_img_2./max(max(avg_img_2));
    im3c(:,:,1)=(avg_img_2./max(max(avg_img_2))).*(1-msk);
    % image2
    subplot('Position',[0.3 0.8 0.2 0.14]);
    set(gca,'FontSize',6);
    %axis square;
    image(avg_img_2*64/max(max(avg_img_2)));
    image(im3c);

    %cells
    h=subplot('Position',[0.55 0.8 0.2 0.14]);
    set(gca,'FontSize',6);
    %axis square
    map(1:sz(2),1:3)=0.5;
    map((i-1)*n_perpage+1:min(i*n_perpage,sz(2)),1:3)=1;
    lblim=label2rgb(labelimg,map,'k');
    image(lblim);

    hold on;

    for j=1:n_perpage
        if n>sz(2)
            break;
        end

        axes(h);
        [y,x] = find(labelimg==n);
        text(x(1),y(1),num2str(n),'Color','r','FontSize',6);

        %subcomment
        subplot('Position',[0.8 0.804-j*0.056 0.15 0.05]);
        axis off;
        text(0.1,0.8,num2str(n));

        subplot('Position',[0.05 0.804-j*0.056 0.2 0.05]);
        set(gca,'FontSize',6);
        plot(tcnorm(:,n));

        set(gca,'xtick',[0 sz(1)]);
        set(gca,'xgrid','on');

        subplot('Position',[0.3 0.804-j*0.056 0.2 0.05]);
        set(gca,'FontSize',6);
        plot(tcavgnorm(:,n));

        %##################### shading
        hold on
        ax=axis;
        for ngray=1:nstim_per_run
            rectangle('Position',[nframes_per_stim*(ngray-1)+stim_inds(1),ax(3),length(stim_inds),ax(4)-ax(3)],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none')
        end
        %replot
        plot(tcavgnorm(:,n));
        axis(ax);
        %####################
        set(gca, 'Layer', 'top');
        set(gca,'xtick',[0 length(tcavgnorm(:,n))]);
        set(gca,'xgrid','on');

        subplot('Position',[0.55 0.804-j*0.056 0.2 0.05]);
        set(gca,'FontSize',6);
        szlong=size(tcnorm);
        szshort=size(tcavgnorm);

        all_shorts=reshape(tcnorm(:,n),szshort(1),szlong(1)/szshort(1));
        plot(all_shorts);
        %##################### shading
        hold on
        ax=axis;
        for ngray=1:nstim_per_run
            rectangle('Position',[nframes_per_stim*(ngray-1)+stim_inds(1),ax(3),length(stim_inds),ax(4)-ax(3)],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none')
        end
        %replot
        plot(all_shorts);
        axis(ax);
        %####################
        set(gca, 'Layer', 'top');
        set(gca,'xtick',[0 szshort(1)]);
        set(gca,'xgrid','on');

        n=n+1;
    end
    if print_flag_Tc_pages==1
        print;
    end

    % output as a postscript file

    if print_flag_Tc_pages==2
        print  %YC 09/17/04
        print ('-dpsc','-append',[matoutfname,'_tcourse_plot.ps']);
    end
end

warning on all
if close_fig_flag
    close all
end

if breakPoint==7
    return
end

%************************************************************
%do oristat
matFile=[analysisDir '\OriStat.mat'];
statFname=matFile;
if exist(matFile)
else
    disp('Doing OriStat');
    %tc's are already only for run_inds
    doOriStatKO3_SY(analysisDir,nframes_per_stim,nstim_per_run,base_inds,stim_inds,master_alpha);
end

if breakPoint==8
    return
end

disp(['reading file: ',matFile]);
load(matFile,'OriStatKO');

%***********************************
% print fitted curves
% '\estimate_dir_tuning'

disp('Printing fits');
sz=size(tcnorm);
%n_perpage2=[8 4];
n_perpage2=[8 2];
npage=ceil(sz(2)/(n_perpage2(1)*n_perpage2(2)));
n=1;

warning off all
for i=1:npage
    hf=figure;
    set(hf,'PaperUnits','normalized');
    set(hf,'PaperPosition',[0.05 0.05 0.9 0.9]);
    %main comment

    subplot('Position',[0.05 0.95 0.9 0.05]);
    set(gca,'FontSize',6);
    axis off;
    text(0.1,0.8,analysisDir,'FontSize',6);

    text(0.1,0.6,[num2str(i) '/' num2str(npage)],'FontSize',6);
    % image
    subplot('Position',[0.05 0.8 0.2 0.14]);
    set(gca,'FontSize',6);
    %axis square;
    image(FOV*64/max(max(FOV)));
    colormap gray;

    ker=fspecial('gaussian',90,30);
    avg_f=filter2KO(ker,FOV);
    avg_img_2=FOV./avg_f;
    avg_img_2=avg_img_2-min(min(avg_img_2));

    msk=im2bw(labelimg,0.0001);
    im3c=avg_img_2./max(max(avg_img_2));
    im3c(:,:,2)=avg_img_2./max(max(avg_img_2));
    im3c(:,:,3)=avg_img_2./max(max(avg_img_2));
    im3c(:,:,1)=(avg_img_2./max(max(avg_img_2))).*(1-msk);
    % image2
    subplot('Position',[0.3 0.8 0.2 0.14]);
    set(gca,'FontSize',6);
    %axis square;
    image(avg_img_2*64/max(max(avg_img_2)));
    image(im3c);

    %cells
    h=subplot('Position',[0.55 0.8 0.2 0.14]);
    set(gca,'FontSize',6);
    %axis square
    map(1:sz(2),1:3)=0.5;
    map((i-1)*n_perpage2(1)*n_perpage2(2)+1:min(i*n_perpage2(1)*n_perpage2(2),sz(2)),1:3)=1;
    lblim=label2rgb(labelimg,map,'k');
    image(lblim);

    hold on;

    for j=1:n_perpage2(1)
        for k=1:n_perpage2(2)
            if n>sz(2)
                break;
            end

            axes(h);
            [y,x] = find(labelimg==n);
            text(x(1),y(1),num2str(n),'Color','r','FontSize',6);
            %read fitted plot from file and put it on page
            if exist([analysisDir  '\estimate_dir_tuning' num2str(n) '.fig'],'file')
                of=open([analysisDir  '\estimate_dir_tuning' num2str(n) '.fig']);
                figure(of);
                a=gca;
                anew=copyobj(a,hf);
                set(anew,'Position',[0.05+(k-1)*(0.9/n_perpage2(2)-0.01) 0.75-j*(0.75/n_perpage2(1)-0.005) 0.45/n_perpage2(2)-0.03 0.8/n_perpage2(1)-0.03]);
                set(anew,'FontSize',6);
                if show_R2_flag
                    title(anew,[num2str(n)  ' DirTW=' num2str(OriStatKO(n).dir_tuning_width) ' R2=' num2str(OriStatKO(n).dir_fit_R2, '%10.2f')]);
                else
                    title(anew,[num2str(n)  ' DirTW=' num2str(OriStatKO(n).dir_tuning_width)]);
                end
                set(anew,'xgrid','on');
                v=axis(anew);
                v(1)=0;
                v(2)=360;
                axis(anew,v);
                set(anew,'YTick',sort(unique([v(3) 0  v(4)])))
                if j~=n_perpage2(1)
                    set(anew,'XTickLabel',[])
                end
                close(of);
            end

            if exist([analysisDir  '\estimate_ori_tuning' num2str(n) '.fig'],'file')
                of=open([analysisDir  '\estimate_ori_tuning' num2str(n) '.fig']);%PO 8/24/17
                figure(of);
                a=gca;
                anew=copyobj(a,hf);
                set(anew,'Position',[0.05+((k-1)+0.5)*(0.9/n_perpage2(2)-0.01) 0.75-j*(0.75/n_perpage2(1)-0.005) 0.45/n_perpage2(2)-0.03 0.8/n_perpage2(1)-0.03]);
                set(anew,'FontSize',6);
                if show_R2_flag
                    title(anew,[num2str(n) '  CV=' num2str(OriStatKO(n).cirVar, '%10.2f') ' OriTW=' num2str(OriStatKO(n).ori_tuning_width) ' R2=' num2str(OriStatKO(n).ori_fit_R2, '%10.2f')] );
                else
                    title(anew,[num2str(n) '  CV=' num2str(OriStatKO(n).cirVar, '%10.2f') ' OriTW=' num2str(OriStatKO(n).ori_tuning_width)] );
                end
                set(anew,'xgrid','on');
                v=axis(anew);
                v(1)=0;
                v(2)=180;
                axis(anew,v);
                set(anew,'YTick',sort(unique([v(3) 0  v(4)])))
                if j~=n_perpage2(1)
                    set(anew,'XTickLabel',[])
                end

                close(of);
            end

            n=n+1;
        end
    end

    if print_flag_fit==1
        print;
    end
    if print_flag_fit==2
        print 
        print ('-dpsc','-append',[matoutfname,'_fit_plot.ps']);
    end
end

warning on all
if close_fig_flag
    close all
end

%************************

disp(['reading file: ',matFile]);
load(matFile,'OriStatKO');

matFile=[analysisDir '\tcourse_sliding_bkg.mat'];
disp(['reading file: ',matFile]);
load(matFile);
disp('Processing T courses _sliding_bkg');
[tcavg_sliding_bkg,tcavgnorm_sliding_bkg,tcnorm_sliding_bkg,cellavg_sliding_bkg,tcnormstd_sliding_bkg,tcavgnormstd_sliding_bkg,...
    tcavgnorm_sliding_bkg2,tcnorm_sliding_bkg2,tcnormstd_sliding_bkg2,tcavgnormstd_sliding_bkg2]=tcourseprocCRSY(tc_sliding_bkg, tc_sliding_bkg_sub_bkg, nframes_per_run);
save(matFile,'tc_sliding_bkg','tc_sliding_bkg_sub_bkg','tc_lowcut','tc_lowcut_dF_F','tcavg_sliding_bkg','tcavgnorm_sliding_bkg','tcavgnorm_sliding_bkg2','tcnorm_sliding_bkg',...
    'tcnorm_sliding_bkg2','cellavg_sliding_bkg','tcnormstd_sliding_bkg','tcavgnormstd_sliding_bkg');

%save each Tc_sliding in separete fig. file
sz=size(tcavgnorm_sliding_bkg);

alltc_sliding_bkg=mean(tcavgnorm_sliding_bkg ,2);
h5=figure;
plot(alltc_sliding_bkg);
hold on;
title([strrep(analysisDir,'\','\\') '  sliding_bkg']);
for r=1:nstim_per_run
    plot((r-1)*nframes_per_stim+base_inds,alltc_sliding_bkg((r-1)*nframes_per_stim+base_inds)+0.01,'ok');
    plot((r-1)*nframes_per_stim+stim_inds,alltc_sliding_bkg((r-1)*nframes_per_stim+stim_inds)+0.01,'or');
end
grid;
%##################### shading
hold on
ax=axis;
for ngray=1:nstim_per_run
    rectangle('Position',[nframes_per_stim*(ngray-1)+stim_inds(1),ax(3),length(stim_inds),ax(4)-ax(3)],'FaceColor',[0.95 0.95 0.95])
end
%replot
plot(alltc_sliding_bkg);
for r=1:nstim_per_run
    plot((r-1)*nframes_per_stim+base_inds,alltc_sliding_bkg((r-1)*nframes_per_stim+base_inds)+0.01,'ok');
    plot((r-1)*nframes_per_stim+stim_inds,alltc_sliding_bkg((r-1)*nframes_per_stim+stim_inds)+0.01,'or');
end
axis(ax);
%####################
set(gca, 'Layer', 'top');

if print_flag_allTc_with_dots_sliding_bkg==1
    print;
end
saveas(h5,[analysisDir '\allTc_with_dots_sliding_bkg.fig']);



%****************************************
% print t courses sliding bkg
disp('Printing T courses sliding bkg');

sz=size(tcnorm_sliding_bkg);
n_perpage=14;
npage=ceil(sz(2)/n_perpage);
n=1;

warning off all
% COLUMN LAYOUT
% ================================
ncols = 5;  % full trace | orient avg | trials | mean stim | stats
left_margin = 0.05;
right_margin = 0.02;
col_spacing = 0.02;
usable_width = 1-left_margin-right_margin-(ncols-1)*col_spacing;
col_width = usable_width/ncols;
col = left_margin + (0:ncols-1)*(col_width+col_spacing);
% ROW LAYOUT (AUTO FIT)
% ================================
top_of_rows = 0.75; 
bottom_margin = 0.05;
usable_height = top_of_rows - bottom_margin;
row_spacing = usable_height / n_perpage;
row_height = row_spacing * 0.8;
row_start = top_of_rows;

% ================================
for i=1:npage
    h=figure;
    set(h,'Units','normalized','Position',[0.02 0.05 0.96 0.9])
    set(h,'PaperUnits','normalized');
    set(h,'PaperPosition',[0.05 0.05 0.9 0.9]);
    % HEADER
    axes('Position',[0.05 0.94 0.9 0.05]);
    axis off
    set(gca,'FontSize',6);
    text(0.1,0.8,analysisDir,'FontSize',6);
    text(0.1,0.5,[num2str(i) '/' num2str(npage) '  sliding bkg'],'FontSize',6);

    % TOP IMAGES (SQUARE)
    img_w = 0.18;
    img_h = 0.18;
    % Image 1
    axes('Position',[0.05 0.78 img_w img_h]);
    image(FOV*64/max(max(FOV)));
    axis image
    colormap gray

    % Image 2
    ker=fspecial('gaussian',90,30);
    avg_f=filter2KO(ker,FOV);
    avg_img_2=FOV./avg_f;
    avg_img_2=avg_img_2-min(min(avg_img_2));
    msk=im2bw(labelimg,0.0001);
    im3c=avg_img_2./max(max(avg_img_2));
    im3c(:,:,2)=im3c(:,:,1);
    im3c(:,:,3)=im3c(:,:,1);
    im3c(:,:,1)=im3c(:,:,1).*(1-msk);
    axes('Position',[0.28 0.78 img_w img_h]);
    image(im3c)
    axis image
    
    % Cell map
    hh=axes('Position',[0.51 0.78 img_w img_h]);
    map(1:sz(2),1:3)=0.5;
    map((i-1)*n_perpage+1:min(i*n_perpage,sz(2)),1:3)=1;
    lblim=label2rgb(labelimg,map,'k');
    image(lblim);
    axis image
    hold on
   
    % NEURON LOOP
    % ==========================
    for j=1:n_perpage
        if n>sz(2)
            break
        end
        row_y = row_start - j*row_spacing;
        axes(hh)
        [y,x] = find(labelimg==n);
        text(x(1),y(1),num2str(n),'Color','r','FontSize',6);

        % COLUMN 5 : statistics
        subplot('Position',[col(5) row_y col_width row_height]);
        axis off
        best_dir_stdmean = OriStatKO(n).dir_ratio_change_STD(OriStatKO(n).best_dir) / OriStatKO(n).dir_ratio_change(OriStatKO(n).best_dir);
        text(0.05,0.8,[num2str(n) ' best dir std/mean:' num2str(best_dir_stdmean)],'FontSize',6);
        text(0.0,0.55,'p resp','FontSize',6);
        text(0.2,0.55,'psel','FontSize',6);
        text(0.4,0.55,'pdF/F','FontSize',6);
        text(0.6,0.55,'pttsep','FontSize',6);
        text(0.8,0.55,'pttmean','FontSize',6);
        pvals = [OriStatKO(n).p_value_resp OriStatKO(n).p_value_sel OriStatKO(n).p_value_resp_dF_F OriStatKO(n).p_value_resp_ttest_sep OriStatKO(n).p_value_resp_ttest_mean ];
        xpos=[0 0.2 0.4 0.6 0.8];

        for k=1:5
            p=pvals(k);
            if p < master_alpha(min(k,2))
                colr='r';
            else
                colr='k';
            end
            if p < 0.001
                txt='0';
            else
                txt=sprintf('%6.4f',p);
            end
            text(xpos(k),0.3,txt,'FontSize',6,'Color',colr)
        end
        
        % COLUMN 1 : full trace
        subplot('Position',[col(1) row_y col_width row_height]);
        plot(tcnorm_sliding_bkg2(:,n))
        set(gca,'FontSize',6,'xgrid','on')
        set(gca,'xtick',[0 sz(1)])
        if j~=n_perpage && n~=sz(2)
            set(gca,'XTick',[])
        end
        
        % COLUMN 2 : orientation avg
        subplot('Position',[col(2) row_y col_width row_height]);
        plot(tcavgnorm_sliding_bkg2(:,n))
        hold on
        ax=axis;
        for ngray=1:nstim_per_run
            rectangle('Position',[nframes_per_stim*(ngray-1)+stim_inds(1),ax(3),length(stim_inds),ax(4)-ax(3)],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none')
        end
        plot(tcavgnorm_sliding_bkg2(:,n))
        axis(ax)
        set(gca,'Layer','top','FontSize',6,'xgrid','on')
        set(gca,'xtick',[0 length(tcavgnorm_sliding_bkg2(:,n))])
        if j~=n_perpage && n~=sz(2)
            set(gca,'XTick',[])
        end

        % COLUMN 3 : trial overlays
        subplot('Position',[col(3) row_y col_width row_height]);
        szlong=size(tcnorm_sliding_bkg2);
        szshort=size(tcavgnorm_sliding_bkg2);
        all_shorts_sliding_bkg = reshape(tcnorm_sliding_bkg2(:,n),szshort(1),szlong(1)/szshort(1));
        plot(all_shorts_sliding_bkg)
        hold on
        ax=axis;
        for ngray=1:nstim_per_run
            rectangle('Position',[nframes_per_stim*(ngray-1)+stim_inds(1),ax(3),length(stim_inds),ax(4)-ax(3)],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none')
        end
        plot(all_shorts_sliding_bkg)
        axis(ax)
        set(gca,'Layer','top','FontSize',6,'xgrid','on')
        set(gca,'xtick',[0 szshort(1)])
        if j~=n_perpage && n~=sz(2)
            set(gca,'XTick',[])
        end

        % COLUMN 4 : MEAN STIM ± SEM
        subplot('Position',[col(4) row_y col_width row_height]);
        stim_all = reshape(all_shorts_sliding_bkg,nframes_per_stim,[]);% frames per trial   % N orientations  % N trials        
        mean_response = mean(stim_all,2);
        sem_response  = std(stim_all,0,2) ./ sqrt(size(stim_all,2));
        t = 1:length(mean_response);
        fill([t fliplr(t)],[mean_response'+sem_response' fliplr(mean_response'-sem_response')],[0.8 0.8 0.8],'EdgeColor','none');
        hold on
        %plot(t,mean_response,'k','LineWidth',1)
        ax=axis;
        rectangle('Position',[stim_inds(1),ax(3),length(stim_inds),ax(4)-ax(3)],'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
        fill([t fliplr(t)],[mean_response'+sem_response' fliplr(mean_response'-sem_response')],[0.6 0.6 0.6],'EdgeColor','none');
        plot(t,mean_response,'k','LineWidth',1);
        axis(ax);
        set(gca,'Layer','top','FontSize',6,'xgrid','on')
        set(gca,'xtick',[0 szshort(1)])
        if j~=n_perpage && n~=sz(2)
            set(gca,'XTick',[])
        end
        
        n=n+1;
    end

    % SAVE
    if print_flag_Tc_pages_sliding_bkg==1
        print;
    end
    saveas(h,[analysisDir '\Tc_multiplot' num2str(i) '.fig']);
    saveas(h,[analysisDir '\Tc_multiplot' num2str(i) '.jpg']);
    if print_flag_Tc_pages_sliding_bkg==2
        print;
        print ('-dpsc','-append',[matoutfname,'_tcourse_sliding_bkg_plot.ps']);
    end
end

%plot to look at more of the t-course variants to see what's going on in the analysis
nrows=5; ncols=1; cellperpage=nrows*ncols; 
c1=1; fignum=1; %legstr='';
f(1)=figure('Position',[100 100 600 1000]); hold on; %subplot(nrows,ncols,1); hold on
for ii = 1:sz(2)    
    if c1<=cellperpage && ii<=sz(2)%c1 <=cellperpage*2 && ii<=size(timecourses,2)
        subplot(nrows,ncols,c1);
        yyaxis left
        plot(tcnorm_sliding_bkg2(:,ii));
        yyaxis right
        plot(tc(:,ii)); hold on;
        plot(tc_lowcut(:,ii),'g');
        title(num2str(ii));        
        
        c1=c1+1;                                     
    end
    if c1>cellperpage                        
        c1=1;fignum=fignum+1;        
        f(fignum)=figure('Position',[100 100 600 1000]); subplot(nrows,ncols,1); hold on
    end   
end
for ii = 1:fignum
    f(ii).Renderer='painters'; 
    %saveas(f(ii),[saveDir '\Tc_plot_compare' num2str(ii) '.eps'], 'psc2');%
    saveas(f(ii),[analysisDir '\Tc_plot_compare' num2str(ii) '.fig']);
    saveas(f(ii),[analysisDir '\Tc_plot_compare' num2str(ii) '.jpg'], 'jpeg');
end
%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning on all
if close_fig_flag
    close all
end
%***************************

if breakPoint==9
    return
end

background = adjust_contrastKO(FOV);

outformat = '.jpg';
alpha = master_alpha;
DI_threshold = 0;
response_threshold = 0.00;
showflag=1;

disp('Creating cell maps');
cell_mapsKO_SY(OriStatKO, labelimg, background, analysisDir, outformat, alpha, DI_threshold, response_threshold, showflag);

disp('Plotting cell maps');
doPlotCellMapKO_SY(analysisDir,print_flag_cell_maps, alpha);

disp('Plotting scatter plots');
%dooridistScatterplotCRYC_SY(analysisDir,sFOV,sMatsize);
%doBandwidthScatterPlots(analysisDir)

doV1ScatterPlots(analysisDir,sFOV,sMatsize,print_flag_scatter,alpha);

disp('Count cells');
docountcellsYCKO_SY(analysisDir,alpha(1),alpha(2),DI_threshold);

PValueCellMaps(OriStatKO,labelimg,analysisDir,alpha);
try
    HistFromV1oristat(statFname,pct_for_hist,print_flag_hist)
catch
end

if close_fig_flag
    close all;
end

clear variables;
% load handel;
% sound(y,Fs);
end