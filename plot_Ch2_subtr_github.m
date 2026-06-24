
% This script can be used for red channel subtraction, neuropil subtraction,
%  FRET and 810 nm subtraction. (Actually use ratio rather than subtraction.)

%%% Sample Inputs for Optogenetic Stimulation Runs %%%
ori_run=0; %set to 1 for ori runs, 0 for opto constrict;
dataDir='E:\Analysis\mouse292\aligned-mouse292-020\set4\analyzed\'; %primary data folder
inertDir='E:\Analysis\mouse292\aligned-mouse292-024\set4\analyzed\'; %data to use for correction
saveDir='E:\Analysis\mouse292\aligned-mouse292-020\set4\analyzed\RatioCorr'; %folder to save outputs
pulses=[301 406]; %start and ending frames of optogenetic stimulation


%%% Sample Inputs for Visual Stimulation Runs (Ori8) %%%
ori_run=1; %set to 1 for ori runs, 0 for opto constrict;
dataDir='E:\Analysis\mouse292\aligned-mouse292-016\set4\analyzed\'; %primary data folder
inertDir='E:\Analysis\mouse292\aligned-mouse292-018\set4\analyzed\'; %data to use for correction
saveDir='E:\Analysis\mouse292\aligned-mouse292-016\set4\analyzed\RatioCorr'; %folder to save outputs
numstim=8; %number of stimuli (8 orientations)
BGframes=50; %number of blank frames
stimframes=10; %number of stimulus frames
frames_per_stim=BGframes+stimframes; %number of total frames a for a single stimulus presentation
frameshift=30; %shift for plotting
base_frames=[41:50]; %frames to use for computing baseline (F0)
resp_frames=[51:60]; %frames to use for computing stimulus response (F)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

base_frac=0.5; % The percentage of the lowest frames to be used for the baseline for global scale with ori8
pct_base_opto=25; % Percentile of lowest frames used for baseline in opto runs
pct_violin=2; % Percentage of data (outliers) to exclude for violin plots
global_base=0; % Set to one to normalize the ori8 baseline to the lowest 'base_frac' percent of pixels; 0 will do a sliding baseline over each blank interval
if ~exist(saveDir),    mkdir(saveDir); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotpopavg=1; %set to 1 to do population averaging and stats - need population of cells - single big mask over whole window will crash
avg_all_ori=1;%set to 1 to average across ori's to get a single stimulus response 

set(groot,...
    'defaultAxesFontName','Arial',...
    'defaultTextFontName','Arial',...
    'defaultLegendFontName','Arial',...
    'defaultColorbarFontName','Arial');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ori_run
    if strcmp(dataDir,inertDir) %for neuropil correction data - use the SimpleTc raw data
        load([dataDir 'tc_data.mat'])
        numrep=floor(size(outstruct.tc,1)/frames_per_stim/numstim);
        nTotFrames=numrep*frames_per_stim*numstim;
        tc=outstruct.tc(1:nTotFrames,:);
        inert_tcdata.tc=outstruct.tcS(1:nTotFrames,:);
        
    else
        load([dataDir 'tcourse.mat']); %runThisFunction from ori8 code cuts off tc at exact rep-frame number
        inert_tcdata=load([inertDir '\tcourse.mat']);
        if size(tc,1)~= size(inert_tcdata.tc,1)%for 810 corrrection, it's different runs
            use_length=min([size(tc,1); size(inert_tcdata.tc,1)]);
            tc=tc(1:use_length,:);
            inert_tcdata.tc=inert_tcdata.tc(1:use_length,:);
        end
        
    end
    
    numrep=floor(size(tc,1)/frames_per_stim/numstim);
    nTrials=numrep*numstim;
    tcnormA=zeros(size(tc)); tcnormB=zeros(size(tc)); tcnormC=zeros(size(tc));
    for ii=1:size(tc,2)
        tcC(:,ii)=tc(:,ii)./inert_tcdata.tc(:,ii);
        if global_base    %use the dimmest X% of pixels across the run as the base for all stims
            bkg_win=1:round(size(tc,1)*base_frac);
            tc_lowcut1A= low_cutKO(tc(:,ii), size(tc,1), 'gaussian', 1);%use this here to remove slow trends (fading during run) so lowest pixels aren't biased by that - not sure how well this works
            tc_sortA=sort(tc_lowcut1A);
            bkgA=mean(tc_sortA(bkg_win));
            tcnormA(:,ii)=tc_lowcut1A/bkgA;
            
            tc_lowcut1B= low_cutKO(inert_tcdata.tc(:,ii), size(tc,1), 'gaussian', 1);%use this here to remove slow trends (fading during run) so lowest pixels aren't biased by that - not sure how well this works
            tc_sortB=sort(tc_lowcut1B);
            bkgB=mean(tc_sortB(bkg_winB));
            tcnormB(:,ii)=tc_lowcut1B/bkgB;
            
            tc_lowcut1C= low_cutKO(tcC(:,ii), size(tcC,1), 'gaussian', 1);%use this here to remove slow trends (fading during run) so lowest pixels aren't biased by that - not sure how well this works
            tc_sortC=sort(tc_lowcut1C);
            bkgC=mean(tc_sortC(bkg_winC));
            tcnormC(:,ii)=tc_lowcut1C/bkgC;
                        
        else %
            baselinePointsA = zeros(1, nTrials); baselinePointsB = zeros(1, nTrials); baselinePointsC = zeros(1, nTrials);
            timePoints = zeros(1, nTrials);  % for interpolation (center frame of averaging window)
            
            for jj = 1:nTrials
                trialStart = (jj - 1) * frames_per_stim + 1;
                blankStart = trialStart;
                blankEnd = blankStart + BGframes - 1;
                timePoints(jj) = blankEnd - (length(base_frames)-1)/2;  % center of frames 41–50
                
                baselineWindowA = tc(blankEnd - (length(base_frames)-1) : blankEnd,ii);  % frames 41–50
                baselinePointsA(jj) = mean(baselineWindowA);
                baselineWindowB = inert_tcdata.tc(blankEnd - (length(base_frames)-1) : blankEnd,ii);  % frames 41–50
                baselinePointsB(jj) = mean(baselineWindowB);
                baselineWindowC = tcC(blankEnd - (length(base_frames)-1) : blankEnd,ii);  % frames 41–50
                baselinePointsC(jj) = mean(baselineWindowC);
                
            end
            fullTime = 1:length(tcC(:,ii));
            baselineCurveA = interp1(timePoints, baselinePointsA, fullTime, 'linear', 'extrap');
            baselineA = movquant(baselineCurveA, 0.2, 100,[],'omitnan','truncate');%added this to reduce impact of spike transients during blank window
            tcnormA(:,ii) = tc(:,ii)./baselineA';%use this one - F/F0
            
            baselineCurveB = interp1(timePoints, baselinePointsB, fullTime, 'linear', 'extrap');
            baselineB = movquant(baselineCurveB, 0.2, 100,[],'omitnan','truncate');%added this to reduce impact of spike transients during blank window
            tcnormB(:,ii) = inert_tcdata.tc(:,ii)./baselineB';%use this one - F/F0
            
            baselineCurveC = interp1(timePoints, baselinePointsC, fullTime, 'linear', 'extrap');
            baselineC = movquant(baselineCurveC, 0.2, 100,[],'omitnan','truncate');%added this to reduce impact of spike transients during blank window
            tcnormC(:,ii) = tcC(:,ii)./baselineC';%use this one - F/F0
            
            %figure; hold on; plot(tcC(:,ii)); plot(baselineCurve); plot(timePoints,baselinePoints,'o-'); plot(baseline);
            %tcnormC2(:,ii)=tcC(:,ii)./baseline;
            % Compute normalization - DF/F
            %tcnormCa(:,ii) = (tcC(:,ii)./baselineCurve');
            %tcnormCb(:,ii) = (tcC(:,ii)-baselineCurve');
            %tcnormC2(:,ii) = (tcC(:,ii)-baselineCurve')./baselineCurve';%use this one I think - true DF/F
            
            %figure; hold on; plot(tcnormC(:,ii)); plot(tcnormC2(:,ii));
            %plot(tcnormC3(:,ii)-1); plot(zeros(length(baseline),1),'k--'); plot(timePoints,zeros(length(timePoints),1),'ko');
            
        end
        tcnorm_sm(:,ii)=movmean(tcnormA(:,ii), 3);%these are the normalized uncorrected t-courses
        inert_tcdata.tcnorm_sm(:,ii)=movmean(tcnormB(:,ii), 3);
        tcnormC_sm(:,ii)=movmean(tcnormC(:,ii), 3);
    end
    
    adata = ([0:numstim-1]/ numstim)*2*pi; adata(end+1)=adata(1);
    tcavg=zeros(size(tcnormA,1)/numrep,size(tcnormA,2));
    tcavgC=zeros(size(tcnormA,1)/numrep,size(tcnormA,2));
    tcavg_sm=zeros(size(tcnormA,1)/numrep,size(tcnormA,2));
    tcavgC_sm=zeros(size(tcnormA,1)/numrep,size(tcnormA,2));
    x1=1:size(tcavg,1);
    x_orig=[x1 x1(end:-1:1) x1(1)];
    for ii=1:size(tc,2)
        reshaped_data1 = reshape(tcnormA(:,ii), [], numrep);
        reshaped_data2 = reshape(tcnormB(:,ii), [], numrep);
        reshaped_data3 = reshape(tcnormC(:,ii), [], numrep);
        reshaped_data1_shift = circshift(reshaped_data1, -frameshift,1);
        reshaped_data2_shift = circshift(reshaped_data2, -frameshift,1);
        reshaped_data3_shift = circshift(reshaped_data3, -frameshift,1);
        tcavg(:,ii)=mean(reshaped_data1_shift,2);
        inert_tcdata.tcavg(:,ii)=mean(reshaped_data2_shift,2);
        tcavgC(:,ii)=mean(reshaped_data3_shift,2);
        % for smoothed %%%%%%%%%%%%%
        reshaped_data1sm = reshape(tcnorm_sm(:,ii), [], numrep);
        reshaped_data2sm = reshape(inert_tcdata.tcnorm_sm(:,ii), [], numrep);
        reshaped_data3sm = reshape(tcnormC_sm(:,ii), [], numrep);
        reshaped_data1sm_shift = circshift(reshaped_data1sm, -frameshift,1);
        reshaped_data2sm_shift = circshift(reshaped_data2sm, -frameshift,1);
        reshaped_data3sm_shift = circshift(reshaped_data3sm, -frameshift,1);
        tcavg_sm(:,ii)=mean(reshaped_data1sm_shift,2);
        inert_tcdata.tcavg_sm(:,ii)=mean(reshaped_data2sm_shift,2);
        tcavgC_sm(:,ii)=mean(reshaped_data3sm_shift,2);
        
        
        cellsem1=std(reshaped_data1_shift,[],2,'omitnan')./sqrt(numrep);
        y1A=tcavg(:,ii)-cellsem1;
        y2A=tcavg(:,ii)+cellsem1;
        y_origA(:,ii)=[y1A' y2A(end:-1:1)' y1A(1)];
        cellsem2=std(reshaped_data2_shift,[],2,'omitnan')./sqrt(numrep);
        y1B=inert_tcdata.tcavg(:,ii)-cellsem2;
        y2B=inert_tcdata.tcavg(:,ii)+cellsem2;
        y_origB(:,ii)=[y1B' y2B(end:-1:1)' y1B(1)];
        cellsem3=std(reshaped_data3_shift,[],2,'omitnan')./sqrt(numrep);
        y1C=tcavgC(:,ii)-cellsem3;
        y2C=tcavgC(:,ii)+cellsem3;
        y_origC(:,ii)=[y1C' y2C(end:-1:1)' y1C(1)];
    end
    
    % Plotting
    numPerFigure = 10;  % 10 rows per figure
    numFigures = ceil(size(tc,2) / numPerFigure);
    
    for ii=1:numFigures
        f1(ii)=figure('Position',[100 100 600 1000]); hold on; 
        ind1=1;        
        for jj = 1:numPerFigure
            seriesIdx = (ii - 1) * numPerFigure + jj;
            if seriesIdx > size(tc,2)
                break;
            end
            row = jj;
            % Left column: raw
            subplot(numPerFigure, 2, 2*row - 1); hold on;
            ymax=max([tc(:,seriesIdx); inert_tcdata.tc(:,seriesIdx)]);
            ymin=min([tc(:,seriesIdx); inert_tcdata.tc(:,seriesIdx)]);
            plot(tc(:,seriesIdx),'Color',[0.3 0.7 0.3]);
            plot(inert_tcdata.tc(:,seriesIdx),'r');
            title(sprintf('Raw #%d', seriesIdx));
            if row < numPerFigure
                set(gca, 'XTickLabel', []);
            end
            
            % Right column: normalized averaged across trials
            subplot(numPerFigure, 2, 2*row); hold on;
            ymin1=floor(min([y_origA(:,seriesIdx) y_origB(:,seriesIdx) y_origC(:,seriesIdx)])*20)/20; ymax1=ceil(max([y_origA(:,seriesIdx) y_origB(:,seriesIdx) y_origC(:,seriesIdx)])*20)/20;
            ymin2=floor(min([tcavg(:,seriesIdx) inert_tcdata.tcavg(:,seriesIdx) tcavgC(:,seriesIdx)])*20)/20; ymax2=ceil(max([tcavg(:,seriesIdx) inert_tcdata.tcavg(:,seriesIdx) tcavgC(:,seriesIdx)])*20)/20;
            ymin=min([ymin1 ymin2]); ymax=max([ymax1 ymax2]);
            for n=1:numstim
                rectangle('Position',[(n-1)*frames_per_stim+BGframes-frameshift,ymin,stimframes,ymax-ymin],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none');%[0.7 0.7 0.7]
            end
            fill(x_orig,y_origA(:,seriesIdx),[0.6 1 0.6],'EdgeColor','none');
            fill(x_orig,y_origB(:,seriesIdx),[1 0.6 0.6],'EdgeColor','none');
            fill(x_orig,y_origC(:,seriesIdx),[0.6 0.6 1],'EdgeColor','none');
            plot(tcavg(:,seriesIdx),'Color',[0.3 0.7 0.3]);
            plot(inert_tcdata.tcavg(:,seriesIdx),'r');
            plot(tcavgC(:,seriesIdx),'b');
            title(sprintf('Normalized trial average #%d', seriesIdx));
            if row < numPerFigure
                set(gca, 'XTickLabel', []);
            end
        end
    end
    
    for jj=1:size(f1,2)
        %print(f1(jj),'-painters','-depsc2',[saveDir '\indiv_cells' num2str(jj) '.eps'])%
        savefig(f1(jj),[saveDir '\indiv_cells' num2str(jj) '.fig']);
        exportgraphics(f1(jj),[saveDir '\indiv_cells' num2str(jj) '.jpg']);
    end        
    %close all;

    %Pop Avg
    if ~exist('excludecells','var') 
        excludecells = [];
    end    
    good_cells = setdiff(1:size(tcavgC,2), excludecells);

    cellavgC = mean(tcavgC(:,good_cells),2); %
    cellavg_func=mean(tcavg(:,good_cells),2);
    cellavg_inert=mean(inert_tcdata.tcavg(:,good_cells),2);
    
    cellsemC=std(tcavgC(:,good_cells),0,2,'omitnan')./sqrt(size(tcavgC,2));
    y1C=cellavgC-cellsemC;    y2C=cellavgC+cellsemC;    yC=[y1C' y2C(end:-1:1)' y1C(1)];
    cellsemF=std(tcavg(:,good_cells),0,2,'omitnan')./sqrt(size(tcavgC,2));
    y1F=cellavg_func-cellsemF;    y2F=cellavg_func+cellsemF;    yF=[y1F' y2F(end:-1:1)' y1F(1)];
    cellsemI=std(inert_tcdata.tcavg(:,good_cells),0,2,'omitnan')./sqrt(size(tcavgC,2));
    y1I=cellavg_inert-cellsemI;    y2I=cellavg_inert+cellsemI;    yI=[y1I' y2I(end:-1:1)' y1I(1)];
    x1=1:size(tcavgC,1);    x=[x1 x1(end:-1:1) x1(1)];
    ymin=floor(min([y1C y1F y1I],[],'all')*20)/20; ymax=ceil(max([y2C y2F y2I],[],'all')*20)/20;    
    
    f4=figure; hold on    
    for n=1:numstim
        rectangle('Position',[(n-1)*frames_per_stim+BGframes-frameshift,ymin,stimframes,ymax-ymin],'FaceColor',[0.7 0.7 0.7]);
    end
    
    fill(x,yC,[0.6 0.6 1],'EdgeColor','none');
    h1=plot(cellavgC,'b','LineWidth',2); hold on;
    fill(x,yF,[0.6 1 0.6],'EdgeColor','none');
    h2=plot(cellavg_func,'Color',[0.3 0.7 0.3],'LineWidth',2); hold on;
    fill(x,yI,[1 0.6 0.6],'EdgeColor','none');
    h3=plot(cellavg_inert,'r','LineWidth',2); hold on;
    x_limits=xlim;
    x_pos=(x_limits(2)-x_limits(1))/1.7;
    ylim([ymin ymax])
    titlestr=['Pop avg'];
    title(titlestr);
    legend([h1 h2 h3],{'Corr','Func','Inert'});
    print(f4,'-painters','-depsc2',[saveDir '\' titlestr '.eps'])%
    %saveas(f4,[saveDir '\' titlestr '.eps'], 'psc2');
    saveas(f4,[saveDir '\' titlestr '.fig']);
    saveas(f4,[saveDir '\' titlestr '.tif']);
    
    if avg_all_ori
        f4b=figure; hold on  
        stim_allC = reshape(cellavgC,frames_per_stim,[]);% frames per trial   % N orientations  
        mean_responseC = mean(stim_allC,2);
        sem_responseC  = std(stim_allC,0,2) ./ sqrt(size(stim_allC,2));
        
        stim_allF = reshape(cellavg_func,frames_per_stim,[]);% frames per trial   % N orientations  
        mean_responseF = mean(stim_allF,2);
        sem_responseF  = std(stim_allF,0,2) ./ sqrt(size(stim_allF,2));
        
        stim_allI = reshape(cellavg_inert,frames_per_stim,[]);% frames per trial   % N orientations  
        mean_responseI = mean(stim_allI,2);
        sem_responseI  = std(stim_allI,0,2) ./ sqrt(size(stim_allI,2));        
        
        t = 1:length(mean_responseC);
        fill([t fliplr(t)],[mean_responseC'+sem_responseC' fliplr(mean_responseC'-sem_responseC')],[0.6 0.6 1],'EdgeColor','none');
        fill([t fliplr(t)],[mean_responseF'+sem_responseF' fliplr(mean_responseF'-sem_responseF')],[0.6 1 0.6],'EdgeColor','none');
        fill([t fliplr(t)],[mean_responseI'+sem_responseI' fliplr(mean_responseI'-sem_responseI')],[1 0.6 0.6],'EdgeColor','none');        
        hold on
        ax=axis;
        rectangle('Position',[BGframes-frameshift+1,ax(3),stimframes,ax(4)-ax(3)],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none');
        fill([t fliplr(t)],[mean_responseC'+sem_responseC' fliplr(mean_responseC'-sem_responseC')],[0.6 0.6 1],'EdgeColor','none');
        fill([t fliplr(t)],[mean_responseF'+sem_responseF' fliplr(mean_responseF'-sem_responseF')],[0.6 1 0.6],'EdgeColor','none');
        fill([t fliplr(t)],[mean_responseI'+sem_responseI' fliplr(mean_responseI'-sem_responseI')],[1 0.6 0.6],'EdgeColor','none');                
        h1=plot(t,mean_responseC,'b','LineWidth',2);
        h2=plot(t,mean_responseF,'Color',[0.3 0.7 0.3],'LineWidth',2);
        h3=plot(t,mean_responseI,'r','LineWidth',2);   
        legend([h1 h2 h3],{'Corr','Func','Inert'});
        titlestr2=['Pop avg stim avg'];
        title(titlestr2);
        print(f4b,'-painters','-depsc2',[saveDir '\' titlestr2 '.eps'])%
        saveas(f4b,[saveDir '\' titlestr2 '.fig']);
        saveas(f4b,[saveDir '\' titlestr2 '.tif']);
        pop_params.mean_responseC=mean_responseC; pop_params.mean_responseF=mean_responseF; pop_params.mean_responseI=mean_responseI;
        pop_params.sem_responseC=sem_responseC; pop_params.sem_responseF=sem_responseF; pop_params.sem_responseI=sem_responseI;        
    end    
    
    if plotpopavg
        %Compute DF/F for each ori across population
        [nFrames,nCells] = size(tcnormA);
        nStimuli = 8;    framesPerStim = 60;
        dffTable_orig = zeros(nCells, numstim);
        dffTable_corr = zeros(nCells, numstim);
        for ii=1:size(tcnormA,2)
            for stimIdx = 1:numstim
                stimStart = (stimIdx - 1) * frames_per_stim + 1;
                baselineFrames = stimStart + (base_frames) - 1;  %
                responseFrames = stimStart + (resp_frames) - 1;  %
                % Compute mean baseline and response per cell
                F_base_orig = mean(tcnormA(baselineFrames,ii));     %
                F_resp_orig = mean(tcnormA(responseFrames,ii));     %
                F_base_corr = mean(tcnormC(baselineFrames,ii));     %
                F_resp_corr = mean(tcnormC(responseFrames,ii));     %
                % Compute DF/F
                dff_orig = (F_resp_orig - F_base_orig) ./ F_base_orig;
                dffTable_orig(ii, stimIdx) = dff_orig;
                dff_corr = (F_resp_corr - F_base_corr) ./ F_base_corr;
                dffTable_corr(ii, stimIdx) = dff_corr;
            end            
        end
        mean_resp_orig=mean(dffTable_orig,2);        
        mean_resp_corr=mean(dffTable_corr,2);
        max_resp_orig=max(dffTable_orig,[],2);
        max_resp_corr=max(dffTable_corr,[],2);

        orig_mean_mean=mean(mean_resp_orig(good_cells));    corr_mean_mean=mean(mean_resp_corr(good_cells));
        orig_mean_max=mean(max_resp_orig(good_cells));    corr_mean_max=mean(max_resp_corr(good_cells));
        orig_std_mean=std(mean_resp_orig(good_cells));    corr_Std_mean=std(mean_resp_corr(good_cells));
        orig_std_max=std(max_resp_orig(good_cells));    corr_Std_max=std(max_resp_corr(good_cells));
        means1 = [orig_mean_mean, corr_mean_mean];    stds1 = [orig_std_mean, corr_Std_mean];
        means2 = [orig_mean_max, corr_mean_max];    stds2 = [orig_std_max, corr_Std_max];
        %
        f5=figure; hold on;
        swarmchart(zeros(length(good_cells),1)+1,mean_resp_orig(good_cells),[],'g');%,'XJitterWidth',0.2);
        swarmchart(zeros(length(good_cells),1)+2,mean_resp_corr(good_cells),[],'b');%,'XJitterWidth',0.2);
        data = {mean_resp_orig(good_cells), mean_resp_corr(good_cells)};
        violin(data);
        xticks([1 2]);    xticklabels({'Orig', 'Corrected'});
        ylabel('Value');    title('Mean ± Std Dev of mean DF/F');
        
        clear data_trimmed
        f5b=figure; hold on;
        swarmchart(zeros(length(good_cells),1)+1,mean_resp_orig(good_cells),[],'g');%,'XJitterWidth',0.2);
        swarmchart(zeros(length(good_cells),1)+2,mean_resp_corr(good_cells),[],'b');%,'XJitterWidth',0.2);
        data = {mean_resp_orig(good_cells), mean_resp_corr(good_cells)};
        for kk=1:length(data)
            prc = prctile(data{kk}, [pct_violin 100-pct_violin]);
            data_trimmed{1,kk} = data{kk}(data{kk} >= prc(1) & data{kk} <= prc(2));
        end
        violin(data_trimmed);
        
        xticks([1 2]);    xticklabels({'Orig', 'Corrected'});
        ylabel('Value');    title('Mean ± Std Dev of mean DF/F trimmed');
        %saveas(f5,[saveDir '\Pop distribution meanDF.eps'], 'psc2');
        print(f5,'-vector','-depsc2',[saveDir '\Pop distribution meanDF.eps'])%
        saveas(f5,[saveDir '\Pop distribution meanDF.fig']);
        saveas(f5,[saveDir '\Pop distribution meanDF.tif']);
        print(f5b,'-vector','-depsc2',[saveDir '\Pop distribution trimmed meanDF.eps'])%
        saveas(f5b,[saveDir '\Pop distribution trimmed meanDF.fig']);
        saveas(f5b,[saveDir '\Pop distribution trimmed meanDF.tif']);
        
        [h1a, p1a] = lillietest(mean_resp_orig(good_cells));    [h2a, p2a] = lillietest(mean_resp_corr(good_cells));
        [h_ttest_a, p_ttest_a, ci, stats] = ttest(mean_resp_orig(good_cells), mean_resp_corr(good_cells));    p_non_param_a = signrank(mean_resp_orig(good_cells), mean_resp_corr(good_cells));
        
        f6=figure; hold on;
        swarmchart(zeros(length(good_cells),1)+1,max_resp_orig(good_cells),[],'g');%,'XJitterWidth',0.2);
        swarmchart(zeros(length(good_cells),1)+2,max_resp_corr(good_cells),[],'b');%,'XJitterWidth',0.2);
        data = {max_resp_orig(good_cells), max_resp_corr(good_cells)};
        violin(data);
        xticks([1 2]);    xticklabels({'Orig', 'Corrected'});
        ylabel('Value');    title('Mean ± Std Dev of max DF/F');
        
        f6b=figure; hold on;
        swarmchart(zeros(length(good_cells),1)+1,max_resp_orig(good_cells),[],'g');%,'XJitterWidth',0.2);
        swarmchart(zeros(length(good_cells),1)+2,max_resp_corr(good_cells),[],'b');%,'XJitterWidth',0.2);
        data = {max_resp_orig(good_cells), max_resp_corr(good_cells)};
        for kk=1:length(data)
            prc = prctile(data{kk}, [pct_violin 100-pct_violin]);
            data_trimmed{1,kk} = data{kk}(data{kk} >= prc(1) & data{kk} <= prc(2));
        end
        violin(data_trimmed);
        
        xticks([1 2]);    xticklabels({'Orig', 'Corrected'});
        ylabel('Value');    title('Mean ± Std Dev of max DF/F trimmed');
        %saveas(f6,[saveDir '\Pop distribution maxDF.eps'], 'psc2');
        print(f6,'-painters','-depsc2',[saveDir '\Pop distribution maxDF.eps'])%
        saveas(f6,[saveDir '\Pop distribution maxDF.fig']);
        saveas(f6,[saveDir '\Pop distribution maxDF.tif']);
        print(f6b,'-painters','-depsc2',[saveDir '\Pop distribution trimmed maxDF.eps'])%
        saveas(f6b,[saveDir '\Pop distribution trimmed maxDF.fig']);
        saveas(f6b,[saveDir '\Pop distribution trimmed maxDF.tif']);
        
        [h1b, p1b] = lillietest(max_resp_orig(good_cells));    [h2b, p2b] = lillietest(max_resp_corr(good_cells));
        [h_ttest_b, p_ttest_b, ci, stats] = ttest(max_resp_orig(good_cells), max_resp_corr(good_cells));    p_non_param_b = signrank(max_resp_orig(good_cells), max_resp_corr(good_cells));
        pop_params.mean_resp_orig=mean_resp_orig; pop_params.mean_resp_corr=mean_resp_corr; pop_params.max_resp_orig=max_resp_orig; pop_params.max_resp_corr=max_resp_corr;
        pop_params.orig_mean_p1=p1a; pop_params.corr_mean_p=p2a; pop_params.orig_max_p1=p1b; pop_params.corr_max_p=p2b;
        pop_params.ttest_mean_p=p_ttest_a; pop_params.non_param_mean_p=p_non_param_a; pop_params.ttest_max_p=p_ttest_b; pop_params.non_param_max_p=p_non_param_b;
        
    end
    run_params.numstim=numstim; run_params.BGframes=BGframes; run_params.stimframes=stimframes; run_params.frames_per_stim=frames_per_stim; run_params.frameshift=frameshift;%
    run_params.base_frames=base_frames; run_params.resp_frames=resp_frames;
    pop_params.cellavgC = cellavgC;     pop_params.cellavg_func=cellavg_func;    pop_params.cellavg_inert=cellavg_inert; % some of these may overlap the parameters above if plotpopavg is run
    if plotpopavg
        save([saveDir '\tcourse.mat'],'tcC','tcnormC','tcavgC','tcavg','tcnormA','tcnormB','inert_tcdata','run_params','pop_params','dffTable_orig','dffTable_corr');
    else
        save([saveDir '\tcourse.mat'],'tcC','tcnormC','tcavgC','tcavg','tcnormA','tcnormB','inert_tcdata','run_params','pop_params');
    end
else %opto constrict run
    load([dataDir 'tc_data.mat'])
    if strcmp(dataDir,inertDir) %for neuropil correction data
        use_data=outstruct.tc;
        use_inert_data=outstruct.tcS;
    else
        inert_tcdata=load([inertDir '\tc_data.mat']);
        use_data=timecoursesraw;
        use_inert_data=inert_tcdata.timecoursesraw;
        if size(use_data,1)~= size(use_inert_data,1)%for 810 corrrection, it's different runs
            use_length=min([size(use_data,1); size(use_inert_data,1)]);
            use_data=timecoursesraw(1:use_length,:);
            use_inert_data=inert_tcdata.timecoursesraw(1:use_length,:);
        end
    end
    for ii=1:size(use_data,2)
        
        tcC(:,ii)=use_data(:,ii)./use_inert_data(:,ii);
        %%% shouldn't need de-trending with the ratio signal
        %tc_lowcut1= low_cutKO(tcC(:,ii), size(tcC,1), 'gaussian', 1);%use this here to remove slow trends (fading during run) so lowest pixels aren't biased by that - not sure how well this works
        %tcnormC_sm(:,ii)=tcnorm_sm(:,ii)./inert_tcdata.tcnorm_sm(:,ii);
        %%% try median of final 75 frames of blank %%%
        %baseline =median(tcC(pulses(1)-base_frames-1:pulses(1)-1,ii));
        %%% try lowest pct_base_opto (25%) of all blank window pixels after ratio and smoothing %%%
        tcCsm(:,ii)=movmean(tcC(:,ii),10);
        baseline=prctile(tcCsm(1:pulses(1)-1,ii),pct_base_opto);
        tcnormC(:,ii) = (tcC(:,ii)-baseline)./baseline;% DF/F
        tcnormC_sm(:,ii)=movmean(tcnormC(:,ii), 3);
        
        tcAsm(:,ii)=movmean(use_data(:,ii),10);
        baselineA=prctile(tcAsm(1:pulses(1)-1,ii),pct_base_opto);
        tcnormA(:,ii) = (use_data(:,ii)-baselineA)./baselineA;% DF/F
        tcnorm_sm(:,ii)=movmean(tcnormA(:,ii), 3);
        
        tcBsm(:,ii)=movmean(use_inert_data(:,ii),10);
        baselineB=prctile(tcBsm(1:pulses(1)-1,ii),pct_base_opto);
        tcnormB(:,ii) = (use_inert_data(:,ii)-baselineB)./baselineB;% DF/F
        inert_tcdata.tcnorm_sm(:,ii)=movmean(tcnormB(:,ii), 3);
        
    end
    %figure; hold on; plot(timecoursesUnCorr(:,ii)); plot(timecourses(:,ii)); plot(timecoursesraw(:,ii));
    %figure; hold on; plot(tcC(:,ii)); plot(tcnormC(:,ii)); plot(tcCsm(:,ii));
    %plot(zeros(length(tcC(:,ii)),1),'k--');
    
    % Plotting
    numPerFigure = 10;  % 10 rows per figure
    numFigures = ceil(size(use_data,2) / numPerFigure);
    
    for ii=1:numFigures
        f1(ii)=figure('Position',[100 100 600 1000]); hold on; %subplot(nrows,ncols,1); hold on;
        ind1=1;
        
        for jj = 1:numPerFigure
            seriesIdx = (ii - 1) * numPerFigure + jj;
            if seriesIdx > size(use_data,2)
                break;
            end
            row = jj;
            % Left column: raw
            subplot(numPerFigure, 2, 2*row - 1); hold on;
            ymax=max([use_data(:,seriesIdx); use_inert_data(:,seriesIdx)]);
            ymin=min([use_data(:,seriesIdx); use_inert_data(:,seriesIdx)]);
            rectangle('Position',[pulses(1),ymin,pulses(2)-pulses(1),ymax-ymin],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none');%[0.7 0.7 0.7]
            plot(use_data(:,seriesIdx),'Color',[0.3 0.7 0.3]);
            plot(use_inert_data(:,seriesIdx),'r');
            
            title(sprintf('Raw #%d', seriesIdx));
            if row < numPerFigure
                set(gca, 'XTickLabel', []);
            end
            
            % Right column: smoothed
            subplot(numPerFigure, 2, 2*row); hold on;
            ymax=max([tcnorm_sm(:,seriesIdx); inert_tcdata.tcnorm_sm(:,seriesIdx); tcnormC_sm(:,seriesIdx)]);
            ymin=min([tcnorm_sm(:,seriesIdx); inert_tcdata.tcnorm_sm(:,seriesIdx); tcnormC_sm(:,seriesIdx)]);
            rectangle('Position',[pulses(1),ymin,pulses(2)-pulses(1),ymax-ymin],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none');%[0.7 0.7 0.7]
            plot(tcnorm_sm(:,seriesIdx),'Color',[0.3 0.7 0.3]);
            plot(inert_tcdata.tcnorm_sm(:,seriesIdx),'r');
            plot(tcnormC_sm(:,seriesIdx),'b');
            title(sprintf('Smooth Normalized #%d', seriesIdx));
            if row < numPerFigure
                set(gca, 'XTickLabel', []);
            end
        end
    end
    
    
    for jj=1:size(f1,2)
        %print(f1(jj),'-painters','-depsc2',[saveDir '\avgTcourse' num2str(jj) '.eps'])%
        savefig(f1(jj),[saveDir '\avgTcourse' num2str(jj) '.fig']);
        exportgraphics(f1(jj),[saveDir '\avgTcourse' num2str(jj) '.jpg']);
    end
    
    %Pop Avg
    if ~exist('excludecells','var') 
        excludecells = [];
    end    
    good_cells = setdiff(1:size(tcnormC,2), excludecells);

    cellavgC = mean(tcnormC(:,good_cells),2); %
    cellavg_func=mean(tcnormA(:,good_cells),2);
    cellavg_inert=mean(tcnormB(:,good_cells),2);
    
    cellsemC=std(tcnormC(:,good_cells),0,2,'omitnan')./sqrt(length(good_cells));
    y1C=cellavgC-cellsemC;    y2C=cellavgC+cellsemC;    yC=[y1C' y2C(end:-1:1)' y1C(1)];
    cellsemF=std(use_data(:,good_cells),0,2,'omitnan')./sqrt(length(good_cells));
    y1F=cellavg_func-cellsemF;    y2F=cellavg_func+cellsemF;    yF=[y1F' y2F(end:-1:1)' y1F(1)];
    cellsemI=std(use_inert_data(:,good_cells),0,2,'omitnan')./sqrt(length(good_cells));
    y1I=cellavg_inert-cellsemI;    y2I=cellavg_inert+cellsemI;    yI=[y1I' y2I(end:-1:1)' y1I(1)];
    x1=1:size(tcnormC,1);    x=[x1 x1(end:-1:1) x1(1)];
    ymin=floor(min([y1C y1F y1I],[],'all')*20)/20; ymax=ceil(max([y2C y2F y2I],[],'all')*20)/20;    %ymin=.85; ymax=1.4;
    f4=figure; hold on    %a1=axes('Parent',f1); hold(a1,'on');
    rectangle('Position',[pulses(1),ymin,pulses(2)-pulses(1),ymax-ymin],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none');%[0.7 0.7 0.7]
    fill(x,yC,[0.6 0.6 1],'EdgeColor','none');
    h1=plot(cellavgC,'b','LineWidth',2); hold on;
    fill(x,yF,[0.6 1 0.6],'EdgeColor','none');
    h2=plot(cellavg_func,'Color',[0.3 0.7 0.3],'LineWidth',2); hold on;
    fill(x,yI,[1 0.6 0.6],'EdgeColor','none');
    h3=plot(cellavg_inert,'r','LineWidth',2); hold on;
    x_limits=xlim;
    x_pos=(x_limits(2)-x_limits(1))/1.7;
   ylim([ymin ymax])
    titlestr=['Pop avg'];
    title(titlestr)
    legend([h1 h2 h3],{'Corr','Func','Inert'});
    %saveas(f4,[saveDir '\' titlestr '.eps'], 'psc2');
    print(f4,'-painters','-depsc2',[saveDir '\' titlestr '.eps'])%
    saveas(f4,[saveDir '\' titlestr '.fig']);
    saveas(f4,[saveDir '\' titlestr '.tif']);
    
    if plotpopavg
        for ii=1:size(tcC,2)
            orig_dur_amp(ii,1)=mean(tcnormA(pulses(1):pulses(2),ii));
            corr_dur_amp(ii,1)=mean(tcnormC(pulses(1):pulses(2),ii));
            orig_bef_amp(ii,1)=mean(tcnormA(1:pulses(1)-1,ii));
            corr_bef_amp(ii,1)=mean(tcnormC(1:pulses(1)-1,ii));
            orig_aft_amp(ii,1)=mean(tcnormA(pulses(2)+1:end,ii));
            corr_aft_amp(ii,1)=mean(tcnormC(pulses(2)+1:end,ii));
        end
        orig_bef_mean=mean(orig_bef_amp(good_cells));
        corr_bef_mean=mean(corr_bef_amp(good_cells));
        orig_bef_std=std(orig_bef_amp(good_cells));
        corr_bef_Std=std(corr_bef_amp(good_cells));
        means_bef = [orig_bef_mean, corr_bef_mean];
        stds_bef = [orig_bef_std, corr_bef_Std];
        [h1_bef, p1_bef] = lillietest(orig_bef_amp(good_cells));    [h2_bef, p2_bef] = lillietest(corr_bef_amp(good_cells));
        [h_ttest_bef, p_ttest_bef, ci_bef, stats_bef] = ttest(orig_bef_amp(good_cells), corr_bef_amp(good_cells));    p_non_param_bef = signrank(orig_bef_amp(good_cells), corr_bef_amp(good_cells));
        
        orig_dur_mean=mean(orig_dur_amp(good_cells));
        corr_dur_mean=mean(corr_dur_amp(good_cells));
        orig_dur_std=std(orig_dur_amp(good_cells));
        corr_dur_Std=std(corr_dur_amp(good_cells));
        means_dur = [orig_dur_mean, corr_dur_mean];
        stds_dur = [orig_dur_std, corr_dur_Std];
        [h1_dur, p1_dur] = lillietest(orig_dur_amp(good_cells));    [h2_dur, p2_dur] = lillietest(corr_dur_amp(good_cells));
        [h_ttest_dur, p_ttest_dur, ci_dur, stats_dur] = ttest(orig_dur_amp(good_cells), corr_dur_amp(good_cells));    p_non_param_dur = signrank(orig_dur_amp(good_cells), corr_dur_amp(good_cells));
        
        orig_aft_mean=mean(orig_aft_amp(good_cells));
        corr_aft_mean=mean(corr_aft_amp(good_cells));
        orig_aft_std=std(orig_aft_amp(good_cells));
        corr_aft_Std=std(corr_aft_amp(good_cells));
        means_aft = [orig_aft_mean, corr_aft_mean];
        stds_aft = [orig_aft_std, corr_aft_Std];
        [h1_aft, p1_aft] = lillietest(orig_aft_amp(good_cells));    [h2_aft, p2_aft] = lillietest(corr_aft_amp(good_cells));
        [h_ttest_aft, p_ttest_aft, ci_aft, stats_aft] = ttest(orig_aft_amp(good_cells), corr_aft_amp(good_cells));    p_non_param_aft = signrank(orig_aft_amp(good_cells), corr_aft_amp(good_cells));
        
        f5=figure; hold on;
        swarmchart(zeros(size(tcC,2),1)+1,orig_bef_amp(good_cells),[],'g');%,'XJitterWidth',0.2);
        swarmchart(zeros(size(tcC,2),1)+2,corr_bef_amp(good_cells),[],'b');%,'XJitterWidth',0.2);
        swarmchart(zeros(size(tcC,2),1)+3,orig_dur_amp(good_cells),[],'g');%,'XJitterWidth',0.2);
        swarmchart(zeros(size(tcC,2),1)+4,corr_dur_amp(good_cells),[],'b');%,'XJitterWidth',0.2);
        swarmchart(zeros(size(tcC,2),1)+5,orig_aft_amp(good_cells),[],'g');%,'XJitterWidth',0.2);
        swarmchart(zeros(size(tcC,2),1)+6,corr_aft_amp(good_cells),[],'b');%,'XJitterWidth',0.2);
        data = {orig_bef_amp(good_cells), corr_bef_amp(good_cells), orig_dur_amp(good_cells), corr_dur_amp(good_cells), orig_aft_amp(good_cells), corr_aft_amp(good_cells)};
        violin(data);
        xticks([1.5 3.5 5.5]);    xticklabels({'Bef', 'Dur', 'Aft'});
        legend({'Orig'; 'Corrected'});  ylabel('Value');    title('population Mean ± Std Dev');
        
        clear data_trimmed
        f5b=figure; hold on;
        swarmchart(zeros(size(tcC,2),1)+1,orig_bef_amp(good_cells),[],'g');%,'XJitterWidth',0.2);
        swarmchart(zeros(size(tcC,2),1)+2,corr_bef_amp(good_cells),[],'b');%,'XJitterWidth',0.2);
        swarmchart(zeros(size(tcC,2),1)+3,orig_dur_amp(good_cells),[],'g');%,'XJitterWidth',0.2);
        swarmchart(zeros(size(tcC,2),1)+4,corr_dur_amp(good_cells),[],'b');%,'XJitterWidth',0.2);
        swarmchart(zeros(size(tcC,2),1)+5,orig_aft_amp(good_cells),[],'g');%,'XJitterWidth',0.2);
        swarmchart(zeros(size(tcC,2),1)+6,corr_aft_amp(good_cells),[],'b');%,'XJitterWidth',0.2);
        data = {orig_bef_amp(good_cells), corr_bef_amp(good_cells), orig_dur_amp(good_cells), corr_dur_amp(good_cells), orig_aft_amp(good_cells), corr_aft_amp(good_cells)};
        for kk=1:length(data)
            prc = prctile(data{kk}, [pct_violin 100-pct_violin]);
            data_trimmed{1,kk} = data{kk}(data{kk} >= prc(1) & data{kk} <= prc(2));
        end
        violin(data_trimmed);
        xticks([1.5 3.5 5.5]);    xticklabels({'Bef', 'Dur', 'Aft'});
        legend({'Orig'; 'Corrected'});  ylabel('Value');    title('trimmed population Mean ± Std Dev');
        
        %saveas(f5,[saveDir '\Pop distribution .eps'], 'psc2');
        print(f5,'-painters','-depsc2',[saveDir '\Pop distribution .eps'])%
        saveas(f5,[saveDir '\Pop distribution .fig']);
        saveas(f5,[saveDir '\Pop distribution .tif']);
        print(f5b,'-painters','-depsc2',[saveDir '\Pop distribution trimmed.eps'])%
        saveas(f5b,[saveDir '\Pop distribution trimmed.fig']);
        saveas(f5b,[saveDir '\Pop distribution trimmed.tif']);

        pop_params.orig_dur_amp=orig_dur_amp; pop_params.corr_dur_amp=corr_dur_amp; pop_params.orig_bef_amp=orig_bef_amp; 
        pop_params.corr_bef_amp=corr_bef_amp; pop_params.orig_aft_amp=orig_aft_amp; pop_params.corr_aft_amp=corr_aft_amp;
        
        pop_params.orig_bef_mean=orig_bef_mean; pop_params.corr_bef_mean=corr_bef_mean; pop_params.orig_bef_norm_p=p1_bef; pop_params.corr_bef_norm_p=p2_bef;
        pop_params.ttest_p_bef=p_ttest_bef; pop_params.non_param_p_bef=p_non_param_bef;
        
        pop_params.orig_dur_mean=orig_dur_mean; pop_params.corr_dur_mean=corr_dur_mean; pop_params.orig_dur_norm_p=p1_dur; pop_params.corr_dur_norm_p=p2_dur;
        pop_params.ttest_p_dur=p_ttest_dur; pop_params.non_param_p_dur=p_non_param_dur;
        
        pop_params.orig_aft_mean=orig_aft_mean; pop_params.corr_aft_mean=corr_aft_mean; pop_params.orig_aft_norm_p=p1_aft; pop_params.corr_aft_norm_p=p2_aft;
        pop_params.ttest_p_aft=p_ttest_aft; pop_params.non_param_p_aft=p_non_param_aft;
    end
    run_params.pulses=pulses; run_params.pct_base_opto=pct_base_opto;
    pop_params.cellavgC=cellavgC; pop_params.cellavg_func=cellavg_func; pop_params.cellavg_inert=cellavg_inert; 
    pop_params.cellsemC=cellsemC; pop_params.cellsemF=cellsemF; pop_params.cellsemI=cellsemI; 
    save([saveDir '\tcourse.mat'],'tcC','tcnormC','tcnormA','tcnormB','inert_tcdata','run_params','pop_params');
        
end



