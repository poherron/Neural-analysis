

% To plot the population summary across animals for the different corrections for DVAA. This is to plot the corrections
% for the prolonged opto data. Original data was generated in plot_Ch2_subtr.m 

saveDir='E:\Analysis\population_analysis\DVAApaper_corrected';
title_cat='810nm_corr'; %'RedCh_corr'; 'NPcorr'; 'FRET'; '810nm_corr';

switch title_cat
    case 'RedCh_corr' %%%% mRuby correction dataset prolonged opto
        dirlist= {   %'Opto_constrict'         
            %'E:\Analysis\mouse400\aligned-mouse400-004\set4\analyzed\RatioCorr';
            % put folder names here
            };
        dirlist2= {   %'Ori8'         
            %'E:\Analysis\mouse400\aligned-mouse400-003\set4\analyzed\RatioCorr';
            % put folder names here
            };
        
    case 'FRET'
        dirlist= {        %'Opto_constrict'     
            %'E:\Analysis\mouse304\aligned-mouse304-009\set2\analyzed\RatioCorr';
            % put folder names here
            };
        dirlist2= {        %'Ori8'     
            %'E:\Analysis\mouse304\aligned-mouse304-006\set2\analyzed\RatioCorr';
            % put folder names here
            };
        
    case 'NPcorr'
        dirlist= {       %'Opto_constrict'      
            %'E:\Analysis\mouse292\aligned-mouse292-004\set4\analyzed\RatioCorr';
            % put folder names here
            };
        
        dirlist2= {     %'Ori8'
            %'E:\Analysis\mouse292\aligned-mouse292-001\set4\analyzed\RatioCorr';            
            % put folder names here
            };
        
    case '810nm_corr'
        dirlist= {       %'Opto_constrict'      
            %'E:\Analysis\mouse439\aligned-mouse439-008\set4\analyzed\RatioCorr';
            % put folder names here            
            };
        dirlist2= {       %'Ori8'      
            %'E:\Analysis\mouse400\aligned-mouse400-014\set4\analyzed\RatioCorr';
            % put folder names here
            };
        
    otherwise
        error('Unknown category: %s', category);
end

set(groot,...
    'defaultAxesFontName','Arial',...
    'defaultTextFontName','Arial',...
    'defaultLegendFontName','Arial',...
    'defaultColorbarFontName','Arial');

clear mouse orig_bef corr_bef orig_dur corr_dur orig_aft corr_aft indiv_bef_diff indiv_dur_diff indiv_aft_diff
run_type='Opto constrict';%
for ii=1:length(dirlist)
    mouse(ii).data=load([dirlist{ii} '\tcourse.mat']);
    idx = strfind(dirlist{ii}, 'mouse');
    mouse(ii).runID=dirlist{ii}(idx(2):idx(2) + length('mouse') + 6);
end

for ii=1:length(mouse)
    orig_bef(ii,1)=mouse(ii).data.pop_params.orig_bef_mean;
    corr_bef(ii,1)=mouse(ii).data.pop_params.corr_bef_mean;
    orig_dur(ii,1)=mouse(ii).data.pop_params.orig_dur_mean;
    corr_dur(ii,1)=mouse(ii).data.pop_params.corr_dur_mean;
    orig_aft(ii,1)=mouse(ii).data.pop_params.orig_aft_mean;
    corr_aft(ii,1)=mouse(ii).data.pop_params.corr_aft_mean;
    
    orig_bef_std(ii,1)=std(mouse(ii).data.pop_params.orig_bef_amp);
    corr_bef_std(ii,1)=std(mouse(ii).data.pop_params.corr_bef_amp);    
    orig_dur_std(ii,1)=std(mouse(ii).data.pop_params.orig_dur_amp);
    corr_dur_std(ii,1)=std(mouse(ii).data.pop_params.corr_dur_amp);
    orig_aft_std(ii,1)=std(mouse(ii).data.pop_params.orig_aft_amp);
    corr_aft_std(ii,1)=std(mouse(ii).data.pop_params.corr_aft_amp);

    indiv_bef_diff{ii}=mouse(ii).data.pop_params.orig_bef_amp - mouse(ii).data.pop_params.corr_bef_amp;
    indiv_dur_diff{ii}=mouse(ii).data.pop_params.orig_dur_amp - mouse(ii).data.pop_params.corr_dur_amp;
    indiv_aft_diff{ii}=mouse(ii).data.pop_params.orig_aft_amp - mouse(ii).data.pop_params.corr_aft_amp;    

    diff_summary(ii,:)=[sum(indiv_bef_diff{ii}>0) sum(indiv_bef_diff{ii}<0) sum(indiv_bef_diff{ii}==0) numel(indiv_bef_diff{ii}) ...
                        sum(indiv_dur_diff{ii}>0) sum(indiv_dur_diff{ii}<0) sum(indiv_dur_diff{ii}==0) numel(indiv_dur_diff{ii}) ...
                        sum(indiv_aft_diff{ii}>0) sum(indiv_aft_diff{ii}<0) sum(indiv_aft_diff{ii}==0) numel(indiv_aft_diff{ii})];
end

pop_orig_bef=mean(orig_bef);    pop_corr_bef=mean(corr_bef);
pop_orig_bef_std=std(orig_bef);    pop_corr_bef_std=std(corr_bef);
pop_means_bef = [pop_orig_bef, pop_corr_bef];    pop_stds_bef = [pop_orig_bef_std, pop_corr_bef_std];

pop_orig_dur=mean(orig_dur);    pop_corr_dur=mean(corr_dur);
pop_orig_dur_std=std(orig_dur);    pop_corr_dur_std=std(corr_dur);
pop_means_dur = [pop_orig_dur, pop_corr_dur];    pop_stds_dur = [pop_orig_dur_std, pop_corr_dur_std];

pop_orig_aft=mean(orig_aft);    pop_corr_aft=mean(corr_aft);
pop_orig_aft_std=std(orig_aft);    pop_corr_aft_std=std(corr_aft);
pop_means_aft = [pop_orig_aft, pop_corr_aft];    pop_stds_aft = [pop_orig_aft_std, pop_corr_aft_std];


%%% Make Table %%%
mytabledata=[   orig_bef orig_bef_std corr_bef corr_bef_std ...
            orig_dur orig_dur_std corr_dur corr_dur_std ...
            orig_aft orig_aft_std corr_aft corr_aft_std];
tabeltitle={'orig_bef' 'orig_bef_std' 'corr_bef' 'corr_bef_std' 'orig_dur' 'orig_dur_std' 'corr_dur' 'corr_dur_std' 'orig_aft' 'orig_aft_std' 'corr_aft' 'corr_aft_std'};
mytabledata(end+1,:)=[  pop_orig_bef pop_orig_bef_std pop_corr_bef pop_corr_bef_std ...
                    pop_orig_dur pop_orig_dur_std pop_corr_dur pop_corr_dur_std ...
                    pop_orig_aft pop_orig_aft_std pop_corr_aft pop_corr_aft_std];

rowtitles={mouse(:).runID};
rowtitles{end+1}='PopAvg';

myTable = array2table(mytabledata,'RowNames', rowtitles,'VariableNames', tabeltitle);
xlsfilename=[saveDir '\' title_cat '_' run_type '_SummaryTable.xlsx'];
writetable(myTable, xlsfilename, 'WriteRowNames', true)

difftabeltitle={'# cells decreased bef' '# cells increased bef' '# cells no change bef' '# cells total bef'...
                '# cells decreased dur' '# cells increased dur' '# cells no change dur' '# cells total dur'...
                '# cells decreased aft' '# cells increased aft' '# cells no change aft' '# cells total aft'}; 
my_diff_Table = array2table(diff_summary,'RowNames', rowtitles(1:end-1),'VariableNames', difftabeltitle);
xls_diff_filename=[saveDir '\' title_cat '_' run_type '_DiffSummaryTable.xlsx'];
writetable(my_diff_Table, xls_diff_filename, 'WriteRowNames', true)
%%%%%%%%%%%%%%%%%%

f1=figure; hold on;
plot(zeros(length(orig_bef),1)+.95,orig_bef,'go'); plot(zeros(length(orig_bef),1)+1.05,corr_bef,'bo');
plot(zeros(length(orig_dur),1)+1.95,orig_dur,'go'); plot(zeros(length(orig_dur),1)+2.05,corr_dur,'bo');
plot(zeros(length(orig_aft),1)+2.95,orig_aft,'go'); plot(zeros(length(orig_aft),1)+3.05,corr_aft,'bo');

xloc=[0.98,1.08; 1.98,2.08; 2.98,3.08];
colors=[0.3 0.7 0.3; 0 0 1];
for kk=1:size(xloc,2)
    errorbar(xloc(1,kk), pop_means_bef(kk), pop_stds_bef(kk), 's','MarkerSize', 10,'MarkerEdgeColor', colors(kk,:),'Color',colors(kk,:),'LineWidth', 1.5,'CapSize', 10);
    errorbar(xloc(2,kk), pop_means_dur(kk), pop_stds_dur(kk), 's','MarkerSize', 10,'MarkerEdgeColor', colors(kk,:),'Color',colors(kk,:),'LineWidth', 1.5,'CapSize', 10);
    errorbar(xloc(3,kk), pop_means_aft(kk), pop_stds_aft(kk), 's','MarkerSize', 10,'MarkerEdgeColor', colors(kk,:),'Color',colors(kk,:),'LineWidth', 1.5,'CapSize', 10);
end
xlim([0.5 3.5]);    xticks([1 2 3]);    xticklabels({'Bef', 'Dur', 'Aft'});
legend({'Orig'; 'Corrected'});  ylabel('Value');    title(['population Mean ± Std Dev ' title_cat ' ' run_type]);

p_non_param_bef = signrank(orig_bef, corr_bef);
p_non_param_dur = signrank(orig_dur, corr_dur);
p_non_param_aft = signrank(orig_aft, corr_aft);

[h_ttest_bef, p_ttest_bef, ci_bef, stats_bef] = ttest(orig_bef, corr_bef);
[h_ttest_dur, p_ttest_dur, ci_dur, stats_dur] = ttest(orig_dur, corr_dur);
[h_ttest_aft, p_ttest_aft, ci_aft, stats_aft] = ttest(orig_aft, corr_aft);

bef_diff=corr_bef - orig_bef;    dur_diff=corr_dur - orig_dur;    aft_diff=corr_aft - orig_aft;

cohenD_bef=mean(bef_diff)/std(bef_diff);
cohenD_dur=mean(dur_diff)/std(dur_diff);
cohenD_aft=mean(aft_diff)/std(aft_diff);

f2=figure; hold on; plot(zeros(length(bef_diff),1)+1,bef_diff,'go'); plot(zeros(length(bef_diff),1)+2,dur_diff,'go'); plot(zeros(length(bef_diff),1)+3,aft_diff,'go')
title(['population Diff Mean orig vs corr ' title_cat ' ' run_type]);
p_non_param_1 = signrank(bef_diff, dur_diff);

saveas(f1,[saveDir '\Pop distribution ' title_cat ' ' run_type '.eps'], 'psc2');
saveas(f1,[saveDir '\Pop distribution ' title_cat ' ' run_type '.fig']);
saveas(f1,[saveDir '\Pop distribution ' title_cat ' ' run_type '.tif']);
saveas(f2,[saveDir '\Pop distribution Diff ' title_cat ' ' run_type '.eps'], 'psc2');
saveas(f2,[saveDir '\Pop distribution Diff ' title_cat ' ' run_type '.fig']);
saveas(f2,[saveDir '\Pop distribution Diff ' title_cat ' ' run_type '.tif']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% for Ori8 data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear mouse2 orig_mean orig_std corr_mean corr_std orig_max corr_max indiv_mean_diff indiv_max_diff
run_type='Ori8'; 
for ii=1:length(dirlist2)
    mouse2(ii).data=load([dirlist2{ii} '\tcourse.mat']);
    idx = strfind(dirlist2{ii}, 'mouse');
    mouse2(ii).runID=dirlist2{ii}(idx(2):idx(2) + length('mouse') + 6);
end
for ii=1:length(mouse2)
    orig_mean(ii,1)=mean(mouse2(ii).data.pop_params.mean_resp_orig);
    orig_std(ii,1)=std(mouse2(ii).data.pop_params.mean_resp_orig);
    corr_mean(ii,1)=mean(mouse2(ii).data.pop_params.mean_resp_corr);
    corr_std(ii,1)=std(mouse2(ii).data.pop_params.mean_resp_corr);
    orig_max(ii,1)=mean(mouse2(ii).data.pop_params.max_resp_orig);
    orig_max_std(ii,1)=std(mouse2(ii).data.pop_params.max_resp_orig);
    corr_max(ii,1)=mean(mouse2(ii).data.pop_params.max_resp_corr);
    corr_max_std(ii,1)=std(mouse2(ii).data.pop_params.max_resp_corr);

    indiv_mean_diff{ii}=mouse2(ii).data.pop_params.mean_resp_orig - mouse2(ii).data.pop_params.mean_resp_corr;
    indiv_max_diff{ii}=mouse2(ii).data.pop_params.max_resp_orig - mouse2(ii).data.pop_params.max_resp_corr;
    
    diff_summary2(ii,:)=[sum(indiv_mean_diff{ii}>0) sum(indiv_mean_diff{ii}<0) sum(indiv_mean_diff{ii}==0) numel(indiv_mean_diff{ii}) ...
                        sum(indiv_max_diff{ii}>0) sum(indiv_max_diff{ii}<0) sum(indiv_max_diff{ii}==0) numel(indiv_max_diff{ii})];
end

pop_orig_mean=mean(orig_mean);    pop_corr_mean=mean(corr_mean);
pop_orig_std=std(orig_mean);    pop_corr_std=std(corr_mean);   
means = [pop_orig_mean, pop_corr_mean];    stds_mean = [pop_orig_std, pop_corr_std];

pop_orig_max=mean(orig_max);    pop_corr_max=mean(corr_max);
pop_orig_std_max=std(orig_max);    pop_corr_std_max=std(corr_max);   
maxs = [pop_orig_max, pop_corr_max];    stds_max = [pop_orig_std_max, pop_corr_std_max];


%%% Make Table %%%
mytabledata2=[  orig_mean orig_std corr_mean corr_std ...
                orig_max orig_max_std corr_max corr_max_std];
tabeltitle2={'orig_mean' 'orig_std' 'corr_mean' 'corr_std' 'orig_max' 'orig_max_std' 'corr_max' 'corr_max_std'};
mytabledata2(end+1,:)=[ pop_orig_mean pop_orig_std pop_corr_mean pop_corr_std ...
                        pop_orig_max pop_orig_std_max pop_corr_max pop_corr_std_max];
rowtitles2={mouse2(:).runID};
rowtitles2{end+1}='PopAvg';

myTable2 = array2table(mytabledata2,'RowNames', rowtitles2,'VariableNames', tabeltitle2);
xlsfilename2=[saveDir '\' title_cat '_' run_type '_SummaryTable.xlsx'];
writetable(myTable2, xlsfilename2, 'WriteRowNames', true)

difftabeltitle2={'# cells decreased mean' '# cells increased mean' '# cells no change mean' '# cells total mean'...
                '# cells decreased max' '# cells increased max' '# cells no change max' '# cells total max'}; 
my_diff_Table2 = array2table(diff_summary2,'RowNames', rowtitles2(1:end-1),'VariableNames', difftabeltitle2);
xls_diff_filename2=[saveDir '\' title_cat '_' run_type '_DiffSummaryTable.xlsx'];
writetable(my_diff_Table2, xls_diff_filename2, 'WriteRowNames', true)

%%%%%%%%%%%%%%%%%%%

f3=figure; hold on; subplot(1,2,1); hold on;
plot(zeros(length(orig_mean),1)+1,orig_mean,'go'); plot(zeros(length(corr_mean),1)+2,corr_mean,'bo');
xloc=[1.1,2.1];
colors=[0.3 0.7 0.3; 0 0 1];
for kk=1:size(xloc,2)
    errorbar(xloc(1,kk), means(kk), stds_mean(kk), 's','MarkerSize', 10,'MarkerEdgeColor', colors(kk,:),'Color',colors(kk,:),'LineWidth', 1.5,'CapSize', 10);
end
xlim([0.5 2.5]);    xticks([1 2]);    xticklabels({'Orig', 'Corr'});
ylabel('Value');    title(['Pop Mean DF/F ' title_cat ' ' run_type]);
subplot(1,2,2); hold on; 
plot(zeros(length(orig_max),1)+1,orig_max,'go'); plot(zeros(length(corr_max),1)+2,corr_max,'bo');
colors=[0.3 0.7 0.3; 0 0 1];
for kk=1:size(xloc,2)
    errorbar(xloc(1,kk), maxs(kk), stds_max(kk), 's','MarkerSize', 10,'MarkerEdgeColor', colors(kk,:),'Color',colors(kk,:),'LineWidth', 1.5,'CapSize', 10);
end
xlim([0.5 2.5]);    xticks([1 2]);    xticklabels({'Orig', 'Corr'});
ylabel('Value');    title(['Pop Max DF/F ' title_cat ' ' run_type]);


p_non_param_mean = signrank(pop_orig_mean, pop_corr_mean);
[h_ttest_mean, p_ttest_mean, ci_mean, stats_mean] = ttest(pop_orig_mean, pop_corr_mean);
mean_diff=corr_mean - orig_mean;    
cohenD_mean=mean(mean_diff)/std(mean_diff);

p_non_param_max = signrank(pop_orig_max, pop_corr_max);
[h_ttest_max, p_ttest_max, ci_max, stats_max] = ttest(pop_orig_max, pop_corr_max);
max_diff=corr_max - orig_max;    
cohenD_max=mean(max_diff)/std(max_diff);

f4=figure; hold on; subplot (1,2,1); hold on;
plot(zeros(length(mean_diff),1)+1,mean_diff,'go');
title(['Pop Diff Mean orig vs corr ' title_cat ' ' run_type]);
ax=gca; currentLimits = ax.YLim;
ymin1=min([0 currentLimits(1)]); ymax1=max([0 currentLimits(2)]);
ax.YLim = [ymin1, ymax1];
subplot (1,2,2); hold on;
plot(zeros(length(max_diff),1)+1,max_diff,'go');
title(['Pop Diff Max orig vs corr ' title_cat ' ' run_type]);
ax2=gca; currentLimits2 = ax2.YLim;
ymin2=min([0 currentLimits2(1)]); ymax2=max([0 currentLimits2(2)]);
ax2.YLim = [ymin2, ymax2];

saveas(f3,[saveDir '\Pop distribution ' title_cat ' ' run_type '.eps'], 'psc2');
saveas(f3,[saveDir '\Pop distribution ' title_cat ' ' run_type '.fig']);
saveas(f3,[saveDir '\Pop distribution ' title_cat ' ' run_type '.tif']);
saveas(f4,[saveDir '\Pop distribution Diff ' title_cat ' ' run_type '.eps'], 'psc2');
saveas(f4,[saveDir '\Pop distribution Diff ' title_cat ' ' run_type '.fig']);
saveas(f4,[saveDir '\Pop distribution Diff ' title_cat ' ' run_type '.tif']);













