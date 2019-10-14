clear;
addpath('.\Utils');
filepath='D:\Miniscope_Chenhaoshan\all_animal';
filenames=dir('D:\Miniscope_Chenhaoshan\all_animal\processed*.mat');
animal=cell(numel(filenames),1);
for n=1:numel(filenames)
    animal{n}=load(fullfile(filepath,filenames(n).name));
end

%%
accuracyA=zeros(numel(filenames),1);
accuracyB=zeros(numel(filenames),1);
for n=1:numel(filenames)
    sig=animal{n}.ms_dff.S_dff;
    
    [accuracyA(n),accuracyB(n)]=predict_ctx2(sig,animal{n}.session_start,animal{n}.protocol,animal{n}.ms.ms_ts);
end
%

%% Combine all shock response 
shock_responsemean=[];
shock_responsetrial=[];
SR_ids=cell(numel(filenames),1);
pre_dur=100; %100 frames= 10s 
post_dur=600;
shockdayresponse=[];
preresponse=[];
postresponse=[];
for n=1:numel(filenames)
    sigt=zscore(animal{n}.ms.sigraw',[],2);
    [shock_id,meanresponse_shock,shock_response]=find_SRcells(animal{n},pre_dur,post_dur);
    SR_ids{n}=shock_id;
    shock_responsetrial=[shock_responsetrial;shock_response];
    shock_responsemean=[shock_responsemean;meanresponse_shock];
    shockdayresponse=[shockdayresponse; sigt(shock_id,animal{n}.session_start(4)+animal{n}.frame_shock(1)-600:animal{n}.session_start(4)+animal{n}.frame_shock(1)+4000)];
    
end

binsize=10;
p_threshold=0.01;
zscorethreshold=2.56;
t=-pre_dur/10:0.1:post_dur/10-0.1;
%responsetemp=shock_responsemean(:,pre_dur:end);
responsetemp=shock_responsetrial(:,pre_dur:end,1);

responsetemp(responsetemp>zscorethreshold)=zscorethreshold;
[~,indmax]=max(responsetemp,[],2);
[~,indsort]=sort(indmax);


f=figure;
for s=1:3
    ax=subplot(1,4,s);
    imagesc('XData',t,'CData',shock_responsetrial(indsort,:,s),[0 zscorethreshold]);
    axis tight;
    ax.YLim=[0.5 length(indsort)+0.5];    
    hold on;
    area([0; 2],[ax.YLim; ax.YLim],'FaceAlpha',0.2,'FaceColor','r','LineStyle','none')
    title(['Shock #' num2str(s)]);
end
ax=subplot(1,4,4);
imagesc('XData',t,'CData',shock_responsemean(indsort,:),[0 zscorethreshold]);
axis tight;
ax.YLim=[0.5 length(indsort)+0.5]; 
hold on;
area([0; 2],[ax.YLim; ax.YLim],'FaceAlpha',0.2,'FaceColor','r','LineStyle','none')
title('Mean Response');

%%
FRpreA=[];
FRpostA=[];
FRpreB=[];
FRpostB=[];
pre_dur=100; %100 frames= 10s 
post_dur=100;
calcmethod=@calcAUC;
for n=1:numel(filenames)
    sig2=animal{n}.ms_dff.S_dff;
    %ids=1:size(sig2,1);
    ids=find_SRcells(animal{n},pre_dur,post_dur);
    session_start=animal{n}.session_start;
    session_end=animal{n}.session_end;
    protocol=animal{n}.protocol;
    FRpreA=[FRpreA; calcmethod(sig2(ids,session_start(strcmp(protocol,'preA')):session_end(strcmp(protocol,'preA'))))];
    FRpostA=[FRpostA; calcmethod(sig2(ids,session_start(strcmp(protocol,'postA')):session_end(strcmp(protocol,'postA'))))];
    FRpreB=[FRpreB; calcmethod(sig2(ids,session_start(strcmp(protocol,'preB')):session_end(strcmp(protocol,'preB'))))];
    FRpostB=[FRpostB; calcmethod(sig2(ids,session_start(strcmp(protocol,'postB')):session_end(strcmp(protocol,'postB'))))];
end
Y=[FRpreB FRpreA FRpostA FRpostB];
g={'preB','preA','postA','postB'};
%%% 1-way ANOVA

% [~,~,stats]=anova1(Y,g);
% [c,~,~,gnames] = multcompare(stats);

% paired T-test
% [~,p1]=ttest(FRpreA,FRpostA);
% [~,p2]=ttest(FRpreB,FRpostB);
% [~,p3]=ttest(FRpreA,FRpreB);
% [~,p4]=ttest(FRpostA,FRpostB);

[p1,h1]=ranksum(FRpreA,FRpostA);
[p2,h2]=ranksum(FRpreB,FRpostB);
[p3,h3]=ranksum(FRpreA,FRpreB);
[p4,h4]=ranksum(FRpostA,FRpostB);
group_pair={[2,3],[1,4],[2,1],[3,4]};
p_pair=[p1 p2 p3 p4];

figure
boxplot(Y,'Notch','on','Labels',g);
hold on;
sigstar(group_pair((p_pair<=0.05)),p_pair(p_pair<=0.05));
ax=gca;
% ax.YLim=[0 0.3];

%%
function [shock_id,meanresponse_shock,shock_response]=find_SRcells(animaldata,pre_dur,post_dur)
    sig=animaldata.ms.sigraw';
    shockts=animaldata.shockts;
    ms_start=animaldata.ms_start;
    shock_start=animaldata.shock_start;
    ms_ts=animaldata.ms.ms_ts;
    session_dur=cellfun(@length,ms_ts);
    session_start=1;
    for i=1:length(session_dur)-1
        starttemp=session_start(i)+session_dur(i);
        session_start=[session_start; starttemp];
    end
    session_end=session_start+session_dur'-1;
    
    % Get shock time
    shock_ts2=shockts.time-shockts.time(ms_start(4));
    shock_bvt=shock_ts2(shock_start);
    shock_mst=ms_ts{4}/1000;
    temp1=abs(double(repmat(shock_mst,length(shock_bvt))')-shock_bvt');
    [minval,frame_shock]=min(temp1,[],1);
    

    binsize=10;
    p_threshold=0.01;
    zscorethreshold=2.56;
    t=-pre_dur/10:0.1:post_dur/10-0.1;

    cdn_sigtemp=sig(:,session_start(4):session_end(4));
    for n=1:size(cdn_sigtemp,1)
        cdn_sig(n,:)=zscore(cdn_sigtemp(n,:));
    end

    %zscore across all days
    % for n=1:size(cdn_sigtemp,1)
    %     sigz(n,:)=zscore(sig(n,:));
    % end
    % cdn_sig=sigz(:,session_start(4):session_end(4));

    shock_response=zeros(size(cdn_sig,1),pre_dur+post_dur,3);
    shock_responsez=zeros(size(cdn_sig,1),pre_dur+post_dur,3);
    postshockresponsez=zeros(size(cdn_sig,1),post_dur,3);
    binresponsepre=[];
    binresponsepost=[];
    for s=1:numel(frame_shock)
        shock_response(:,:,s)=cdn_sig(:,frame_shock(s)-pre_dur:frame_shock(s)+post_dur-1);
        preshocktemp=cdn_sig(:,frame_shock(s)-pre_dur:frame_shock(s)-1);
        postshocktemp=cdn_sig(:,frame_shock(s):frame_shock(s)+post_dur-1);
        postshockresponsez(:,:,s)=postshocktemp;
        binresponsepre=[binresponsepre avgbinresponse(pre_dur,10,preshocktemp)];
        binresponsepost=[binresponsepost avgbinresponse(post_dur,10,postshocktemp)];
    end
    shock_responsemean=mean(shock_response,3);


    % find significantly excited units 
    p=zeros(size(sig,1),1);
    for n=1:size(sig,1)
        p(n)=ranksum(binresponsepre(n,:),binresponsepost(n,:),'tail','left');
    end
    idtemp=find(p<=0.05);
    responsetemp=shock_responsemean(idtemp,pre_dur:end);
    responsetemp(responsetemp>zscorethreshold)=zscorethreshold;
    [~,indmax]=max(responsetemp,[],2);
    [~,indsort]=sort(indmax);
    shock_id=idtemp(indsort);
    meanresponse_shock=shock_responsemean(shock_id,:);
    shock_response=shock_response(shock_id,:,:);
    
end

