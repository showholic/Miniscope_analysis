clear;
filepath='E:\Miniscope_Chenhaoshan\all_animal';
filenames=dir('E:\Miniscope_Chenhaoshan\all_animal\*.mat');
animal=cell(numel(filenames),1);
for n=1:numel(filenames)
    animal{n}=load(fullfile(filepath,filenames(n).name));
end

%%
shock_responsemean=[];
shock_responsetrial=[];
for n=1:numel(filenames)
    [meanresponse_shock,shock_response]=find_SRcells(animal{n});
    shock_responsetrial=[shock_responsetrial;shock_response];
    shock_responsemean=[shock_responsemean;meanresponse_shock];
end
pre_dur=100; %100 frames= 10s 
post_dur=100;
binsize=10;
p_threshold=0.01;
zscorethreshold=2.56;
t=-pre_dur/10:0.1:post_dur/10-0.1;
responsetemp=shock_responsemean(:,pre_dur:end);
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
idtemp=find_SRcells(animal{1});
function [meanresponse_shock,shock_response]=find_SRcells(animaldata)
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
    
    pre_dur=100; %100 frames= 10s 
    post_dur=100;
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

function binresponse=avgbinresponse(dur,binsize,response)
    bins=discretize(1:dur,0:binsize:dur);
    bincount=dur/binsize;
    binresponse=zeros(size(response,1),bincount);
    for i = 1:bincount
        binresponse(:,i)=mean(response(:,bins==i),2);
    end
end