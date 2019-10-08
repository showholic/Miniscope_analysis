clear;
addpath('.\Utils');
filepath='E:\Miniscope_Chenhaoshan\all_animal';
filenames=dir('E:\Miniscope_Chenhaoshan\all_animal\*.mat');
animal=cell(numel(filenames),1);
for n=1:numel(filenames)
    animal{n}=load(fullfile(filepath,filenames(n).name));
end
%% FLM without motion restriction 
accuracyA=zeros(numel(filenames),1);
accuracyB=zeros(numel(filenames),1);
for n=1:numel(filenames)
    sig=animal{n}.ms.sigraw';
    
    [accuracyA(n),accuracyB(n)]=predict_ctx(sig,animal{n}.session_start,animal{n}.protocol,animal{n}.ms.ms_ts);
end
%% Combine all shock response 
shock_responsemean=[];
shock_responsetrial=[];
SR_ids=cell(numel(filenames),1);
pre_dur=50; %100 frames= 10s 
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
figure;
imagesc(shockdayresponse(indsort,:),[0,3]);


%% compare pre and post context preference 
ctxp1=[];
ctxp2=[];
for n=1:numel(filenames)
    animaldata=animal{n};
    sig2=animaldata.ms.sigdeconvolved';
    %ids=find_SRcells(animaldata,pre_dur,post_dur); %all SR cells 
    ids=1:size(sig2,1); %all cells
    session_start=animaldata.session_start;
    session_end=animaldata.session_end;
    protocol=animaldata.protocol;
    FRpreA=calcFR(sig2(ids,session_start(strcmp(protocol,'preA')):session_end(strcmp(protocol,'preA'))));
    FRpostA=calcFR(sig2(ids,session_start(strcmp(protocol,'postA')):session_end(strcmp(protocol,'postA'))));
    FRpreB=calcFR(sig2(ids,session_start(strcmp(protocol,'preB')):session_end(strcmp(protocol,'preB'))));
    FRpostB=calcFR(sig2(ids,session_start(strcmp(protocol,'postB')):session_end(strcmp(protocol,'postB'))));
    Y=[FRpreB FRpreA FRpostA FRpostB];
    ctxpref1=(FRpreA-FRpreB)./(FRpreA+FRpreB);
    ctxpref2=(FRpostA-FRpostB)./(FRpostA+FRpostB);
    ctxp1=[ctxp1; abs(ctxpref1)];
    ctxp2=[ctxp2; abs(ctxpref2)];
end

figure;
cdfplot(ctxp1);
hold on;
cdfplot(ctxp2);
legend('Pre-conditioning','Post-conditioning');
[h,p] = kstest2(ctxp1,ctxp2,'Tail','larger');
%%
ctxp1=[];
ctxp2=[];
sigall=[];
for n=1:numel(filenames)
    animaldata=animal{n};
    sig2=animaldata.ms.sigdeconvolved';
    ids=find_SRcells(animaldata,pre_dur,post_dur); %all SR cells 
    %ids=1:size(sig2,1); %all cells
    session_start=animaldata.session_start;
    session_end=animaldata.session_end;
    protocol=animaldata.protocol;
    FRpreA=calcFR(sig2(ids,session_start(strcmp(protocol,'preA')):session_end(strcmp(protocol,'preA'))));
    FRpostA=calcFR(sig2(ids,session_start(strcmp(protocol,'postA')):session_end(strcmp(protocol,'postA'))));
    FRpreB=calcFR(sig2(ids,session_start(strcmp(protocol,'preB')):session_end(strcmp(protocol,'preB'))));
    FRpostB=calcFR(sig2(ids,session_start(strcmp(protocol,'postB')):session_end(strcmp(protocol,'postB'))));
    Y=[FRpreB FRpreA FRpostA FRpostB];
    ctxpref1=(FRpreA-FRpreB)./(FRpreA+FRpreB);
    ctxpref2=(FRpostA-FRpostB)./(FRpostA+FRpostB);
    ctxp1=[ctxp1; (ctxpref1)];
    ctxp2=[ctxp2; (ctxpref2)];
    sig1temp=zscore(animaldata.ms.sigraw(:,ids)',[],2);
    sig1=[];
    for s=1:numel(animaldata.session_start)
        sig1=[sig1 sig1temp(:,animaldata.session_start(s):animaldata.session_start(s)+2900)];
    end
    sigall=[sigall;sig1];
end
%%
t=1:2900:size(sigall,2);
figure;

imagesc(sigall(inddiff,:),[0,3]);
hold on;
ax=gca;
for i=1:length(t)-1
    plot([t(i+1),t(i+1)],ax.YLim,'--k');
end

figure;
indtemp=inddiff(350:end);
sigt=sigall(indtemp,:);
sigt=sigt./max(sigt,[],2);
plot((sigt+(1:size(sigt,1))')');
hold on;
ax=gca;

for i=1:length(t)-1
    plot([t(i+1),t(i+1)],ax.YLim,'--k');
end
%% Non-SR cells
ctxp1=[];
ctxp2=[];
for n=1:numel(filenames)
    animaldata=animal{n};
    sig2=animaldata.ms.sigdeconvolved';
    SRcell=find_SRcells(animaldata); %all SR cells 
    allcell=(1:size(sig2,1))'; %all cells
    ids=setdiff(allcell,SRcell);
    session_start=animaldata.session_start;
    session_end=animaldata.session_end;
    protocol=animaldata.protocol;
    FRpreA=calcFR(sig2(ids,session_start(strcmp(protocol,'preA')):session_end(strcmp(protocol,'preA'))));
    FRpostA=calcFR(sig2(ids,session_start(strcmp(protocol,'postA')):session_end(strcmp(protocol,'postA'))));
    FRpreB=calcFR(sig2(ids,session_start(strcmp(protocol,'preB')):session_end(strcmp(protocol,'preB'))));
    FRpostB=calcFR(sig2(ids,session_start(strcmp(protocol,'postB')):session_end(strcmp(protocol,'postB'))));
    Y=[FRpreB FRpreA FRpostA FRpostB];
    ctxpref1=(FRpreA-FRpreB)./(FRpreA+FRpreB);
    ctxpref2=(FRpostA-FRpostB)./(FRpostA+FRpostB);
    ctxp1=[ctxp1; abs(ctxpref1)];
    ctxp2=[ctxp2; abs(ctxpref2)];
end

figure;
cdfplot(ctxp1);
hold on;
cdfplot(ctxp2);
legend('Pre-conditioning','Post-conditioning');
[h,p] = kstest2(ctxp1,ctxp2,'Tail','larger');
%% Identify pre-conditioning context cells through bootstrap (SLOW !!!!!)
ctx_dur=3000-120;
shocksession_dur=3500;
n_shuffle=100;
calcmethod=@calcAUC;
splitsize=6;
ctxpref_pre=cell(numel(filenames),1);
ctxpref_post=cell(numel(filenames),1);
sigpre=cell(numel(filenames),1);
sigpost=cell(numel(filenames),1);
ctxA_ind=cell(numel(filenames),1);
ctxB_ind=cell(numel(filenames),1);

for y=1:numel(filenames)
    sig2=animal{y}.ms.sigraw';   
    ts = tinv([0.01  0.99],size(sig2,1)-1);
    session_start=animal{y}.session_start;
    protocol=animal{y}.protocol;
    ids=1:size(sig2,1);
    sigpreA=sig2(ids,session_start(strcmp(protocol,'preA')):session_start(strcmp(protocol,'preA'))+ctx_dur-1);
    sigpostA=sig2(ids,session_start(strcmp(protocol,'postA')):session_start(strcmp(protocol,'postA'))+ctx_dur-1);
    %sigshock=sig2(ids,session_start(strcmp(protocol,'conditioning')):session_start(strcmp(protocol,'conditioning'))+shocksession_dur-1);
    sigpreB=sig2(ids,session_start(strcmp(protocol,'preB')):session_start(strcmp(protocol,'preB'))+ctx_dur-1);
    sigpostB=sig2(ids,session_start(strcmp(protocol,'postB')):session_start(strcmp(protocol,'postB'))+ctx_dur-1);
    FRpreA=calcmethod(sigpreA);
    FRpostA=calcmethod(sigpostA);
    FRpreB=calcmethod(sigpreB);
    FRpostB=calcmethod(sigpostB);
    ctxpref_pre{y}=(FRpreA-FRpreB)./(FRpreA+FRpreB);
    ctxpref_post{y}=(FRpostA-FRpostB)./(FRpostA+FRpostB);
    sigpre{y}=[sigpreB sigpreA];
    sigpost{y}=[sigpostA sigpostB];    
    ctxprefpre_shuffle=zeros(size(sig2,1),n_shuffle);
    for x=1:n_shuffle
        randshift=randi(500,1);
        sigshift=circshift(sigpre{y},randshift,2); %shuffle pre-conditioning transients 
        randsplit=randperm(splitsize);
        sigsplit=reshape(sigshift,size(sigshift,1),[],splitsize);
        sigshuffle=reshape(sigsplit(:,:,randsplit),size(sigshift,1),[]);
        frb=calcmethod(sigshuffle(:,1:ctx_dur));
        fra=calcmethod(sigshuffle(:,ctx_dur+1:end));
        ctxprefpre_shuffle(:,x)=(fra-frb)./(fra+frb);
    end
    ctxsigpref=zeros(size(ctxprefpre_shuffle,1),1);
    
    for n=1:size(ctxprefpre_shuffle,1)
        ctxpreftemp=ctxprefpre_shuffle(n,:);
        SEM = std(ctxpreftemp)/sqrt(length(ctxpreftemp));               % Standard Error
        
        CI = mean(ctxpreftemp) + ts*SEM;
        if ctxpref_pre{y}(n)>=CI(2)
            ctxsigpref(n)=1;
        elseif ctxpref_pre{y}(n)<=CI(1)
            ctxsigpref(n)=-1;
        else
            ctxsigpref(n)=0;
        end
    end
    ctxA_ind{y}=find(ctxsigpref==1);
    ctxB_ind{y}=find(ctxsigpref==-1);
end
%%
figure;
hold on;
ctxprefpreA_all=[];
ctxprefpostA_all=[];
ctxprefpreB_all=[];
ctxprefpostB_all=[];
for n=1:6
    plot_pair(ctxpref_pre{n},ctxpref_post{n},ctxA_ind{n},'r')
    plot_pair(ctxpref_pre{n},ctxpref_post{n},ctxB_ind{n},'b')
    ctxprefpreA_all=[ctxprefpreA_all; ctxpref_pre{n}(ctxA_ind{n})];
    ctxprefpostA_all=[ctxprefpostA_all; ctxpref_post{n}(ctxA_ind{n})];
    ctxprefpreB_all=[ctxprefpreB_all; ctxpref_pre{n}(ctxB_ind{n})];
    ctxprefpostB_all=[ctxprefpostB_all; ctxpref_post{n}(ctxB_ind{n})];
end
[~,p1]=ttest(ctxprefpreA_all,ctxprefpostA_all,'tail','right');
[~,p2]=ttest(ctxprefpreB_all,ctxprefpostB_all,'tail','left');
%% 
figure;
hold on;
ctxprefpreA_all=[];
ctxprefpostA_all=[];
ctxprefpreB_all=[];
ctxprefpostB_all=[];
for n=[3,5,6]
    plot_pair(ctxpref_pre{n},ctxpref_post{n},ctxA_ind{n},'r')
    plot_pair(ctxpref_pre{n},ctxpref_post{n},ctxB_ind{n},'b')
    ctxprefpreA_all=[ctxprefpreA_all; ctxpref_pre{n}(ctxA_ind{n})];
    ctxprefpostA_all=[ctxprefpostA_all; ctxpref_post{n}(ctxA_ind{n})];
    ctxprefpreB_all=[ctxprefpreB_all; ctxpref_pre{n}(ctxB_ind{n})];
    ctxprefpostB_all=[ctxprefpostB_all; ctxpref_post{n}(ctxB_ind{n})];
end
[~,p1]=ttest(ctxprefpreA_all,ctxprefpostA_all,'tail','right');
[~,p2]=ttest(ctxprefpreB_all,ctxprefpostB_all,'tail','left');
%% Context A cells, 2 different post-conditioning sequences 
figure;
i=1;
for n=[1,2,4]
    subplot(3,2,(i-1)*2+1)
    imagesc(zscore(sigpre{n}(ctxA_ind{n},:),[],2),[1.65 6]);
    subplot(3,2,(i-1)*2+2)
    imagesc(zscore(sigpost{n}(ctxA_ind{n},:),[],2),[1.65 6]);
    i=i+1;
end

figure;
i=1;
for n=[3,5,6]
    subplot(3,2,(i-1)*2+1)
    imagesc(zscore(sigpre{n}(ctxA_ind{n},:),[],2),[1.65 6]);
    subplot(3,2,(i-1)*2+2)
    imagesc(zscore(sigpost{n}(ctxA_ind{n},:),[],2),[1.65 6]);
    i=i+1;
end
%% Context B cells
figure;
for n=1:numel(filenames)
    subplot(6,2,(n-1)*2+1)
    imagesc(zscore(sigpre{n}(ctxB_ind{n},:),[],2),[1.65 6]);
    subplot(6,2,(n-1)*2+2)
    imagesc(zscore(sigpost{n}(ctxB_ind{n},:),[],2),[1.65 6]);
end

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

