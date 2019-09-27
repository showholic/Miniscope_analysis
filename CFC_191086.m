clear;
close all;
filepath='E:\Miniscope_Chenhaoshan\Results_191086\20190912_155106';
load(fullfile(filepath,'ms.mat'));
shockts=readtable(fullfile(filepath,'shock_behavts.csv'));
ms_start=[20;19;20;29;28;20];
shock_start=[5434;7236;9038];
protocol={'Home','preB','preA','conditioning','postA','postB'};
shock_dur=2/(1/10); %in frames
%% Overview
figure;
roi=ms.SFP;
sig=ms.sigraw';
ms_ts=ms.ms_ts;
imax=max(ms.SFP,[],3);
subplot(121)
imshow(imax,[0 0.4])
subplot(122)
numplot=10;
sigt=sig(1:numplot,:);
sigt=sigt./max(sigt,[],2);
plot((sigt+(1:size(sigt,1))')');
hold on;
session_dur=cellfun(@length,ms_ts);
session_start=1;
for i=1:length(session_dur)-1
    starttemp=session_start(i)+session_dur(i);
    session_start=[session_start; starttemp];
    plot([starttemp starttemp],ylim,'--k');
end
session_end=session_start+session_dur'-1;

% Get shock time
shock_ts2=shockts.time-shockts.time(ms_start(4));
shock_bvt=shock_ts2(shock_start);
shock_mst=ms_ts{4}/1000;
temp1=abs(double(repmat(shock_mst,length(shock_bvt))')-shock_bvt');
[minval,frame_shock]=min(temp1,[],1);

%verify shock end time 
shock_bvt2=shock_ts2(shock_start)+2;
temp1=abs(double(repmat(shock_mst,length(shock_bvt2))')-shock_bvt2');
[minval,frame_shock_end]=min(temp1,[],1);

for s=1:numel(frame_shock)
    area([session_start(4)+frame_shock(s); session_start(4)+frame_shock(s)+shock_dur],[ylim; ylim],'FaceAlpha',1,'FaceColor','r','LineStyle','none')
end
    
%% Conditioning day
pre_dur=100; %100 frames= 10s 
post_dur=100;

cdn_sig=sig(:,session_start(4):session_end(4));

% sigt=cdn_sig(1:numplot,:);
% sigt=sigt./max(sigt,[],2);
% figure;
% plot((sigt+(1:size(sigt,1))')');
% hold on;
% for s=1:numel(frame_shock)
%     area([frame_shock(s); frame_shock(s)+shock_dur],[ylim; ylim],'FaceAlpha',0.8,'FaceColor','r','LineStyle','none')
% end

baseline=cdn_sig(:,frame_shock(1)-pre_dur:frame_shock(1)-1);
baseline_mean=mean(baseline,2);
baseline_std=std(baseline,0,2);
shock_response=zeros(size(baseline,1),post_dur,3);
shock_responsez=zeros(size(baseline,1),post_dur,3);
for s=1:numel(frame_shock)
    shock_response(:,:,s)=cdn_sig(:,frame_shock(s):frame_shock(s)+post_dur-1);
    shock_responsez(:,:,s)=(shock_response(:,:,s)-baseline_mean)./baseline_std;
end
shock_responsemean=mean(shock_responsez,3);
% mean shock response 
[valmax,indmax]=max(shock_responsemean,[],2);
[~,indsort]=sort(indmax);

figure;
for s=1:numel(frame_shock)
    subplot(1,4,s);
    imagesc(shock_responsez(indsort,:,s),[1.65,6]);
    hold on;
    area([0; 0+shock_dur],[ylim; ylim],'FaceAlpha',0.2,'FaceColor','r','LineStyle','none')
    title(['Shock #' num2str(s)]);
end
subplot(1,4,4);
imagesc(shock_responsemean(indsort,:),[1.65,6]);
hold on;
area([0; 0+shock_dur],[ylim; ylim],'FaceAlpha',0.2,'FaceColor','r','LineStyle','none')
title('Mean Response');

figure;
plot(sum(shock_responsemean(indsort,:),1)/size(shock_responsemean,1));

%% Shock responsive cell (rank-sum test) 
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

shock_response=zeros(size(baseline,1),pre_dur+post_dur,3);
shock_responsez=zeros(size(baseline,1),pre_dur+post_dur,3);
postshockresponsez=zeros(size(baseline,1),post_dur,3);
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
SR_id_increase=idtemp(indsort);
f=figure;
for s=1:numel(frame_shock)
    ax=subplot(1,4,s);
    imagesc('XData',t,'CData',shock_response(SR_id_increase,:,s),[0 zscorethreshold]);
    axis tight;
    ax.YLim=[0.5 length(SR_id_increase)+0.5];    
    hold on;
    area([0; 2],[ax.YLim; ax.YLim],'FaceAlpha',0.2,'FaceColor','r','LineStyle','none')
    title(['Shock #' num2str(s)]);
end
ax=subplot(1,4,4);
imagesc('XData',t,'CData',shock_responsemean(SR_id_increase,:),[0 zscorethreshold]);
axis tight;
ax.YLim=[0.5 length(SR_id_increase)+0.5]; 
hold on;
area([0; 2],[ax.YLim; ax.YLim],'FaceAlpha',0.2,'FaceColor','r','LineStyle','none')
title('Mean Response');

% find significantly inhibited units
for n=1:size(sig,1)
    p(n)=ranksum(binresponsepre(n,:),binresponsepost(n,:),'tail','right');
end
idtemp=find(p<=p_threshold);
responsetemp=shock_responsemean(idtemp,1:pre_dur-1);
responsetemp(responsetemp>zscorethreshold)=zscorethreshold;
[valmax,indmax]=max(responsetemp,[],2);
[~,indsort]=sort(indmax);
SR_id_decrease=idtemp(indsort);
figure;
for s=1:numel(frame_shock)
    ax=subplot(1,4,s);
    imagesc('XData',t,'CData',shock_response(SR_id_decrease,:,s),[0 zscorethreshold]);
    axis tight;
    ax.YLim=[0.5 length(SR_id_decrease)+0.5];    
    hold on;
    area([0; 2],[ax.YLim; ax.YLim],'FaceAlpha',0.2,'FaceColor','r','LineStyle','none')
    title(['Shock #' num2str(s)]);
end
ax=subplot(1,4,4);
imagesc('XData',t,'CData',shock_responsemean(SR_id_decrease,:),[0 zscorethreshold]);
axis tight;
ax.YLim=[0.5 length(SR_id_decrease)+0.5];  
hold on;
area([0; 2],[ax.YLim; ax.YLim],'FaceAlpha',0.2,'FaceColor','r','LineStyle','none')
title('Mean Response');
saveas(f,fullfile(filepath,'shock_response.png'))
save('E:\Miniscope_Chenhaoshan\all_animal\processed_191086.mat','shockts','ms_start', ...
    'shock_start','protocol','ms','shock_response','shock_responsemean');

%% Define shock responsive cells using shuffling 
cdn_sig = sig(:,session_start(4):session_start(4)+3599);
cdn_sig(cdn_sig<0)=0;
pre_dur=100;
post_dur=100;
n_shuffle=1000;
AUC_shock=zeros(size(cdn_sig,1),numel(frame_shock));
AUC_shock_shuffle=zeros(n_shuffle,size(cdn_sig,1),numel(frame_shock));
for s=1:numel(frame_shock)
    postshocktemp=cdn_sig(:,frame_shock(s):frame_shock(s)+post_dur-1);
    for n=1:size(cdn_sig)
        AUC_shock(n,s)=trapz(postshocktemp(n,:));
    end
end
% Shuffle signals (not recommended)
for x= 1: n_shuffle
    for n=1:size(cdn_sig)
        sigshuffle=shuffle_sig(cdn_sig(n,:),500,6);
        for s=1:numel(frame_shock)
            postshocktemp=sigshuffle(:,frame_shock(s):frame_shock(s)+post_dur-1);
            AUC_shock_shuffle(x,n,s)=trapz(postshocktemp);
        end
    end
end
%% Visualize bootstrap shock
shocksigre=zeros(size(cdn_sig,1),numel(frame_shock));
%ts = tinv([0.01  0.99],length(AUCshuffle)-1);      % T-Score
ts=[-1.65 1.65];
for n=1:size(cdn_sig,1)
    for s=1:numel(frame_shock)
        AUCshuffle=AUC_shock_shuffle(:,n,s);
        AUC=AUC_shock(n,s);
        SEM = std(AUCshuffle);%/sqrt(length(AUCshuffle));               % Standard Error
       
        CI = mean(AUCshuffle) + ts*SEM;    
        if AUC>=CI(2)
            shocksigre(n,s)=1;
        elseif AUC<=CI(1)
            shocksigre(n,s)=-1;
        else
            shocksigre(n,s)=0;
        end
    end
end
idtemp=find(sum(shocksigre,2)>=1); %Chnage this number to select robustness 
f=figure;
responsetemp=shock_responsemean(idtemp,pre_dur+1:end);
[valmax,indmax]=max(responsetemp,[],2);
[~,indsort]=sort(indmax);
SR_shuffle_id=idtemp(indsort);
for s=1:numel(frame_shock)
    ax=subplot(1,4,s);
    imagesc('XData',t,'CData',shock_response(SR_shuffle_id,:,s),[0 zscorethreshold]);
    axis tight;
    ax.YLim=[0.5 length(SR_shuffle_id)+0.5];    
    hold on;
    area([0; 2],[ax.YLim; ax.YLim],'FaceAlpha',0.2,'FaceColor','r','LineStyle','none')
    title(['Shock #' num2str(s)]);
end
ax=subplot(1,4,4);
imagesc('XData',t,'CData',shock_responsemean(SR_shuffle_id,:),[0 zscorethreshold]);
axis tight;
ax.YLim=[0.5 length(SR_shuffle_id)+0.5]; 
hold on;
area([0; 2],[ax.YLim; ax.YLim],'FaceAlpha',0.2,'FaceColor','r','LineStyle','none')
title('Mean Response');
%% Use the same baseline 
baseline_dur=200;
baseline=avgbinresponse(baseline_dur,binsize,cdn_sig(:,frame_shock(1)-baseline_dur:frame_shock(1)-1));
for s=1:numel(frame_shock)
    for n=1:size(sig,1)
        %p_ishock(n,s)=ranksum(baseline(n,:),binresponsepost(n,(s-1)*(post_dur/binsize)+1:s*(post_dur/binsize)),'tail','left');
        maxz_ishock(n,s)=max(postshockresponsez(n,:,s));
    end
end
id_ishock=cell(numel(frame_shock),1);
figure;
for s=1:numel(frame_shock)
    %id_ishock{s}=find(p_ishock(:,s)<=0.01);
    id_ishock{s}=find(maxz_ishock(:,s)>=3);
    ax=subplot(1,3,s);
    imagesc('XData',t,'CData',shock_response(id_ishock{s},:,s),[0 6]);
    axis tight;
    ax.YLim=[0.5 length(id_ishock{s})+0.5];    
    hold on;
    area([0; 2],[ax.YLim; ax.YLim],'FaceAlpha',0.2,'FaceColor','r','LineStyle','none');
    title(['Shock #' num2str(s)]);
end

figure;

Boundary={};
for s=1:numel(frame_shock)
    idtemp=id_ishock{s};
    roitemp=ms.SFP(:,:,idtemp);
    subplot(1,3,s)
    hold on;
    for n=1:numel(idtemp)
        roitempp=roitemp(:,:,n);
        roitempp((roitempp<=(1- 0.6)*max(roitempp,[],'all')))=0;
        B=bwboundaries(roitempp);
        Boundary{n,s}=B{1};
        plot(Boundary{n,s}(:,2),Boundary{n,s}(:,1));
    end
    axis equal;
    ax=gca;
    ax.XLim=[0 ms.width];
    ax.YLim=[0 ms.height];
end

%% Video 
v=load(fullfile(filepath,'motioncorrected.mat'));
%%
% frameno=session_start(4)+frame_shock(1);
% v.CurrentTime=frameno*(1/v.FrameRate);
pre_dur=20;
post_dur=100;
vid_dur=pre_dur+post_dur;
im=zeros(ms.height,ms.width,vid_dur,3);

for s=1:3    
    for i=1:vid_dur
        frameno=session_start(4)+frame_shock(s)-pre_dur+i-1;
        im(:,:,i,s)=squeeze(v.vid(frameno,:,:));
    end
end

vid_shock=[im(:,:,:,1),im(:,:,:,2),im(:,:,:,3)];
v_out=VideoWriter(fullfile(filepath,'Shock_video.avi'),'Uncompressed AVI');
v_out.FrameRate=10;
open(v_out);
figure;
img=imshow(vid_shock(:,:,1),[0 250]);
axis equal;
axis tight;
for i=1:vid_dur
    img.CData=vid_shock(:,:,i);
    if i==pre_dur+1
        hold on;
        sh=scatter(ms.width*1.5,230,10,'filled','MarkerFaceColor','r');
%         for s=1:3
%             for n=1:numel(id_ishock{s})
%                 plot(Boundary{n,s}(:,2)+(s-1)*double(ms.width),Boundary{n,s}(:,1));
%             end
%         end
    elseif i==pre_dur+shock_dur+1
        delete(sh);
    end
    writeVideo(v_out,getframe(gca));
end
close(v_out);
%%
max_proj=max(vid_shock(:,:,pre_dur+1:end),[],3);
figure;
imshow(max_proj,[0 200]);
hold on;
for s=1:3
    for n=1:numel(id_ishock{s})
        plot(Boundary{n,s}(:,2)+(s-1)*double(ms.width),Boundary{n,s}(:,1));
    end
end 
%% Compare FR 
sig2=ms.sigdeconvolved';
% All cells
compare_session_FR(1:size(sig2,1),sig2,session_start,session_end,protocol)
% Shock activated cells 
compare_session_FR(SR_id_increase,sig2,session_start,session_end,protocol)
% Shock decreased cells 
compare_session_FR(SR_id_decrease,sig2,session_start,session_end,protocol)





%% Define pre-conditioning context cells using shuffling
sig2=ms.sigraw';
calcmethod=@calcAUC;
ids=1:size(sig2,1);

ctx_dur=3000-60;
shocksession_dur=3500;
sigpreA=sig2(ids,session_start(strcmp(protocol,'preA')):session_start(strcmp(protocol,'preA'))+ctx_dur-1);
sigpostA=sig2(ids,session_start(strcmp(protocol,'postA')):session_start(strcmp(protocol,'postA'))+ctx_dur-1);
sigshock=sig2(ids,session_start(strcmp(protocol,'conditioning')):session_start(strcmp(protocol,'conditioning'))+shocksession_dur-1);
sigpreB=sig2(ids,session_start(strcmp(protocol,'preB')):session_start(strcmp(protocol,'preB'))+ctx_dur-1);
sigpostB=sig2(ids,session_start(strcmp(protocol,'postB')):session_start(strcmp(protocol,'postB'))+ctx_dur-1);
FRpreA=calcmethod(sigpreA);
FRpostA=calcmethod(sigpostA);
FRpreB=calcmethod(sigpreB);
FRpostB=calcmethod(sigpostB);
ctxpref_pre=(FRpreA-FRpreB)./(FRpreA+FRpreB);
ctxpref_post=(FRpostA-FRpostB)./(FRpostA+FRpostB);

calcmethod=@calcTF;
sigpre=[sigpreB sigpreA];
sigpost=[sigpostA sigpostB];
%% shuffling process 
n_shuffle=100;
ctxprefpre_shuffle=zeros(size(sigpre,1),n_shuffle);
for x=1:n_shuffle
    randshift=randi(500,1);
    sigshift=circshift(sigpre,randshift,2);
    randsplit=randperm(splitsize);
    sigsplit=reshape(sigshift,size(sigshift,1),[],splitsize);
    sigshuffle=reshape(sigsplit(:,:,randsplit),size(sigshift,1),[]);
    frb=calcmethod(sigshuffle(:,1:ctx_dur));
    fra=calcmethod(sigshuffle(:,ctx_dur+1:end));
    ctxprefpre_shuffle(:,x)=(fra-frb)./(fra+frb);
end
%% 
ctxsigpref=zeros(size(ctxprefpre_shuffle,1),1);
ts = [-2.56 2.56];      % T-Score
for n=1:size(ctxprefpre_shuffle,1)
    ctxpreftemp=ctxprefpre_shuffle(n,:);
    SEM = std(ctxpreftemp); %/sqrt(length(ctxpreftemp));               % Standard Error
    
    CI = mean(ctxpreftemp) + ts*SEM;
    if ctxpref_pre(n)>=CI(2)
        ctxsigpref(n)=1;
    elseif ctxpref_pre(n)<=CI(1)
        ctxsigpref(n)=-1;
    else
        ctxsigpref(n)=0;
    end
end
ctxA_ind=find(ctxsigpref==1);
ctxB_ind=find(ctxsigpref==-1);
figure;
hold on;
plot_pair(ctxpref_pre,ctxpref_post,ctxA_ind,'r');
plot_pair(ctxpref_pre,ctxpref_post,ctxB_ind,'b');
[~,p1]=ttest(ctxpref_pre(ctxA_ind),ctxpref_post(ctxA_ind));
[~,p2]=ttest(ctxpref_pre(ctxB_ind),ctxpref_post(ctxB_ind));
xlim([0 3]);
ax=gca;
ax.XTick=[1,2];
ax.XTickLabel={'Pre','Post'};
title('Context preference calculated by Firing Rate');
%%
figure;
subplot(121)
imagesc(zscore([sigpre(ctxA_ind,:); sigpre(ctxB_ind,:)],[],2),[1.65,6]);
subplot(122)
imagesc(zscore([sigpost(ctxA_ind,:); sigpost(ctxB_ind,:)],[],2),[1.65,6]);

%% FR 
sig2=ms.sigdeconvolved';
ids=SR_id_increase;
%ids=1:size(sig2,1);
FRpreA=calcFR(sig2(ids,session_start(strcmp(protocol,'preA')):session_end(strcmp(protocol,'preA'))));
FRpostA=calcFR(sig2(ids,session_start(strcmp(protocol,'postA')):session_end(strcmp(protocol,'postA'))));
FRpreB=calcFR(sig2(ids,session_start(strcmp(protocol,'preB')):session_end(strcmp(protocol,'preB'))));
FRpostB=calcFR(sig2(ids,session_start(strcmp(protocol,'postB')):session_end(strcmp(protocol,'postB'))));
Y=[FRpreB FRpreA FRpostA FRpostB];
ctxpref_pre=(FRpreA-FRpreB)./(FRpreA+FRpreB);
ctxpref_post=(FRpostA-FRpostB)./(FRpostA+FRpostB);

ctxA_ind=find(ctxpref_pre>=0.1);
ctxB_ind=find(ctxpref_pre<=-0.1);

figure;
hold on;
plot_pair(ctxpref_pre,ctxpref_post,ctxA_ind,'r');
plot_pair(ctxpref_pre,ctxpref_post,ctxB_ind,'b');
[~,p1]=ttest(ctxpref_pre(ctxA_ind),ctxpref_post(ctxA_ind));
[~,p2]=ttest(ctxpref_pre(ctxB_ind),ctxpref_post(ctxB_ind));
xlim([0 3]);
ax=gca;
ax.XTick=[1,2];
ax.XTickLabel={'Pre','Post'};
title('Context preference calculated by Firing Rate');

%% AUC
sig2=ms.sigraw';
ids=SR_id_increase;
%ids=1:size(sig2,1);
ctx_dur=2900;
shocksession_dur=3500;
sigpreA=sig2(ids,session_start(strcmp(protocol,'preA')):session_start(strcmp(protocol,'preA'))+ctx_dur-1);
sigpostA=sig2(ids,session_start(strcmp(protocol,'postA')):session_start(strcmp(protocol,'postA'))+ctx_dur-1);
sigshock=sig2(ids,session_start(strcmp(protocol,'conditioning')):session_start(strcmp(protocol,'conditioning'))+shocksession_dur-1);
sigpreB=sig2(ids,session_start(strcmp(protocol,'preB')):session_start(strcmp(protocol,'preB'))+ctx_dur-1);
sigpostB=sig2(ids,session_start(strcmp(protocol,'postB')):session_start(strcmp(protocol,'postB'))+ctx_dur-1);
FRpreA=calcAUC(sig2(ids,session_start(strcmp(protocol,'preA')):session_start(strcmp(protocol,'preA'))+ctx_dur-1));
FRpostA=calcAUC(sig2(ids,session_start(strcmp(protocol,'postA')):session_start(strcmp(protocol,'postA'))+ctx_dur-1));
FRpreB=calcAUC(sig2(ids,session_start(strcmp(protocol,'preB')):session_start(strcmp(protocol,'preB'))+ctx_dur-1));
FRpostB=calcAUC(sig2(ids,session_start(strcmp(protocol,'postB')):session_start(strcmp(protocol,'postB'))+ctx_dur-1));

ctxpref_pre=(FRpreA-FRpreB)./(FRpreA+FRpreB);
ctxpref_post=(FRpostA-FRpostB)./(FRpostA+FRpostB);
ctxA_ind=find(ctxpref_pre>=0.1);
ctxB_ind=find(ctxpref_pre<=-0.1);

figure;
hold on;
plot_pair(ctxpref_pre,ctxpref_post,ctxA_ind,'r');
plot_pair(ctxpref_pre,ctxpref_post,ctxB_ind,'b');
[~,p1]=ttest(ctxpref_pre(ctxA_ind),ctxpref_post(ctxA_ind));
[~,p2]=ttest(ctxpref_pre(ctxB_ind),ctxpref_post(ctxB_ind));
xlim([0 3]);
ax=gca;
ax.XTick=[1,2];
ax.XTickLabel={'Pre','Post'};
title('Context preference calculated by AUC');
%% Significant transient frequency 
sig2=ms.sigraw';
ids=SR_id_increase;
%ids=1:size(sig2,1);
ctx_dur=2900;
shocksession_dur=3500;
sigpreA=sig2(ids,session_start(strcmp(protocol,'preA')):session_start(strcmp(protocol,'preA'))+ctx_dur-1);
sigpostA=sig2(ids,session_start(strcmp(protocol,'postA')):session_start(strcmp(protocol,'postA'))+ctx_dur-1);
sigshock=sig2(ids,session_start(strcmp(protocol,'conditioning')):session_start(strcmp(protocol,'conditioning'))+shocksession_dur-1);
sigpreB=sig2(ids,session_start(strcmp(protocol,'preB')):session_start(strcmp(protocol,'preB'))+ctx_dur-1);
sigpostB=sig2(ids,session_start(strcmp(protocol,'postB')):session_start(strcmp(protocol,'postB'))+ctx_dur-1);
FRpreA=calcTF(sigpreA);
FRpostA=calcTF(sigpostA);
FRpreB=calcTF(sigpreB);
FRpostB=calcTF(sigpostB);

ctxpref_pre=(FRpreA-FRpreB)./(FRpreA+FRpreB);
ctxpref_post=(FRpostA-FRpostB)./(FRpostA+FRpostB);
ctxA_ind=find(ctxpref_pre>=0.1);
ctxB_ind=find(ctxpref_pre<=-0.1);

figure;
hold on;
plot_pair(ctxpref_pre,ctxpref_post,ctxA_ind,'r');
plot_pair(ctxpref_pre,ctxpref_post,ctxB_ind,'b');
[~,p1]=ttest(ctxpref_pre(ctxA_ind),ctxpref_post(ctxA_ind));
[~,p2]=ttest(ctxpref_pre(ctxB_ind),ctxpref_post(ctxB_ind));
xlim([0 3]);
ax=gca;
ax.XTick=[1,2];
ax.XTickLabel={'Pre','Post'};
title('Context preference calculated by Significant Transient Frequency');

%% Visualize 
[diffval,diffind]=sort(ctxpref_post-ctxpref_pre);
sigtemp=[sigpreA, sigpreB, sigshock, sigpostA, sigpostB];
figure;
subplot(121)
hold on;
imagesc(zscore(sigtemp(diffind,:),[],2),[1.65,3]);
session_start2=[1;ctx_dur+1;ctx_dur*2+1;ctx_dur*2+shocksession_dur+1;ctx_dur*3+shocksession_dur+1;ctx_dur*4+shocksession_dur+1];
axis tight;
ax=gca;
for i=1:length(session_start2)
    plot([session_start2(i) session_start2(i)],ax.YLim,'--k');
end
for s=1:numel(frame_shock)
    area([session_start2(3)+frame_shock(s); session_start2(3)+frame_shock(s)+shock_dur],[ylim; ylim],'FaceAlpha',0.4,'FaceColor','r','LineStyle','none')
end

subplot(122)
scatter(ctxpref_pre(ctxA_ind),ctxpref_post(ctxA_ind),'r');
hold on;
scatter(ctxpref_pre(ctxB_ind),ctxpref_post(ctxB_ind),'b');
xlim([-1,1]);
ylim([-1,1]);
plot([0 0],[-1 1],'--k')
plot([-1 1],[0 0],'--k')
axis square
xlabel('Pre-conditioning');
ylabel('Post-conditioning');



%%
function plot_pair(ctxpref_pre,ctxpref_post,id,color)
    for n=1:length(id)
        plot([1,2],[ctxpref_pre(id(n)),ctxpref_post(id(n))],'color',color);
    end
    scatter(ones(length(id),1),ctxpref_pre(id),'filled',color);
    
    scatter(2*ones(length(id),1),ctxpref_post(id),'filled',color);

end








