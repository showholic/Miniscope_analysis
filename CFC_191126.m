clear;
close all;
filepath='D:\Miniscope_Chenhaoshan\Results_191126\20190921_170501';
load(fullfile(filepath,'ms.mat'));
load(fullfile(filepath,'ms_dff.mat'));
shockts=readtable(fullfile(filepath,'shock_behavts.csv'));
ms_start=[1;1;1;20;20;28];
shock_start=[5426;7227;9029];
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
%% Behavior analysis
load(fullfile(filepath,'arena.mat'));
body_part='Body';
behavpos_file=dir(fullfile(filepath,'behave_video*.csv'));
behavpos_filename=cell(length(behavpos_file),1);

behavdatatable=cell(length(behavpos_file),1);
for n=1:length(behavpos_file)
    behavpos_filename{n}=fullfile(filepath,behavpos_file(n).name);
    behavdatatable{n}=readtable(behavpos_filename{n},'ReadRowNames',true);    
end

datatemp=behavdatatable{1};
temp=datatemp{1,:};
tempind=zeros(length(temp),1);
for i=1:length(temp)
    tempind(i)=strcmp(temp{i},body_part);
end
behavdata=cell(length(behavpos_file),1);
body_ind=find(tempind==1);

figure;
for n=1:length(behavpos_file)
    datatemp=behavdatatable{n};
    behavdata{n}.c=str2double(datatemp{3:end,body_ind(3)});
    behavdata{n}.x=str2double(datatemp{3:end,body_ind(1)});
    behavdata{n}.y=str2double(datatemp{3:end,body_ind(2)});
    behavdata{n}.y(behavdata{n}.y<0)=1;
    behavdata{n}.x(behavdata{n}.x<0)=1;
    for i=1:length(behavdata{n}.x)
        if bw{n}(abs(round(behavdata{n}.y(i))),abs(round(behavdata{n}.x(i))))==0
            try
                behavdata{n}.x(i)=behavdata{n}.x(i-1);
                behavdata{n}.y(i)=behavdata{n}.y(i-1);
            catch
                behavdata{n}.x(i)=behavdata{n}.x(i+1);
                behavdata{n}.y(i)=behavdata{n}.y(i+1);
            end
        end
    end
    subplot(1,length(behavpos_file),n)
    plot(behavdata{n}.x,behavdata{n}.y);
    axis equal;
    axis tight;
end

%% Method 1: align ms by behavior timestamps 
behav2=cell(length(behavpos_file),1);
for ns=1:length(behavpos_file)
    behavts=(0:1/30:(length(behavdata{ns}.x)-1)*(1/30))';
    behavts=behavts-behavts(ms_start(ns));
    behavts(behavts<0)=[];
    mst=double(ms_ts{ns})'/1000;
    behavts(behavts>mst(end)+1/30)=[];
    behavts_ds=behavts(1:3:end);
    siginterp=zeros(size(sig,1),length(behavts_ds));
    spkinterp=zeros(size(sig,1),length(behavts_ds));
    xtemp=behavdata{ns}.x(ismember(behavts_ds,behavts));
    ytemp=behavdata{ns}.y(ismember(behavts_ds,behavts));
    stemp=[0;sqrt(diff(ytemp).^2+diff(xtemp).^2)./diff(behavts_ds)];
    windowSize = 5; 
    b = (1/windowSize)*ones(1,windowSize);
    a = 1;
    sfilt = filter(b,a,stemp);
    for i=1:size(sig,1)
        sigt=sig(i,session_start(ns):session_end(ns));
        spkt=ms_dff.S_dff(i,session_start(ns):session_end(ns));
        siginterp(i,:)=interp1(mst',sigt,behavts_ds');
        spkinterp(i,:)=interp1(mst',spkt,behavts_ds');
    end
    behav2{ns}.t=behavts_ds;
    behav2{ns}.x=xtemp;
    behav2{ns}.y=ytemp;
    behav2{ns}.s=stemp;
    behav2{ns}.sf=sfilt;
    behav2{ns}.sig=siginterp;
    behav2{ns}.spk=spkinterp;
    
end



%% Method 2: align behavior by ms timestamps 
figure;
behav=cell(length(behavpos_file),1);
for ns=1:length(behavpos_file)
    behavts=(0:1/30:(length(behavdata{ns}.x)-1)*(1/30))';
    behavts=behavts-behavts(ms_start(ns));
    mst=double(ms_ts{ns})'/1000;
    temp1=repmat(behavts,1,length(mst))';
    temp2=abs(temp1-mst);
    [minval,frame_behav]=min(temp2,[],2);
    ctemp=behavdata{ns}.c(frame_behav);
    ttemp=behavts(frame_behav);
    xtemp=behavdata{ns}.x(frame_behav);
    ytemp=behavdata{ns}.y(frame_behav);
    stemp=[0;sqrt(diff(ytemp).^2+diff(xtemp).^2)./diff(ttemp)];
    windowSize = 5; 
    b = (1/windowSize)*ones(1,windowSize);
    a = 1;
    sfilt = filter(b,a,stemp);
    behav{ns}.x=xtemp;
    behav{ns}.y=ytemp;
    behav{ns}.t=ttemp;
    behav{ns}.s=stemp;
    behav{ns}.sf=sfilt;
    subplot(length(behavpos_file),1,ns);
    sigtemp=sig(:,session_start(ns):session_end(ns));
    imagesc(zscore(sigtemp,[],2),[0 6]);
    hold on;
    yyaxis right;
    lh=plot(sfilt);
    lh.Color=[0,1,0,0.3];
    if ns==4
        for s=1:numel(frame_shock)
            area([frame_shock(s); frame_shock(s)+shock_dur],[ylim; ylim],'FaceAlpha',0.4,'FaceColor','r','LineStyle','none')
        end
    end
end

[~,filepart1]=fileparts(fileparts(filepath));
save(['D:\Miniscope_Chenhaoshan\all_animal\processed' filepart1(end-6:end)],'shockts','ms_start', ...
    'shock_start','protocol','ms','ms_dff','session_start','session_end','frame_shock','behav','behav2');
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
img=imshow(vid_shock(:,:,1),[0 130]);
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
%sig2=ms_dff.S_dff;
% All cells
compare_session_FR(1:size(sig2,1),sig2,session_start,session_end,protocol)
% Shock activated cells 
compare_session_FR(SR_id_increase,sig2,session_start,session_end,protocol)
% Shock decreased cells 
compare_session_FR(SR_id_decrease,sig2,session_start,session_end,protocol)
%% Draw Arena
behav_videofile=dir(fullfile(filepath,'behave_video*.mp4'));
figure;
bw=cell(length(behav_videofile),1);
for n=1:length(behav_videofile)
    v=VideoReader(fullfile(filepath,behav_videofile(n).name));
    frame=readFrame(v);
    imshow(frame);
    roih=drawpolygon;
    bw{n}=roipoly(frame,roih.Position(:,1),roih.Position(:,2));
end
save(fullfile(filepath,'arena.mat'),'bw');




