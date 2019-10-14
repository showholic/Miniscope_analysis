clear;
addpath('.\Utils');
filepath='E:\Miniscope_Chenhaoshan\all_animal';
filenames=dir('E:\Miniscope_Chenhaoshan\all_animal\processed*.mat');
animal=cell(numel(filenames),1);
for n=1:numel(filenames)
    animal{n}=load(fullfile(filepath,filenames(n).name));
end
%% Speed correlation examples
animaldata=animal{1};

ns=4;
sig=animaldata.ms.sigraw';
sigtemp=zscore(sig(:,animaldata.session_start(ns):animaldata.session_end(ns)),[],2);
sigtemp(sigtemp<1.65)=0;
stemp=animaldata.behav{ns}.sf';
R=zeros(size(sigtemp,1),1);
for n=1:size(sigtemp,1)
    rtemp=corrcoef(sigtemp(n,:),stemp);
    R(n)=rtemp(1,2);
end
[Rsort,indsort]=sort(R,'descend');
% figure;
% hold on;
% sigt=sigtemp(indsort(1:10),:);
% sigt=sigt./max(sigt,[],2);
% plot((sigt+(1:size(sigt,1))')');
% yyaxis right;
% lh=plot(stemp);
% lh.Color=[0,1,0,0.4];
% ylim([0,200]);


speed_threshold=12;
FI_threshold=10;
binaryVector = stemp < speed_threshold;
[labeledVector, numRegions] = bwlabel(binaryVector);
% Measure lengths of each region and the indexes
measurements = regionprops(labeledVector, stemp, 'Area', 'PixelValues');
% Find regions where the area (length) are 3 or greater and
% put the values into a cell of a cell array
num_epoch=0;
i=1;
freeze_epoch={};

for k = 1 : numRegions
  if measurements(k).Area >= FI_threshold
    freeze_epoch{i}.value = measurements(k).PixelValues;
    freeze_epoch{i}.duration=length(freeze_epoch{i}.value);
    indtemp=find(labeledVector==k);
    freeze_epoch{i}.x1=indtemp(1);
    freeze_epoch{i}.x2=indtemp(end);
    num_epoch=num_epoch+1;
    i=i+1;
  end
end

figure;
hold on;

imagesc(sigtemp(indsort,:),[1.65 5]);
axis tight;
for i=1:length(freeze_epoch)
    area([freeze_epoch{i}.x1  freeze_epoch{i}.x2],[ylim; ylim],'FaceAlpha',0.4,'FaceColor','r','LineStyle','none')
end
yyaxis right
lh=plot(stemp);
plot(xlim,[speed_threshold speed_threshold],'--k');
lh.Color=[0 0 0 0.2];
%% LDA without locomotion
animaldata=animal{1};

sig=animaldata.ms.sigraw';
protocol=animaldata.protocol;
session_start=animaldata.session_start;
ms_ts=animaldata.ms.ms_ts;
sigz=zscore(sig,[],2);
ctx_dur=280;
sig_ctxBpre=sigz(:,session_start((strcmp(protocol,'preB')==1)):session_start((strcmp(protocol,'preB')==1))+ctx_dur*10-1);
sig_ctxApre=sigz(:,session_start((strcmp(protocol,'preA')==1)):session_start((strcmp(protocol,'preA')==1))+ctx_dur*10-1);
sig_ctxBpost=sigz(:,session_start((strcmp(protocol,'postB')==1)):session_start((strcmp(protocol,'postB')==1))+ctx_dur*10-1);
sig_ctxApost=sigz(:,session_start((strcmp(protocol,'postA')==1)):session_start((strcmp(protocol,'postA')==1))+ctx_dur*10-1);


siginterp_ctxApre=interpsig(sig_ctxApre,'preA',protocol,ms_ts,ctx_dur);
siginterp_ctxBpre=interpsig(sig_ctxBpre,'preB',protocol,ms_ts,ctx_dur);
siginterp_ctxApost=interpsig(sig_ctxApost,'postA',protocol,ms_ts,ctx_dur);
siginterp_ctxBpost=interpsig(sig_ctxBpost,'postB',protocol,ms_ts,ctx_dur);
pvd=[siginterp_ctxBpre siginterp_ctxApre]';
pvd_post=[siginterp_ctxBpost siginterp_ctxApost]';
ctx=zeros(size(pvd,1),1);
for i=1:ctx_dur
    ctx(i)=0;
end
for i=ctx_dur+1:size(pvd,1)
    ctx(i)=1;
end

W = LDA(pvd,ctx);

L = [ones(size(pvd,1),1) pvd] * W';
L2=[ones(size(pvd_post,1),1) pvd_post] * W';
P = exp(L2) ./ repmat(sum(exp(L2),2),[1 2]);
P1=P(:,1);
P2=P(:,2);

ci=0.8;
figure;
plot(L(1:ctx_dur,1),L(1:ctx_dur,2),'b');
hold on;
plot(L(ctx_dur:end,1),L(ctx_dur:end,2),'r');

plot(L2(1:ctx_dur,1),L2(1:ctx_dur,2),'b');
plot(L2(ctx_dur:end,1),L2(ctx_dur:end,2),'r');
axis off
%%
function siginterp=interpsig(sigt,trialname,protocol,ms_ts,ctx_dur)

    ms_ctx=double(ms_ts{(strcmp(protocol,trialname)==1)}(1:ctx_dur*10))/1000;
    siginterp=zeros(size(sigt,1),ctx_dur);
    ts_new=1:1:ctx_dur;
    for n=1:size(sigt,1)
        siginterp(n,:)=interp1(ms_ctx,sigt(n,:),ts_new);
    end
end