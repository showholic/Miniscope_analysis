clear;
addpath('.\Utils');
filepath='D:\Miniscope_Chenhaoshan\all_animal';
filenames=dir('D:\Miniscope_Chenhaoshan\all_animal\processed*.mat');
animal=cell(numel(filenames),1);
for n=1:numel(filenames)
    animal{n}=load(fullfile(filepath,filenames(n).name));
end
load(fullfile(filepath,'arena.mat'));
%% Speed correlation examples
animaldata=animal{2};
ns=6;
sig=animaldata.ms.sigraw';
sigtemp=sig(:,animaldata.session_start(ns):animaldata.session_end(ns));
stemp=animaldata.behav{ns}.sf';
R=zeros(size(sigtemp,1),1);
for n=1:size(sigtemp,1)
    rtemp=corrcoef(sigtemp(n,:),stemp);
    R(n)=rtemp(1,2);
end
[Rsort,indsort]=sort(R,'descend');
figure;
hold on;
sigt=sigtemp(indsort(1:10),:);
sigt=sigt./max(sigt,[],2);
plot((sigt+(1:size(sigt,1))')');
yyaxis right;
lh=plot(stemp);
lh.Color=[0,1,0,0.4];

ylim([0,200]);

%%
speed_threshold=10;
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

imagesc(sigtemp(indsort,:),[5 10]);
axis tight;
for i=1:length(freeze_epoch)
    area([freeze_epoch{i}.x1  freeze_epoch{i}.x2],[ylim; ylim],'FaceAlpha',0.4,'FaceColor','r','LineStyle','none')
end
yyaxis right
plot(stemp);
plot(xlim,[speed_threshold speed_threshold],'--k');

%%

sigma=3;
binsize=20;
bwtemp=bw{ns};
[yy,xx]=find(bwtemp==1);
arena.x1=min(xx);
arena.y1=min(yy);
arena.x2=max(xx);
arena.y2=max(yy);
trkxtemp=animaldata.behav{ns}.x;
trkytemp=animaldata.behav{ns}.y;
% figure;
% plot(trkxtemp,trkytemp);


sig=animaldata.ms.sigdeconvolved';
sigtemp=sig(:,animaldata.session_start(ns):animaldata.session_end(ns));
maxbatch=1;
figure('Name',['Session ' num2str(ns)],'NumberTitle','off');
ha = tight_subplot(maxbatch,10,[.01 .03],[.05 .01],[.01 .01]);
nbinx=ceil((arena.x2-arena.x1)/binsize);
nbiny=ceil((arena.y2-arena.y1)/binsize);
[occupancy,Xedges,Yedges,binX,binY]=histcounts2(trkxtemp,trkytemp,[nbinx nbiny]);
oi=occupancy/sum(occupancy(:));
for nbatch=1:maxbatch

    for i=1:10
        sigzs=zscore(sigtemp(((nbatch-1)*10+i),:));
        signm=double(sigtemp(((nbatch-1)*10+i),:));
        signm(signm<0)=0;
        sigmean=mean(signm);

        zsmap=zeros(length(Yedges)-1,length(Xedges)-1);
        simap=zeros(length(Yedges)-1,length(Xedges)-1);
        for y=1:size(zsmap,1)
            for x=1:size(zsmap,2)
                indt=find(binX==x & binY==y);
                if ~isempty(indt)
                    sigt=sigzs(indt);
                    zsmap(y,x)=mean(sigt);
                    ai=mean(signm(indt),'omitnan');
                    simap(y,x)=oi(y,x)*ai*log2(ai/sigmean);
                else
                    zsmap(y,x)=0;
                    simap(y,x)=0;
                end
            end
        end
        zsmapraw=zsmap;
        zsmap(isnan(zsmap))=0;
        zsmapfilt=imgaussfilt(flipud(zsmap),sigma);
        axes(ha((nbatch-1)*10+i));
        img=imagesc('XData',Xedges,'YData',Yedges,'CData',zsmapfilt);
        %img=imagesc('XData',Xedges,'YData',Yedges,'CData',zsmap);
        axis equal;
        axis tight;
        axis off;
        

    end
end

%% Find significant place cells
animaldata=animal{2};
ns=6;
binsize=20;
numshuffle=100;
bwtemp=bw{ns};
[yy,xx]=find(bwtemp==1);
arena.x1=min(xx);
arena.y1=min(yy);
arena.x2=max(xx);
arena.y2=max(yy);
ctx_dur=2800;
trkxtemp=animaldata.behav{ns}.x(1:ctx_dur);
trkytemp=animaldata.behav{ns}.y(1:ctx_dur);
% figure;
% plot(trkxtemp,trkytemp);

sig=animaldata.ms_dff.S_dff;
sigtemp=sig(:,animaldata.session_start(ns):animaldata.session_start(ns)+ctx_dur-1);
maxbatch=1;

nbinx=ceil((arena.x2-arena.x1)/binsize);
nbiny=ceil((arena.y2-arena.y1)/binsize);
[occupancy,Xedges,Yedges,binX,binY]=histcounts2(trkxtemp,trkytemp,[nbinx nbiny]);
oi=occupancy/sum(occupancy(:));
placecell={};
npc=1;
for i=1:size(sigtemp,1)
    sigzs=zscore(sigtemp(i,:));
    signm=double(sigtemp(i,:));
    signm(signm<0)=0;
    sigmean=mean(signm);

    zsmap=zeros(size(occupancy));
    simap=zeros(size(occupancy));
    for y=1:size(zsmap,1)
        for x=1:size(zsmap,2)
            indt=find(binX==x & binY==y);
            if ~isempty(indt)
                sigt=sigzs(indt);
                zsmap(y,x)=mean(sigt);
                ai=mean(signm(indt),'omitnan');
                simap(y,x)=oi(y,x)*ai*log2(ai/sigmean);
            else
                zsmap(y,x)=0;
                simap(y,x)=0;
            end
        end
    end
    zsmapraw=zsmap;
    zsmap(isnan(zsmap))=0;
    %zsmapfilt=imgaussfilt(flipud(zsmap),sigma);


    % Calculate spatial information 
    simap(isnan(simap))=0;
    si=sum(simap(:));
    sishuffle=zeros(numshuffle,1);
    for s=1:numshuffle
        sigshuffle=shuffle_sig(signm,500,10);
        sishuffle(s)=calcSI(sigshuffle,oi,binX,binY,Yedges,Xedges);        
    end
    
    if si>=1.65*std(sishuffle)+mean(sishuffle) %&& si>=0.1
        placecell{npc}.zsmap=zsmap;
        placecell{npc}.si=si;
        placecell{npc}.id=i;
        npc=npc+1;
    end
end


%% Display place cell place fields 
placecell2={};
j=1;
for i=1:length(placecell)
    if placecell{i}.si>=0.5
        placecell2{j}=placecell{i};
        j=j+1;
    end
end

sigma=2.5;
maxbatch=ceil(length(placecell2)/10);
figure('Name',['Session ' num2str(ns)],'NumberTitle','off');
ha = tight_subplot(maxbatch,10,[.01 .01],[.05 .05],[.01 .01]);
for i=1:length(placecell2)
    zsmapfilt=imgaussfilt(flipud(placecell2{i}.zsmap),sigma);
    imagesc(ha(i),zsmapfilt)
    %imagesc(placecell{i}.zsmap)
    axes(ha(i));
    axis square;
    axis off;
end
delete(ha(length(placecell2)+1:end));