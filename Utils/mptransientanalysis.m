function [transient,threshold]=mptransientanalysis(sigt,isfig)

bsline=median(sigt);
stdsigt=mad(sigt);
threshold=bsline+1.65*stdsigt;



[pks,locp]=findpeaks(sigt);
[trf,loct]=findpeaks(-sigt);
if length(trf)<length(pks)
    trf=[trf sigt(end)];
    loct=[loct length(sigt)];
end
ind=find(pks>threshold);
pks=pks(ind);
locp=locp(ind);

if loct(ind(1))>locp(1)
    try
        saddlel=loct(ind-1);
        trf=trf(ind-1);
    catch %case where the signal started with a high peak
        pks(1)=[];
        locp(1)=[];
        ind(1)=[];
        saddlel=loct(ind-1);
        trf=trf(ind-1);
    end
else
    saddlel=loct(ind);
    trf=trf(ind);
end

% figure(5);
% clf
% plot(sigt);
% hold on;
% scatter(locp,pks);
% scatter(saddlel,-trf);
%%

%find the left and right bound of transients detected
idxl=zeros(length(pks),1);
idxr=zeros(length(pks),1);
t=1;
idxl(t)=locp(t);
while idxl(t)>1 && sigt(idxl(t))>threshold
    idxl(t)=idxl(t)-1;
end
idxr(t)=locp(t);
try
    while idxr(t)<saddlel(t+1) && sigt(idxr(t))>threshold
        idxr(t)=idxr(t)+1;
    end
    for t=2:(length(pks)-1)
        idxl(t)=locp(t);
        while idxl(t)>saddlel(t) && sigt(idxl(t))>threshold
            idxl(t)=idxl(t)-1;
        end
        
        idxr(t)=locp(t);
        while idxr(t)<saddlel(t+1) && sigt(idxr(t))>threshold
            idxr(t)=idxr(t)+1;
        end
        
    end
    
    t=length(pks); %special case of the last transient
    idxl(t)=locp(t);
    while idxl(t)>=idxl(t-1) && sigt(idxl(t))>threshold
        idxl(t)=idxl(t)-1;
    end
    idxr(t)=locp(t);
    while idxr(t)<=length(sigt) && sigt(idxr(t))>threshold
        if idxr(t)==length(sigt)
            break;
        else
            idxr(t)=idxr(t)+1;
        end
        
    end
    
catch %only 1 transient
    while idxr(t)<length(sigt) && sigt(idxr(t))>threshold
        idxr(t)=idxr(t)+1;
    end
end




temp=ismember(idxl,idxr);
if ismember(1,temp)==0 %no multipeak transients
    transient.peaks=pks;
    transient.locp=locp;
    transient.onset=idxl;
    transient.offset=idxr;
    transient.dur=idxr-idxl;
    if isfig==1
        %plot left and right bounds of all transients
        
        plot(sigt);
        hold on;
        line(xlim,[threshold threshold])
        scatter(transient.locp,transient.peaks,'*g') %plot all single peak transients
        scatter(transient.onset,sigt(transient.onset),'filled','b') %plot single peak transients onset
        scatter(transient.offset,sigt(transient.offset),'filled','y')
    end
else %there is multi-peak transients
    overlapind=find(ismember(idxl,idxr)==1); %remove these ind you treat multi-peaked transient as single peak
    mpindd(:,1)=find(diff(temp)==1);
    indd=find(diff(temp)==-1);
    if length(indd)~=size(mpindd,1)
        mpindd(:,2)=[indd;length(temp)];
    else
        mpindd(:,2)=indd;
    end
    
    spind=1:length(pks);
    % isolate single peak transients
    for i=1:size(mpindd,1)
        spind(mpindd(i,1):mpindd(i,2))=0;
    end
    spind(spind==0)=[];
    sptransient.peaks=pks(spind);
    sptransient.onset=idxl(spind);
    sptransient.offset=idxr(spind);
    sptransient.dur=sptransient.offset-sptransient.onset;
    sptransient.locp=locp(spind);
    
    % isolate multi-peak transients
    for i=1:size(mpindd,1)
        peaktemp=pks(mpindd(i,1):mpindd(i,2));
        [mptransient.peaks(i),loctemp]=max(peaktemp);
        mptransient.onset(i)=idxl(mpindd(i,1));
        mptransient.offset(i)=idxr(mpindd(i,2));
        mptransient.locp(i)=locp(mpindd(i,1)+loctemp-1);
    end
    mptransient.dur=mptransient.offset-mptransient.onset;
    
    
    transient.peaks=[sptransient.peaks'; mptransient.peaks'];
    transient.onset=[sptransient.onset; mptransient.onset'];
    transient.offset=[sptransient.offset; mptransient.offset'];
    transient.dur=[sptransient.dur; mptransient.dur'];
    transient.locp=[sptransient.locp';mptransient.locp'];
    if isfig==1
        %plot left and right bounds of all transients
        
        plot(sigt);
        hold on;
        for i=1:size(mpindd,1)
            mpinv=idxl(mpindd(i,1)):idxr(mpindd(i,2));
            plot(mpinv,sigt(mpinv),'r');
        end
        line(xlim,[threshold threshold])
        scatter(sptransient.locp,sptransient.peaks,'*g') %plot all single peak transients
        scatter(sptransient.onset,sigt(sptransient.onset),'.g') %plot single peak transients onset
        scatter(sptransient.offset,sigt(sptransient.offset),'.g')
        scatter(mptransient.locp,mptransient.peaks,'*r')
        scatter(mptransient.onset,sigt(mptransient.onset),'.r')
        scatter(mptransient.offset,sigt(mptransient.offset),'.r')
    end
end

end