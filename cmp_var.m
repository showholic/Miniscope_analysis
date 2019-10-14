clear
load('E:\Miniscope_Chenhaoshan\Results_191082\20190920_181558\ms_dff.mat');
load('E:\Miniscope_Chenhaoshan\Results_191082\20190920_181558\ms.mat');
%%
id=10;
sigdff=ms_dff.dff;
sig=ms.sigraw';
figure
lh1=plot(sig(id,:));
lh1.Color=[0.5,1,0.5,0.5];
lh1.LineWidth=2;
hold on;
lh2=plot(sigdff(id,:));



s_dff=ms_dff.S_dff;
s=ms.sigdeconvolved';
% figure
% lh1=plot(s(1,:));
% lh1.Color=[1,1,0.5,0.5];
% lh1.LineWidth=2;
% hold on;
% lh2=plot(s_dff(1,:));
% lh1.Color=[0.5,1,1,1];

%%
id=1;
figure
ax1=subplot(211)
lh1=plot(sig(id,:));
hold on;
lh2=plot(s(id,:));
legend({'C','S'});
axis tight;
ax2=subplot(212)
lh1=plot(sigdff(id,:));
hold on;
lh2=plot(s_dff(id,:));
legend({'df/f0','Sdff'});
axis tight;
linkaxes([ax1 ax2],'xy');

%%
AUC=zeros(size(sig,1),1);
FR=zeros(size(sig,1),1);
for i=1:size(sig,1)
    sigt=sig(i,:);
    sigt(sigt<0)=0;
    AUC(i)=trapz(sigt)/length(sigt);
    
    spkt=s(i,:);
    spkt(spkt>0)=1;
    FR(i)=sum(spkt)/length(spkt);
    
end

figure;
scatter(AUC,FR);
R1=corrcoef(AUC,FR);

%% 
AUC=zeros(size(sig,1),1);
FR=zeros(size(sig,1),1);
for i=1:size(sig,1)
    sigt=sigdff(i,:);
    sigt(sigt<0)=0;
    AUC(i)=trapz(sigt)/length(sigt);
    
    spkt=s_dff(i,:);
    spkt(spkt>0)=1;
    FR(i)=sum(spkt)/length(spkt);
    
end

figure;
scatter(AUC,FR);
R2=corrcoef(AUC,FR);

%% AUC with threshold cutoff vs FR
AUC=zeros(size(sig,1),1);
FR=zeros(size(sig,1),1);
for i=1:size(sig,1)
    sigt=sigdff(i,:);
    baseline=median(sigt)+3*mad(sigt);
    sigt(sigt<baseline)=0;
    AUC(i)=trapz(sigt)/length(sigt);
    
    spkt=s_dff(i,:);
    spkt(spkt>0)=1;
    FR(i)=sum(spkt)/length(spkt);
    
end

figure;
scatter(AUC,FR);
R3=corrcoef(AUC,FR);

%% AUC vs AUC
AUC=zeros(size(sig,1),1);
FR=zeros(size(sig,1),1);
for i=1:size(sig,1)
    sigt=sigdff(i,:);
    baseline=median(sigt)+3*mad(sigt);
    %baseline=0;
    sigt(sigt<baseline)=0;
    AUC(i)=trapz(sigt)/length(sigt);
    
    spkt=s_dff(i,:);
    FR(i)=trapz(spkt)/length(sigt);
    
end

figure;
scatter(AUC,FR);
R3=corrcoef(AUC,FR);