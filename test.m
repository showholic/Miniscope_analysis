[file,path] = uigetfile('*.mat');
app.SelectRegionButton.Enable = 'on';

if isequal(file,0)
    %disp('User selected Cancel');   
    path='E:\Haoshan\Miniscope\191126\20190905\Results';
    app.rawdata=load('E:\Haoshan\Miniscope\191126\20190905\Results\ms.mat');
    try
        app.videoreg=VideoReader(fullfile(path,'ms_mc.avi'));
        app.videoraw=VideoReader(fullfile(path,'msCam_concat.avi'));
    catch
        disp('No video file found');
        app.DisplayButton.Enable='off';
        app.Slidervideo.Enable='off';
        app.NextButton.Enable='off';
        app.PreviousButton.Enable='off';

    end
    path='E:\miniscope\Data\192501_aCA3_VR\04272019';
else
    %disp(['User selected ', fullfile(path,file)]);
    app.rawdata=load(fullfile(path,file));       
    try
        app.videoreg=VideoReader(fullfile(path,'ms_mc.avi'));
        app.videoraw=VideoReader(fullfile(path,'msCam_concat.avi'));
    catch
        disp('No video file found');
        disp('No video file found');
        app.DisplayButton.Enable='off';
        app.Slidervideo.Enable='off';
        app.NextButton.Enable='off';
        app.PreviousButton.Enable='off';
    end
end

app.pathname=path;

try
    app.acceptedPool=app.rawdata.acceptedPool;
    app.deletePool=app.rawdata.deletePool;
catch
    app.acceptedPool=[];
    app.deletePool=[];
end            
app.processedPool=[app.deletePool;app.acceptedPool];
app.rawdata.pixw=app.rawdata.ms.width;
app.rawdata.pixh=app.rawdata.ms.height;
%%
app.img2=app.rawdata.ms.PNR;
idaccept=app.rawdata.ms.idx_accepted+1;
roifn=reshape(app.rawdata.ms.SFP,[],size(app.rawdata.ms.SFP,3));
roifn=roifn(:,idaccept);
sigraw=app.rawdata.ms.sigraw';
sigraw=sigraw(idaccept,:);
%%
% sigfn=app.rawdata.ms.sigdeconvolved';
% sigfn=sigfn(idaccept,:);
%%
figure;
plot(sigfn(1,:));
hold on;
plot(sigraw(1,:));
%%
app.videoreg.CurrentTime=0;
t=[];
while hasFrame(app.videoreg)
    t=[t, app.videoreg.CurrentTime];
    readFrame(app.videoreg);
end
%%
frameno=100;
app.videoraw.CurrentTime=round((1/app.videoraw.FrameRate),5)*(frameno-1);
img=reshape(roifn(:,:)*sigraw(:,frameno),app.rawdata.pixh,app.rawdata.pixw);
img(img<0)=0;
figure;
subplot(122)
imagesc(img)
subplot(121)
imagesc(readFrame(app.videoraw),[0,0.85])
%%
frameno=2;
app.videoreg.CurrentTime=round((1/app.videoreg.FrameRate),5)*(frameno-1);
img=reshape(roifn(:,:)*sigraw(:,frameno),app.rawdata.pixh,app.rawdata.pixw);
figure;
subplot(122)
imagesc(img,[0,50])
subplot(121)
imagesc(readFrame(app.videoreg),[0 0.95]);