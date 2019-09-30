%% LDA with Locomotion
accuracy=fld_withmotion(animal,30);
figure;
boxplot(accuracy([3 5 6],3:5),'Labels',{'A0','A','B'});
ylim([0 1]);
title('Decode accuracy');

%%
function accuracy=fld_withmotion(animal,speed_threshold)

    ctx_durall=[280;280;160;280;280];
    ctx_order={'preB','preA','conditioning','postB','postA'};
    ctx_id=[0 1 1 0 1];
    accuracy=zeros(numel(animal),length(ctx_order)-2);
    for nn=1:numel(animal)
        animaldata=animal{nn};
        protocol=animaldata.protocol;
        sig=animaldata.ms.sigraw';
        sigz=zscore(sig,[],2);
        
        
        sigepoch=cell(length(ctx_order),1);
        
        ctx_pre=[];
        ctx_post=[];
        for c=1:length(ctx_order)
            ctx_dur=ctx_durall(c);
            ns=find(strcmp(protocol,ctx_order{c})==1);
            sigtemp=sigz(:,animaldata.session_start(ns):animaldata.session_start(ns)+ctx_dur*10-1);

            stemp=animaldata.behav{ns}.sf(1:ctx_dur*10)';
            
            FI_threshold=10;
            binaryVector = stemp >= speed_threshold;
            [labeledVector, numRegions] = bwlabel(binaryVector);
            measurements = regionprops(labeledVector, stemp, 'Area', 'PixelValues');
            num_epoch=0;
            i=1;
            freeze_epoch={};
            sigepoch{c}=[];
            binsize=10;
            for k = 1 : numRegions
              if measurements(k).Area >= FI_threshold
                indtemp=find(labeledVector==k);
                epoch_duration=length(indtemp);
                ttemp=1:epoch_duration;
                edges=0:(binsize):ttemp(end)+1;
                [~,edges,bin]=histcounts(ttemp,edges);
                sigepocht=sigtemp(:,indtemp(1):indtemp(end));
                for m=1:length(edges)-1
                    sigepoch{c}=[sigepoch{c} mean(sigepocht(:,bin==m),2)];
                end
                i=i+1;
              end
            end
            if c==1 || c==2
                ctx_pre=[ctx_pre;ones(size(sigepoch{c},2),1)*ctx_id(c)];
            else
                ctx_post=[ctx_post;ones(size(sigepoch{c},2),1)*ctx_id(c)];
            end
        end

        pvd=[sigepoch{1} sigepoch{2}]';
        %pvd_post=[sigepoch{3} sigepoch{4}]';
        % 
        % W = LDA(pvd,ctx_pre);
        % 
        % L = [ones(size(pvd,1),1) pvd] * W';
        % L2=[ones(size(pvd_post,1),1) pvd_post] * W';
        % P = exp(L2) ./ repmat(sum(exp(L2),2),[1 2]);
        % P1=P(:,1);
        % P2=P(:,2);
        % 
        % 
        % figure;
        % plot(L(1:size(sigepoch{1},2),1),L(1:size(sigepoch{1},2),2),'b');
        % hold on;
        % plot(L(size(sigepoch{1},2):end,1),L(size(sigepoch{1},2):end,2),'r');
        % 
        % plot(L2(1:size(sigepoch{3},2),1),L2(1:size(sigepoch{3},2),2),'b');
        % plot(L2(size(sigepoch{3},2):end,1),L2(size(sigepoch{3},2):end,2),'r');
        MdLinear=fitcdiscr(pvd,ctx_pre);
        
        for t=3:length(ctx_order)
            predict_post=zeros(size(sigepoch{t},2),1);
            for n=1:size(sigepoch{t},2)
                prediction=predict(MdLinear,sigepoch{t}(:,n)');
                if prediction==ctx_id(t)
                    predict_post(n)=1;
                else
                    predict_post(n)=0;
                end
            end
            accuracy(nn,t)=length(find(predict_post==1))/size(sigepoch{t},2);
        end


       
    end
end
%%
