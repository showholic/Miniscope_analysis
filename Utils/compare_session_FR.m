function compare_session_FR(ids,sig2,session_start,session_end,protocol)
    FRpreA=calcFR(sig2(ids,session_start(strcmp(protocol,'preA')):session_end(strcmp(protocol,'preA'))));
    FRpostA=calcFR(sig2(ids,session_start(strcmp(protocol,'postA')):session_end(strcmp(protocol,'postA'))));
    FRpreB=calcFR(sig2(ids,session_start(strcmp(protocol,'preB')):session_end(strcmp(protocol,'preB'))));
    FRpostB=calcFR(sig2(ids,session_start(strcmp(protocol,'postB')):session_end(strcmp(protocol,'postB'))));
    Y=[FRpreB FRpreA FRpostA FRpostB];
    g={'preB','preA','postA','postB'};
    %%% 1-way ANOVA
     
    % [~,~,stats]=anova1(Y,g);
    % [c,~,~,gnames] = multcompare(stats);

    % paired T-test
    [~,p1]=ttest(FRpreA,FRpostA);
    [~,p2]=ttest(FRpreB,FRpostB);
    [~,p3]=ttest(FRpreA,FRpreB);
    [~,p4]=ttest(FRpostA,FRpostB);

    % [p1,h1]=ranksum(FRpreA,FRpostA);
    % [p2,h2]=ranksum(FRpreB,FRpostB);
    % [p3,h3]=ranksum(FRpreA,FRpreB);
    % [p4,h4]=ranksum(FRpostA,FRpostB);
    group_pair={[2,3],[1,4],[2,1],[3,4]};
    p_pair=[p1 p2 p3 p4];

    figure
    boxplot(Y,'Notch','on','Labels',g);
    hold on;
    sigstar(group_pair((p_pair<=0.05)),p_pair(p_pair<=0.05));
    ax=gca;
    ax.YLim=[0 0.6];
end