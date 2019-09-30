function plot_pair(ctxpref_pre,ctxpref_post,id,color)
    for n=1:length(id)
        plot([1,2],[ctxpref_pre(id(n)),ctxpref_post(id(n))],'color',color);
    end
    scatter(ones(length(id),1),ctxpref_pre(id),'filled',color);
    
    scatter(2*ones(length(id),1),ctxpref_post(id),'filled',color);

end