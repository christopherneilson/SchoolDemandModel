function expdeltanorm=norm_delta(expdelta,Index)


if isstruct(Index)==0
    expdeltanorm=expdelta./expdelta(Index);
else
    expdeltanorm=NaN(size(expdelta));
    for m=reshape(Index.Markets,1,length(Index.Markets))
        for y=reshape(Index.Years,1,length(Index.Years))
            spots=(Index.SIndex(:,1)==m & Index.SIndex(:,2)==y);
            base=expdelta(Index.normID==1 & spots==1);
            
            if sum(length(base)~=1)
                error('No base school to norm')
            end
            expdeltanorm(spots)=expdelta(spots)/base;
        end
    end   
end
