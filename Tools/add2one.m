function vectorOut=add2one(vector0)

vector0=reshape(vector0,length(vector0),1);

if sum(vector0)==1
    vectorOut=vector0;
elseif abs(sum(vector0)-1)>10^6
    error('Something Else Wrong')
else
    vector1=double(vector0);
    for i=1:100
        vector2=vector1./sum(vector1);
        dif=1-sum(vector2);
        vector2(1)=vector2(1)+dif;
        
      
        if sum(vector2)==1
            vectorOut=vector2;
            %fprintf('Fixed Imprecise Sum \n')
            break
        else
            vector1=vector2;
        end
    end
    
    
end
