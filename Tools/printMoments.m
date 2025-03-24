function M=printMoments(Model,MM)

 wMM=diag(Model.WMM);
 Mm=[];
 mnames={'Average Quality','Average Price','Average Distance'};MMnames={};
 for m=1:3
     
     if m==3 
     pick=(Model.mAllIndex(:,end-1)==3 & Model.mAllIndex(:,2)==2011);
     year=Model.mAllIndex(pick,2);
     pick(year<2011)=0;

     Mm=[Mm; [ Model.MM(pick) MM(pick)]];
     type=Model.mAllIndex(pick,3);
     Market=Model.mAllIndex(pick,1);
     
     MMnames=strvcat(MMnames,strvcat([repmat('Market - ',size(year,1),1) num2str(Market) repmat([' - '],size(year,1),1) repmat(mnames{m},size(year,1),1) repmat([' - '],size(year,1),1) num2str(year) repmat(' - ',size(year,1),1) num2str(type)]));    
         
     else
     pick=(Model.mAllIndex(:,end-1)==m);
     year=Model.mAllIndex(pick,2);
     Market=Model.mAllIndex(pick,1);
     Mm=[Mm; [Model.MM(pick) MM(pick)]];
     type=Model.mAllIndex(pick,3);
     MMnames=strvcat(MMnames,strvcat([repmat('Market - ',size(year,1),1) num2str(Market) repmat([' - '],size(year,1),1) repmat(mnames{m},size(year,1),1) repmat([' - '],size(year,1),1) num2str(year) repmat(' - ',size(year,1),1) num2str(type)]));
     end
 end
 vm.rnames=strvcat({' ', MMnames});
 vm.cnames=strvcat({' Data','Model'});
 mprint(Mm,vm)
 