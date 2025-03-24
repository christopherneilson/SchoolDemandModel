function vNames=format_vNames(nameVar,All)

if nargin==1
    All=0;
end

    nameVar=strvcat(nameVar);
   
    % Fix names to have nice caps
    vNames={};
    for ij=1:size(nameVar,1)
        
    txt=deblank(lower(nameVar(ij,:)));
    %txt=txt{:};
    loc=strfind(txt,'_');
    txt(loc)=' ';
    loc=strfind(txt,' ');
    txt(1)=upper(txt(1));
    for ii=1:length(loc)
        txt(loc(ii)+1);
    txt(loc(ii)+1)=upper(txt(loc(ii)+1));
    end
    txt=strtrim(txt);
    vNames{ij,1}=txt;
    end
   
    
    if All==1
       for ij=1:size(vNames,1) 
       txt=vNames{ij,1};
       dropSpace = isspace(txt); % Find within string blanks
       txt(dropSpace==1)=[];% Get rid of all blanks
       
       r=findstr('.',txt);
       if isempty(r)==0
       txt(r)=[];    
       end
       vNames{ij,1}=txt;
       end
    end