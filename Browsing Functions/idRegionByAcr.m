

function tree = idRegionByAcr(st, acr)

id = st.id(find(strcmp(st.acronym, acr)));

r = st(find(strcmp(st.acronym, acr)),:);

fprintf(1, '%s\n', r.name{1})

id = r.id; 
tree = {acr}; 

while id~=997
    parentID = r.parent_structure_id;
    
    r = st(st.id==parentID,:); 
    
    for q = 1:numel(tree)
        fprintf(1, '-'); 
    end
    
    fprintf(1, '%s\n', r.name{1})
    
    id = r.id;
    tree{end+1} = r.acronym{1};
end