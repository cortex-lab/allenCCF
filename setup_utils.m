%% downloads

% data from: http://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/

% nrrd reader from: https://uk.mathworks.com/matlabcentral/fileexchange/34653-nrrd-format-file-reader

% structure tree from:
% http://api.brain-map.org/api/v2/data/query.csv?criteria=model::Structure,rma::criteria,[ontology_id$eq1],rma::options[order$eq%27structures.graph_order%27][num_rows$eqall]

%% sanitize structure tree

sanitizeStructureTree('query.csv', 'structure_tree_safe_2017.csv');
st = loadStructureTree('structure_tree_safe_2017.csv');


%% convert annotation volume to index version
% This cell can take a while, both for loading and doing the conversion

% av = nrrdread('annotation_10.nrrd');
% st = loadStructureTree('structure_tree_safe.csv');

avI = av;
ind = st.index;
id = st.id;
idS = sparse(double(id), ones(size(id)), double(ind), double(max(id)),1);
avI(av==0) = 997;
tic; avS = idS(avI); toc
avS = full(avS);
avI = reshape(avS, size(av));
avI = uint16(avI+1);

% running this for the 2017 version, I have to also add this line:
avI = permute(avI, [2 1 3]);

% writeNPY(avI, 'annotation_volume_10um_by_index.npy')

%% make colormap from structure_tree based on index

q = st.color_hex_triplet;

q(cellfun(@numel,q)==5) = {'019399'}; % special case where leading zero was evidently dropped
c1 = cellfun(@(x)hex2dec(x(1:2)), q, 'uni', false);
c2 = cellfun(@(x)hex2dec(x(3:4)), q, 'uni', false);
c3 = cellfun(@(x)hex2dec(x(5:6)), q, 'uni', false);
cmap = horzcat(vertcat(c1{:}),vertcat(c2{:}),vertcat(c3{:}))./255;

% save allen_ccf_colormap.mat cmap