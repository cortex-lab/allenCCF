

%% convert annotation volume to index version

avI = av;
ind = st.index;
id = st.id;
idS = sparse(double(id), ones(size(id)), double(ind), double(max(id)),1);
avI(av==0) = 997;
tic; avS = idS(avI); toc
avS = full(avS);
avI = reshape(avS, size(av));
avI = uint16(avI+1);

% writeNPY(avI, 'annotation_volume_10um_by_index.npy')

%% make colormap from structure_tree based on index

q = st.color_hex_triplet;

q(cellfun(@numel,q)==5) = {'019399'}; % special case where leading zero was evidently dropped
c1 = cellfun(@(x)hex2dec(x(1:2)), q, 'uni', false);
c2 = cellfun(@(x)hex2dec(x(3:4)), q, 'uni', false);
c3 = cellfun(@(x)hex2dec(x(5:6)), q, 'uni', false);
cmap = horzcat(vertcat(c1{:}),vertcat(c2{:}),vertcat(c3{:}))./255;

% save allen_ccf_colormap.mat cmap