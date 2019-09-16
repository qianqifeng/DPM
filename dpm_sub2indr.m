function [ind col]  = dpm_sub2indr(sze,vl,vu,dim)

ind0 = [0:sze(dim)-1].*sze(1);
ind = reshape([(vl+ind0); (vu+ind0)],1,numel(ind0)*2);
indstr = num2str(ind');
ind = eval(['[' reshape([indstr repmat(': ',1,length(vl))']',1,numel([indstr repmat(': ',1,length(vl))']')) ']']);
col = ceil((ind)/sze(1));
