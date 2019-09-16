function ndx        = dpm_sub2ind(siz,varargin)

siz = double(siz);

if length(siz) ~= nargin-1
    %Adjust input
    if length(siz)<nargin-1
        %Adjust for trailing singleton dimensions
        siz = [siz ones(1,nargin-length(siz)-1)];
    else
        %Adjust for linear indexing on last element
        siz = [siz(1:nargin-2) prod(siz(nargin-1:end))];
    end
end

%Compute linear indices
k = [1 cumprod(siz(1:end-1))];
ndx = 1;
for i = 1:length(siz),
    v = varargin{i};
    ndx = ndx + (v-1)*k(i);
end
