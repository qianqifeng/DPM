function str       	= dpm_code(s,v,m)
% v = [1 2 3];
% s = 'inp.U{*}';
% str = '';
if ~exist('m','var')
    m='';
end
v    = reshape(v,1,numel(v));
str  = repmat([s m],1,length(v));
str  = str(1:end-length(m));
star = strfind(str, '#');
str(star) = strrep(int2str(v), ' ', '');
