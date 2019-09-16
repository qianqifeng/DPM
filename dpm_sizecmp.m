function c          = dpm_sizecmp(a,b)%比较矩阵a和矩阵b尺寸大小是否相同，相同行数和相同列数才返回1
sa = size(a);
sb = size(b);
c = numel(sa)==numel(sb) & numel(a) == numel(b) & sum(sa==sb);
