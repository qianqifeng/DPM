function c          = dpm_sizecmp(a,b)%�ȽϾ���a�;���b�ߴ��С�Ƿ���ͬ����ͬ��������ͬ�����ŷ���1
sa = size(a);
sb = size(b);
c = numel(sa)==numel(sb) & numel(a) == numel(b) & sum(sa==sb);
