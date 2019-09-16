function sze        = get_size(current_grd)%得到curr grd中的Xcell中的大小
sze = [];
for i=1:length(current_grd.X)
    sze = [sze length(current_grd.X{i})];
end
if length(sze)==1
    sze = [sze 1];
end
