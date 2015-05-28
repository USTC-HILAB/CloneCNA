function corrected_ratio = CloneCNA_ratio_correction(data_ratio, data_gc)
% 30/10/2014 by yzh
% GC correction for read count ratio

corrected_ratio = data_ratio;
m_all_gc = median(data_ratio);
int_gc = floor(data_gc*100);
int_gc_u = unique(int_gc);

for i = 1:length(int_gc_u)
    tv = int_gc == int_gc_u(i);
    m_gc = median(data_ratio(tv));
    if m_gc > 0
        corrected_ratio(tv) = data_ratio(tv)*m_all_gc/m_gc;
    end
end

% corrected_ratio = data_ratio;
% m_all = median(data_ratio);
% int_gc = floor(data_gc*100);
% int_gc_u = unique(int_gc);
% int_map = floor(data_map*100);
% int_map_u = unique(int_map);
% 
% for i = 1:length(int_gc_u)
%     for j = 1:length(int_map_u)
%         tv = int_gc == int_gc_u(i) & int_map == int_map_u(j);
%         m_local = median(data_ratio(tv));
%         if m_local > 0
%             corrected_ratio(tv) = data_ratio(tv)*m_all/m_local;
%         end
%     end
% end