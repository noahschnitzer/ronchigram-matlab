function [columns,means] = column_means(data_x,data_y)
    data_x = round(data_x);
    [sorted_data_x,I] = sort(data_x);
    sorted_data_y = data_y(I);
    columns = [];
    means = [];
    columns(1) = sorted_data_x(1);
    means(1) = sorted_data_y(1);
    ct = 1;
    cv = columns(1);
    for it = 2:length(sorted_data_x)
       if sorted_data_x(it) ==  cv
           ct = ct + 1; 
           means(end) = means(end) + sorted_data_y(it);
       else
           means(end) = means(end)/ct;
           ct = 1;
           cv = sorted_data_x(it);
           means(end+1) = sorted_data_y(it);
           columns(end+1) = sorted_data_x(it);
       end
    end

end