function norm = normalize_data(ins,type)
    if strcmp(type,'max')
    norm = ins - min(ins(:));
    norm = norm ./ max(norm(:));
    elseif strcmp(type,'total')
       norm = ins./ sum(ins(:)); 
    else
        norm = -1;
        
    end
end