function str = dblStructToString(stct)
    str = '';
    flds = fieldnames(stct);
    for ii=1:length(flds)
        val = getfield(stct, flds{ii});
        if isscalar(val)
            valstr = num2str(val, '%16.16g');
        else
            valstr = mat2str(val, 16);
        end
        str = [str '#' flds{ii} '=' valstr];
    end
    if str(1) == '#'
        str = str(2:end);
    end
end