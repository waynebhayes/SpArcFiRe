function stct = stringToDblStruct(str)
    stct = struct();
    elts = regexp(str, '#', 'split');
    for ii=1:length(elts)
        fv = regexp(elts{ii}, '=', 'split');
        assert(length(fv) == 2)
        fldname = fv{1};
        valstr = fv{2};
        assert(ischar(valstr))
        assert(~isempty(valstr));
        if valstr(1) == '['
            value = str2num(valstr);
        else
            value = str2double(fv{2});
        end
        stct = setfield(stct, fldname, value);
    end
end