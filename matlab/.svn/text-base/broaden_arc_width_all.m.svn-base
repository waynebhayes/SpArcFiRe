function broaden_arc_width_all(dir_path, suffix, arcs_colored)

flist = dir(dir_path);
flist = flist(arrayfun(@(x)(~x.isdir), flist));
flist = {flist.name};
sfx_len = length(suffix);
flist = flist(cellfun(@(x)((length(x) >= sfx_len) && strcmp(x(end-(sfx_len-1):end), suffix)), flist));
img_names = flist;

for ii=1:length(img_names)
    in_name = img_names{ii}
    toks = regexp(in_name,'\.','split');
    out_name = [toks{1:end-1} '_arcEnhanced.png']
%     if exist(out_name, 'file')
%         error([out_name ' already exists'])
%     end
    img = imread([dir_path filesep in_name]);
    tol = 100;
    tol = 0;
    if arcs_colored
        img = broaden_arc_width(img, [255, 0, 0], [], tol);
        img = broaden_arc_width(img, [0, 255, 0], [], tol);
        img = broaden_arc_width(img, [0, 0, 255], [0, 255, 255], tol);
        img = broaden_arc_width(img, [255, 0, 255], [], tol);
    else
        img = broaden_arc_width(img, [255, 255, 255], [], tol);
    end
    imwrite(img, [dir_path filesep out_name]);
end

end