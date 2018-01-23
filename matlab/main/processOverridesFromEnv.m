function stgs = processOverridesFromEnv(stgs)
% Modifies fields in the stgs struct according to system environment
% variables. A field is modified if the full capitalization of its name
% matches a system environment variable. The value of the field is set to
% the value of the environment variable, after doing a string-to-number
% conversion.

flds = fieldnames(stgs);
for fld_idx = 1:length(flds)
    cur_fld = flds{fld_idx};
    new_value = getenv(upper(cur_fld));
    new_value = str2num(new_value);
    if ~isempty(new_value)
        old_value = stgs.(cur_fld);
        if (length(size(new_value)) ~= length(size(old_value))) || ...
                any(size(new_value) ~= size(old_value))
            error('cannot change value of "%s" from "%s" to "%s" (size mismatch)',...
                cur_fld, num2str(old_value), num2str(new_value))
        end
        fprintf('applying environment-variable override to "%s": %s -> %s\n',...
            cur_fld, num2str(old_value), num2str(new_value))
        stgs.(cur_fld) = new_value;
    end
end

end