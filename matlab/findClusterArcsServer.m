function findClusterArcsServer(stgs)

server_run_id = sprintf('%12.f', round(now() * 1e10));
fprintf('server run ID is %s\n', server_run_id);

if strcmpi(stgs, 'NONE') || isempty(stgs)
    fprintf('No fit settings specified. Using defaults.\n');
    stgs = getDefaultSettings();
    stgs = processOverridesFromEnv(stgs);
elseif ischar(stgs)
    fprintf('stgs specified as a path string, attempting to load from disk.\n');
    load(stgs, 'stgs');
end
stgs_txtfile_name = sprintf('%s_settings.txt', server_run_id);
diary(stgs_txtfile_name);
fprintf('%s\n', datestr(clock));
fprintf('server run ID: %s\n', server_run_id);
stgs;
diary off
fprintf('settings recorded to %s\n', stgs_txtfile_name);

expected_tokens = 7;
nxtSeqNum = 1;
consecutive_empty_lines = 0;
max_consecutive_empty_lines = 1;
while true
    fprintf('ready for next image (in_dir, image_name, image_suffix, starmask_suffix,guide_dir, output_dir):\n');
    in_str = input(sprintf('%d>> ', nxtSeqNum), 's');
    %in_str = input('', 's');
    %in_str
    %fprintf('"%s"\n', in_str);
    nxtSeqNum = nxtSeqNum + 1;
    if isempty(in_str)
        if ~ischar(in_str)
            consecutive_empty_lines = consecutive_empty_lines + 1;
            if consecutive_empty_lines >= max_consecutive_empty_lines
                fprintf('Exiting...\n');
                break
            else
                fprintf('Enter another empty line to exit.\n');
                continue
            end
        else
            fprintf('Exiting...\n');
            break
        end
    else
        consecutive_empty_lines = 0;
    end
    tokens = regexp(deblank(in_str), '\s', 'split');
    % remove empty tokens because deblank doesn't get rid of all spaces
    tokens = tokens(~cellfun(@isempty, tokens));
    if length(tokens) ~= expected_tokens
        if length(tokens) == 1 && ...
                (strcmpi(tokens{1}, 'exit') || strcmpi(tokens{1}, 'quit'))
            fprintf('Exiting...\n');
            break
        end
        fprintf(2, sprintf('ERROR: expected %d whitespace-separated inputs\n', ...
            expected_tokens));
        fprintf(2, 'expected inputs (whitespace-delimited) are:\n');
        fprintf(2, ['(in_dir, image_name, image_suffix, starmask_suffix, '...
            'elps_fit_params_file,guide_dir, output_dir)\n']);
        continue
    end
    in_dir = tokens{1};
    image_name = tokens{2};
    image_suffix = tokens{3};
    starmask_suffix = tokens{4};
    elps_fit_params_filename = tokens{5};
    guide_dir = tokens{6};
    output_dir = tokens{7};
    fprintf(['running with in_dir=%s, image_name=%s, image_suffix=%s, '...
    'starmask_suffix=%s, output_dir=%s, guide_dir=%s\n'], in_dir, image_name, ...
    image_suffix, starmask_suffix, output_dir, guide_dir)
    image_run_id = sprintf('%s_%07d', server_run_id, nxtSeqNum-1);
    
    try
        batchFindClusterArcs(in_dir, {image_name}, image_suffix, starmask_suffix,...
            [], [], stgs, output_dir, elps_fit_params_filename, guide_dir, [], image_run_id);
    catch ME
        fprintf(2, '%s', sprintf('ERROR: %s\n', ME.message));
    end
end

end
