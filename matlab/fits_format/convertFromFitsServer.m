function convertFromFitsServer()

expected_tokens = 7;
nxtSeqNum = 1;
consecutive_empty_lines = 0;
max_consecutive_empty_lines = 1;
while true
    in_str = input('', 's');
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
        fprintf(2, ['(imgPath, starMaskPath, starMaskVal, qLvl, nReps,'...
                    'blackQlvl, outPath)\n']);
        continue
    end
    img_path = tokens{1};
    starmask_path = tokens{2};
    starmask_val = tokens{3};
    qLvl = tokens{4};
    nReps = tokens{5};
    blackQlvl = tokens{6};
    outPath = tokens{7};
%     image_run_id = sprintf('%s_%07d', server_run_id, nxtSeqNum-1);
%     if strcmpi(elps_fit_params_filename, 'NONE');
%         elps_fit_params = {};
%     else
%         elps_fit_file = fopen(elps_fit_params_filename, 'rt');
%         elps_fit_str = fgetl(elps_fit_file);
%         fclose(elps_fit_file);
%         elps_fit_params = {stringToDblStruct(elps_fit_str)};
%     end
    try
        convertFromFits(img_path, starmask_path, starmask_val, qLvl, ...
                        nReps, blackQlvl, outPath);
    catch ME
        fprintf(2, '%s', sprintf('ERROR: %s\n', ME.message));
    end
end

end
