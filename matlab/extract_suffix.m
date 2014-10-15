function [suffix belong] = extract_suffix(file_name, suffix_model, pos_end)
min_len = pos_end + length(suffix_model);
if (length(file_name) >= min_len)
    suffix = file_name(end - pos_end - length(suffix_model) + 1 : end - pos_end);
    belong = strcmp(suffix, suffix_model);
else
    suffix = ' ';
    belong = 0;
end
