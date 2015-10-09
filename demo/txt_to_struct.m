function [A] = txt_to_struct(filename)
    A = struct;
    fid = fopen(filename);
    line = fgetl(fid);
    while ischar(line)
        splitline = strsplit(line);
        numVals = size(splitline,2);
        vals = [];
        for i=2:numVals
            vals = [vals str2num(splitline{i})];
        end
        A.(splitline{1}) = vals;
        line = fgetl(fid);
    end
end
