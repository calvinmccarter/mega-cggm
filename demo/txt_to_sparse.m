function [A] = txt_to_sparse(filename)
    data = dlmread(filename);
    A = sparse(data(2:end,1), data(2:end,2), data(2:end,3), ...
        data(1,1), data(1,2));
end
