function sparse_to_txt(filename, A)
    [r, c, v] = find(A);
    data = [size(A,1) size(A,2) nnz(A); r c v];
    dlmwrite(filename, data, 'delimiter', ' ', 'precision', 15);
end
