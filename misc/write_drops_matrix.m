function [] = write_drops_matrix (S, f)
% Write a sparse matrix M to the file f. This can be read by DROPS operator>>.

fid= fopen( f, 'w');
if fid == -1
     error( 'Could not open file.');
end

[r, c]= size( S);
nz= nnz( S);
fprintf( fid, '%% %.0f x %.0f %.0f nonzeros\n', r, c, nz);

[row, col, v]= find( S);
fprintf( fid, '%.0f %.0f %.15g\n', [row, col, v]');

fclose( fid);
end
