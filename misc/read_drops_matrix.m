function M = read_drops_matrix (filename)
% Read a sparse matrix from a file written by DROPS operator<<.
% This creates a sparse matrix.

fid= fopen( filename);
if fid == -1
     error( 'Could not open file.');
end

C= textscan( fid, '%n%n%f', 'commentStyle', '%');
fclose( fid);
M= sparse( C{1}, C{2}, C{3});
end
