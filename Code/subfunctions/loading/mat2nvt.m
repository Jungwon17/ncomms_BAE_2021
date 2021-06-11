function vFile = mat2nvt(mFile)
%mat2nvt modifies matrix file to nvt file
narginchk(1,1);

predir = 'C:\\Users\\Lapis\\OneDrive\\project\\workingmemory_interneuron\\data\\';
curdir = 'D:\\Cheetah_data\\workingmemory_interneuron\\';
vFile = [fileparts(regexprep(mFile,predir,curdir)), '\VT1.nvt'];