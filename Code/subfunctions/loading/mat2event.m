function tFile = mat2tt(mFile)
%mat2tt modifies matrix file to t file
narginchk(1,1);

predir = 'C:\\Users\\Lapis\\OneDrive\\project\\workingmemory_interneuron\\data\\';
curdir = 'D:\\Cheetah_data\\workingmemory_interneuron\\';
mFile = regexprep(mFile,predir,curdir);

preext = '.mat';
curext = '.t';

tFile = regexprep(mFile,preext,curext);