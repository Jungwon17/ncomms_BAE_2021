function [T, W, H, P] =  NTT2Mat(fNm)

% fNm = 'D:\Cheetah_data\classical_conditioning\CC-SOM-ChR4\2016-03-31_AD_1.30DV\TT2.ntt';
fID = fopen(fNm);

H = {};
for iL = 1:49
    tline = fgetl(fID);
    H = [H; tline];
end
fseek(fID, 0, 'eof');
endPt = ftell(fID);
nP = floor(endPt - 16*1024)/(8+4+4+4*8+2*32*4);

fseek(fID, 16*1024, 'bof');

T = uint64([]);
S = uint32([]);
C = uint32([]);
P = uint32([]);
W = int16(zeros(0,4,32));

% for iP = 1:nP
for iP = 1:100
    T = [T; fread(fID, 1, 'uint64')];
    S = [S; fread(fID, 1, 'uint32')];
    C = [C; fread(fID, 1, 'uint32')];
    P = [P; fread(fID, [1, 8], 'uint32')];
    W = [W; reshape(fread(fID, [4 32], 'int16'), 1, 4, 32)];
end

fclose(fID);