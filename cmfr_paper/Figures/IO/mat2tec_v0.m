% Simple export of text and MATLAB arrays to native tecplot data format

%% Example data
x = (0:0.01:2*pi)';
y = cos(x) + sin(x);

%% Set write parameters
strVar = 'x, y';
varMat = [x, y];
strFile = 'test1.dat';

%% Write file
strTitle = ['Variables = ', strVar];
dlmwrite(strFile, strTitle, 'delimiter', '');
dlmwrite(strFile, varMat, 'delimiter', '\t', '-append');
