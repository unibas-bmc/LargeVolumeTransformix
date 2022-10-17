function [A] = tparamreplaceline(file1,file2,varargin)
% tparamreplaceline(tparam_old,tparam_new,lineno1,strnew1,lineno2,strnew2,lineno3,strnew3);
        
if nargin < 3
    error('Need to give line number and string...\n')
elseif not(rem(nargin, 2) == 0)
    error('Need to give inputfile, outputfile, and pairs of line number and string...\n')
else
    noLines = (nargin-2)/2;
    linenolist = zeros(1,noLines);
    strlist = {'';''};
    for i = 1:noLines
        linenolist(i) = varargin{2*i-1};
        strlist{i} = varargin{2*i};
    end
end

fid = fopen(file1,'r');
i = 1;
tline = fgetl(fid);
A{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
end
fclose(fid);

% Change cell A
for l = 1:noLines
    A{linenolist(l)} = strlist{l};
end

% Write cell A into txt
fid = fopen(file2, 'w');
for i = 1:numel(A)
    if A{i+1} == -1
        fprintf(fid,'%s', A{i});
        break
    else
        fprintf(fid,'%s\n', A{i});
    end
end

end

