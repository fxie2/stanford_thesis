function lattable=LatexTable(arr,strformat,captionString,colformat,rowformat)
% lattable=LatexTable(arr,colformat)
%
% Produces latex table from matlab cell array
% 
  sarr=size(arr);
  % loop over each colum to make header:
  header = '\\begin{table}[htbp]\n\\centering\n';
  header=[header '\\begin{tabular}{'];
  for j=1:sarr(2)
    if nargin==3 % Is column format defined?
      header=[header 'c|'];
    else
      header=[header colformat{j} ];
    end
  end
  lattable=[header '}'];
  lattable=[lattable ' \n ' ];
  % loop over each row:
  for j=1:sarr(1)
    row=parseentry(arr{j,1},strformat);
    % loop over each column
    for k=2:sarr(2);
      row=[row ' & ' parseentry(arr{j,k},strformat)];
    end
    row=[row ' \\\\'];
    lattable=[lattable ' \n ' row];
    lattable=[lattable ' \n ' rowformat{j}];
  end
  lattable=[lattable ' \n ' '\\end{tabular}'];
  lattable=[lattable '\n' captionString];
  lattable=[lattable ' \n ' '\\end{table}\n'];
end
  
function ret=parseentry(entry,strformat)
  if isempty(entry)
    entry=' ';
  end
  if isstr(entry)
    ret=entry;
  else
    ret=sprintf(strformat,entry);
    end
  end