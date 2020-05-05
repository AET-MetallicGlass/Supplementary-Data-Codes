%%  FileExist %%
% check files whether exist

function Ex = FileExist(FileName)
dirFile = dir(FileName);
if length(dirFile) == 1
    Ex = ~(dirFile.isdir);
else
    Ex = false;
end
end