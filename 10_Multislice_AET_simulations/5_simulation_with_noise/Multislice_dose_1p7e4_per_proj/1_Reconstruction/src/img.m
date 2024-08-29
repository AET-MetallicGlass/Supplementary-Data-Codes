function img(varargin)
%IMG displays the input object with given titles
%   ob: can be a 2D/3D numerica array or 1D cell with each cell as a 2D
%   numeric matrix
%   label: empty [] or a 1D cell same length as the z direction of *ob*
%   
%GUI instruction
%   'd' OR 'right arrow': next image
%   'a' OR 'left arrow': previous image
%   'w' OR 'up arrow': first image
%   's' OR 'down arrow': last image
%   'c': close figure
%
%
%
%Examples: main
%   img(ob1);
%       ob1: 2D/3D array
%   img(ob1, label1); 
%       ob1: 2D/3D array
%       label1:
%           empty [] (if array is 3D, then title of each 2D frame will be 
%           'frame: 1', 'frame: 2',...) OR
%           1D cell (each cell corresponding to each 2D frames) OR
%           a string (for all 2D frames)
%   img(ob1, label1, ob2, label2, ...);
%       ob1, ob2, ... have the same size in 3rd dimension (size(ob,3))
%
%Examples: modes
%   img(___, 'size', [1,3]);
%       generates subplot(1,3,kk)
%   img(___, 'axes', 'off');
%       remove white space between axes
%   img(___, 'position', 'nospace');
%       sets the figure position so that the mode ('axes','off') actually
%       removes all white space.
%       set to anything else but 'nospace' to deactivate this mode
%       only works if 'axes' mode is at 'off'
%   img(___, 'abs', 'on');
%       only plot absolute values of given pixel value
%       set to 'off' if wishes to plot actual real values (negative and
%       positive)
%   img(___, 'colormap', 'jet');
%       uses the colormap 'jet'
%       can be changed to any MATLAB colormap
%
%Remarks
%   If label is a string, then DO NOT USE mode names, e.g., 'size','axes',
%   etc.
%   
%   convert 3D arrays to cell arrays since cell arrays are faster when
%   slicing 2D arrays (it does not need to make a copy every time a 2D 
%   array is sliced)
%   should consider the case where Nz are not the same for all objects
%   
%Custom Functions
%   setAxes
%   cell2var
%   figsize
%
%   last updated 05/14/2018 by Arjun Rana
%   version 1.0
%
%   Author: Andrew Yuan
%   Jianwei (John) Miao Coherent Imaging Group
%   University of California, Los Angeles
%   Copyright (c) 2017, All Rights Reserved
%
%% determine modes
input = varargin;
mode = {'size','axes','position','abs','colormap','caxis'}; % case IN-sensitive
default = {[], 'off','nospace','on','jet','auto'}; % default values for each mode
mode_N = length(mode);
for jj = 1:mode_N
    if any(strcmpi(input, mode{jj}))
        item_idx = find(strcmpi(input, mode{jj})) +1;
        default{jj} = input{item_idx};

        input(item_idx-1:item_idx) = [];
    end
end

[layout, ax_mode, fig_pos, abs_mode, cmap, cax] = cell2var(default);
%% parameter preparation
inputN = length(input);
if inputN == 1
    ob_NN = 1;
    % generate object
    ob{1} = input{1};
    
    if iscell(ob{1}) % do nothing if 1D cell
        Nz = length(ob{1}); % assumed to be the same for all objects
    else % change 2D/3D array to 1D cell
        Nz = size(ob{1}, 3); % Nz >= 1
        buff = cell(1,Nz);
        for kk = 1:Nz
            buff{kk} = ob{1}(:,:,kk);
        end
        
        ob{1} = buff;
    end
    
    % generate labels
    label{1} = cell(1, Nz);
    % if Nz == 1, labels{1} = {[]};
    if Nz > 1
        for kk = 1:Nz
            label{1}{kk} = ['frame = ', int2str(kk)];
        end
    end
    
%     fullscr = 0;
else % inputN >= 2
    ob_NN = floor(inputN/2);
    ob = cell(1, ob_NN);
    label = cell(1, ob_NN);
    for ii = 1:ob_NN
        % generate each object
        ob{ii} = input{2*ii-1};
        if iscell(ob{ii}) % do nothing if 1D cell
            Nz = length(ob{ii}); % assumed to be the same for all objects
        else % change 2D/3D array to 1D cell
            Nz = size(ob{ii}, 3); % Nz > 1
            buff = cell(1,Nz);
            for kk = 1:Nz
                buff{kk} = ob{ii}(:,:,kk);
            end
            ob{ii} = buff;
        end
        
        
        label{ii} = cell(1, Nz);
        % generate labels for each object
        if ~isempty(input{2*ii})
            if ischar(varargin{2*ii}) 
                % if the input label is just a string (character vector),
                % then repeat the string and make a cell          
                if Nz > 1
                    for kk = 1:Nz
                        label{ii}{kk} = ['frame = ', int2str(kk)];
                    end
                    
                    fcnh = @(x) [input{2*ii},', ',x];
                    label{ii} = cellfun(fcnh,label{ii},'UniformOutput',false);
                else
                    label{ii} = {input{2*ii}};
                end
                
            else
                label{ii} = input{2*ii};
            end
        else
            % if Nz == 1, labels{ii} = {[]};
            if Nz > 1
                for kk = 1:Nz
                    label{ii}{kk} = ['frame = ', int2str(kk)];
                end
            end
        end
    end
    
%     % determine whether or not show full screen
%     fullscr = 0;
%     if rem(inputN, 2) == 1 && input{end} == 1
%         fullscr = 1;
%     end
end

if strcmp(abs_mode,'on')
    for ii = 1:ob_NN
        ob{ii} = cellfun(@abs, ob{ii}, 'UniformOutput', false);
    end
end
%% plot images
fr = 1;
if isempty(layout)
    if ob_NN == 1
        imagesc(ob{1}{fr}); axis image; colormap(cmap); caxis(cax);
        title(label{1}{fr});
    else % ob_NN >= 2
        for ii = 1:ob_NN
            subplot(1,ob_NN,ii);
            imagesc(ob{ii}{fr}); axis image; colormap(cmap); caxis(cax);
            title(label{ii}{fr});
        end
    end
else
    for ii = 1:min(ob_NN,prod(layout))
        subplot(layout(1),layout(2),ii);
        imagesc(ob{ii}{fr}); axis image; colormap(cmap); caxis(cax);
        title(label{ii}{fr});
    end
end

% full screen
% if fullscr == 1
%     set(gcf, 'position', fpos);
% %     % generates warning (suppressed)
% %     pause(0.00001);
% %     frame_h = get(handle(gcf),'JavaFrame');
% %     set(frame_h,'Maximized',1);
% end

if strcmpi(ax_mode,'off')
   setAxes([0.05,0,0,0]);
   
   if strcmpi(fig_pos,'nospace')
       fullpos = get(groot,'ScreenSize'); % full screen size
       fpos = get(gcf,'position');
       fsize = figsize(gcf);
       
       fpos(3) = fpos(4)/fsize(1) *0.95 *fsize(2); % set X size based on Y size
       if (fpos(1)+fpos(3)) >= fullpos(3)
           fpos(3) = fullpos(3) - fpos(1); % set X size all the way to the right end of screen
           fpos(4) = fpos(3)/(fsize(2)*0.95) *fsize(1);
           if (fpos(2)+fpos(4)) >= fullpos(4)
               fpos(4) = fullpos(4) - fpos(2);
           end
       end
       fpos=max(fpos,0);
       set(gcf, 'position', fpos);
   end
end

Nz = length(ob{1});
if Nz > 1
    % generate structure
    S.fr = fr;
    S.ax_mode = ax_mode;
    S.abs_mode = abs_mode;
    S.layout = layout;
    S.cmap = cmap;
    S.cax = cax;
    
    guidata(gcf,S);
    set(gcf,'WindowKeyPressFcn',{@fh_wkpfcn, ob, label});
end
end

%% --------------------------press key------------------------------
function fh_wkpfcn(fh, edata, ob, label)
S = guidata(fh);

fr = S.fr;
ax_mode = S.ax_mode;
layout = S.layout;
cmap = S.cmap;
cax = S.cax;
kp = edata.Key;

Nz = length(ob{1});
ob_NN = length(ob);
switch kp
    case {'d','rightarrow'}
        fr = fr+1;
    case {'a','leftarrow'}
        fr = fr-1;
    case {'w','uparrow'}
        fr = 1;
    case {'s','downarrow'}
        fr = Nz;
    case 'f'
        fr = input('Show frame number: ');
end
fr = mod(fr-1,Nz)+1;

if isempty(layout)
    if ob_NN == 1
        imagesc(ob{1}{fr}); axis image; colormap(cmap);caxis(cax);
        title(label{1}{fr});
    else % ob_NN >= 2
        for ii = 1:ob_NN
            subplot(1,ob_NN,ii);
            imagesc(ob{ii}{fr}); axis image; colormap(cmap);caxis(cax);
            title(label{ii}{fr});
        end
    end
else
    for ii = 1:min(ob_NN,prod(layout))
        subplot(layout(1),layout(2),ii);
        imagesc(ob{ii}{fr}); axis image; colormap(cmap);caxis(cax);
        title(label{ii}{fr});
    end
end

if strcmpi(ax_mode,'off')
   setAxes;
end

S.fr = fr;

guidata(gcf, S);

if strcmp(kp,'c')
    close(gcf);
end

end

% ----------------------- custom functions --------------------------------
% =========================================================================
% setAxes
% =========================================================================
function setAxes(varargin)
%SETAXES automatically sets the position and size of the axes (subplots) in
%the given figure so that white spaces are removed
%
%Optional Inputs
%   pad: A vector of size 1,2,3,4. The padding method follows that of HTML
%       If size 1, then each subplot is padded with white space by the 
%       given amount
%       If size 2, then for each subplot, the top and bottom are padded by 
%       pad[0], and the right and left are padded by pad[1]
%       If size 3, then for each subplot, the top is padded by pad[0], the
%       right and left are padded by pad[1], and the bottom is padded by 
%       pad[2]
%       If size 4, then for each subplot, the top is padded by pad[0], the 
%       right by pad[1], the bottom by pad[2], the left by pad[3]
%   fh: figure handle. By default, == gcf
%
%Examples
%   setAxes
%       set axes positions for current figure (gcf)
%   setAxes(0.1)
%       each subplot in current figure has boundary white spaces of 0.1 (in
%       normalized units)
%   setAxes([0.05,0,0,0])
%       each subplot in current figure has top boundary (white space) of 
%       0.05 and 0 in the right, bottom and left boundary
%   setAxes([0.05,0,0,0],fh)
%       set axes for subplots in the figure with figure handle fh
%
%Custom Functions
%   cell2var
%   figsize
%
%   last updated 09/22/2017
%   version 1.0
%
%   Author: Andrew Yuan
%   Jianwei (John) Miao Coherent Imaging Group
%   University of California, Los Angeles
%   Copyright (c) 2017, All Rights Reserved
%
%% optional inputs
input = {[0.05,0,0,0], gcf};
input_N = length(input);
for jj = 1:input_N
    if nargin >= jj && ~isempty(varargin{jj})
        input{jj} = varargin{jj};
    end
end
[pad, fh] = cell2var(input);

%% parameter preparation
switch length(pad)
    case 1
        pad = [pad, pad, pad, pad]; % top right bottom left
    case 2
        pad = [pad(1),pad(2),pad(1),pad(2)]; % top right bottom left
    case 3
        pad = [pad(1),pad(2),pad(3),pad(2)]; % top right bottom left
    case 4
        % do nothing
    otherwise
        error('padding size is invalid');
end

%% main
num = figsize(fh);
ax = findobj(fh,'type','axes');
axN = length(ax);

axX = 1/num(2)-pad(2)-pad(4); % X direction
axY = 1/num(1)-pad(1)-pad(3); % Y direction


for ii = 1:num(1)
    posY = pad(3)+1-ii/num(1);
    for jj = 1:num(2)
        posX = pad(4)+(jj-1)/num(2);
        
        kk = (ii-1)*num(2)+jj;
        set(ax(axN+1-kk),'units','normalized',...
                'position',[posX,posY,axX,axY],...
                'XTick','','YTick','');
    end
end

end

% =========================================================================
% figsize
% =========================================================================
function num = figsize(varargin)
%FIGSIZE returns the number of figures in the given figure handle (usually
%gcf)
%   num: 1-by-2 vector specifying the number of rows and columns of
%   subplots in the given figure
%
%Optional Inputs
%   fh: figure handle. By default, == gcf
%
%   last updated 09/22/2017
%   version 1.0
%
%   Author: Andrew Yuan
%   Jianwei (John) Miao Coherent Imaging Group
%   University of California, Los Angeles
%   Copyright (c) 2017, All Rights Reserved
%
%% optional input
fh = gcf;
if nargin >= 1 && ~isempty(varargin{1})
    fh = varargin{1};
end

%%
ax = findobj(fh,'type','axes');
pos = get(ax,'position');
if iscell(pos)
    pos = cell2mat(pos);
    colNum = numel(unique(pos(:,1))); % same X positions
    rowNum = numel(unique(pos(:,2))); % same Y positions
else
    colNum = 1;
    rowNum = 1;
end

num = [rowNum, colNum];

end

% =========================================================================
% cell2var
% =========================================================================
function [varargout] = cell2var(cellA, varargin)
%CELL2VAR stands for cell to variable. It converts a cell of length N to N
%specified variables
%   cellA: 1D cell
%
%Optional Inputs
%   num: a # specifying which cell to convert to variable
%
%   last updated 06/05/2017
%   version 1.0
%
%   Author: Andrew Yuan
%   Jianwei (John) Miao Coherent Imaging Group
%   University of California, Los Angeles
%   Copyright (c) 2017, All Rights Reserved
%

%% possible case
if nargin >= 2 && ~isempty(varargin{1})
    varargout{1} = cellA{varargin{1}};
    return;
end
%% main
if nargout <= length(cellA)
   for kk = 1:nargout
       varargout{kk} = cellA{kk};
   end
else
    error('number of outputs exceeds length of input cell');
end

end


