function h = irf_customize_subplot(varargin)
%This function draw the skeleton of subplots
% Input:
%  Row = Number of rows 
%  Column = Number of columns
%  x0 = left margin  
%  y0 = bottom margin 
%  xgap = gap between columns
%  ygap = gap between rows
%  ytop = top margin 
%  xright = right margin 
%
% Examples 
% h = irf_customize_subplot0(Row, Column, 'xg',0.07, 'yg', 0.07)
%
% h = irf_customize_subplot(2, 2, 'x0', 0.07, 'y0', 0.7)
% h(1) = irf_plot(x)
% h(2) = irf_plot(y)
% h(3) = irf_plot(z)
% h(4) = irf_plot(q)
%
% The handles of the axis are return in 'h '
% following way axis handle is assign in 2 X 2 subplots arrangement  
%  (1,1) -> h(1), (2,1) -> h(2), (1,2) -> h(3), (2,2) -> h(4)
%
% Author : Ajay Lotekar,  Jun 22, 2021
%%
args = varargin;
nargs = nargin;
orignal_args=args;

if nargs == 0 % show only help
    help irf_subplot;
    return
end

if nargs<2 || ~(mod(nargs,2)==0)
    error('IRFU_MATLAB:irf_customize_subplot:InvalidNumberOfInputs','Incorrect number of input arguments')
end


%% Stack plot genration mesurment
No_of_Row = args{1};            % No. of Row
No_of_column = args{2};        %No. of Column
nR = 1:No_of_Row;        % Row no. matrix
nC = 1:No_of_column;   % Column no. matrix

% ---- default values --------------
x0 = 0.08;    % starting x co-ordinate
y0 = 0.16;      % Starting y co-ordinate
xg = 0.07;    % x gap (g)
yg = 0.00;    % y gap (g)
xR = 0.04;    % x right side space
yR = 0.07;      % y top side space
% -----------------------------------

Lth = []; Hgt=[];
while ~isempty(args)
    x = args{1}; args(1) = [];
    switch lower(x)
        case {'x0'}
            x0 = args{1}; args(1) = [];
        case {'y0'}
            y0 = args{1}; args(1) = [];
        case {'xgap'}
            xg= args{1}; args(1) = [];
        case {'ygap'}
            yg= args{1}; args(1) = [];
        case {'xright'}
            xR = args{1}; args(1) = [];
        case {'ytop'}
            yR = args{1}; args(1) = [];
        case {'lth'}
            Lth = args{1}; args(1) = [];
        case {'hgt'}
            Hgt = args{1}; args(1) = [];
    end
end

xL = 1-x0 - (max(nC) -1)*xg - xR; % Available leght for plots
yL = 1-y0 - (max(nR)-1)*yg - yR;   % Available hight for plots
if isempty(Lth)
    Lth = xL/max(nC);  % Length
end
if isempty(Hgt)
    Hgt = yL/max(nR);  % Hight
end
k=1;


for ju=1:No_of_column
    for iu=1:No_of_Row
        h(k)= subplot('position' ,[(x0+(nC(ju)-1)*xg+(nC(ju)-1)*Lth) (y0+(No_of_Row-iu)*yg+(No_of_Row-iu)*Hgt) Lth Hgt])
        k=k+1; 
    end
end

end