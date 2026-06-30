%{
Map 4: Christmas terrain with 7 threat cylinders
- Dense uniform threats (all R=70m) covering the entire map
%}

function model = CreateModel_map4()

    H = imread('ChrismasTerrain.tif'); % Get elevation data
    H (H < 0) = 0;
    MAPSIZE_X = size(H,2); % x index: columns of H
    MAPSIZE_Y = size(H,1); % y index: rows of H
    [X,Y] = meshgrid(1:MAPSIZE_X,1:MAPSIZE_Y); % Create all (x,y) points to plot

    % Threats as cylinders (7 threats, uniform R=70m)
    R1=70;  % Radius
    x1 = 400; y1 = 700; z1 = 150; % center (upper-left)

    R2=70;  % Radius
    x2 = 350; y2 = 530; z2 = 150; % center (mid-left)

    R3=70;  % Radius
    x3 = 390; y3 = 170; z3 = 150; % center (lower-left)

    R4=70;  % Radius
    x4 = 520; y4 = 360; z4 = 150; % center (center)

    R5=70;  % Radius
    x5 = 600; y5 = 550; z5 = 150; % center (upper-center)

    R6=70;  % Radius
    x6 = 620; y6 = 210; z6 = 150; % center (lower-right)

    R7=70;  % Radius
    x7 = 680; y7 = 700; z7 = 150; % center (upper-right)

    % Map limits
    xmin= 1;
    xmax= MAPSIZE_X;

    ymin= 1;
    ymax= MAPSIZE_Y;

    zmin = 100;
    zmax = 200;

    % Start and end position
    start_location = [200;100;150];
    end_location = [800;800;150];

    % Number of path nodes (not including the start position (start node))
    n=10;

    % Incorporate map and searching parameters to a model
    model.start=start_location;
    model.end=end_location;
    model.n=n;
    model.xmin=xmin;
    model.xmax=xmax;
    model.zmin=zmin;
    model.ymin=ymin;
    model.ymax=ymax;
    model.zmax=zmax;
    model.MAPSIZE_X = MAPSIZE_X;
    model.MAPSIZE_Y = MAPSIZE_Y;
    model.X = X;
    model.Y = Y;
    model.H = H;
    model.threats = [x1 y1 z1 R1; x2 y2 z2 R2; x3 y3 z3 R3; x4 y4 z4 R4; x5 y5 z5 R5; x6 y6 z6 R6; x7 y7 z7 R7];
    PlotModel(model);
end
