%{
Map 3: Christmas terrain with 4 threat cylinders
- Map 2 threats + 1 additional threat (lower-left area)
%}

function model = CreateModel_map3()

    H = imread('ChrismasTerrain.tif'); % Get elevation data
    H (H < 0) = 0;
    MAPSIZE_X = size(H,2); % x index: columns of H
    MAPSIZE_Y = size(H,1); % y index: rows of H
    [X,Y] = meshgrid(1:MAPSIZE_X,1:MAPSIZE_Y); % Create all (x,y) points to plot

    % Threats as cylinders (4 threats = Map 2 + 1 extra)
    R1=80;  % Radius
    x1 = 420; y1 = 490; z1 = 150; % center

    R2=75;  % Radius
    x2 = 580; y2 = 360; z2 = 150; % center

    R3=80;  % Radius
    x3 = 620; y3 = 220; z3 = 150; % center

    R4=60;  % Radius
    x4 = 280; y4 = 310; z4 = 150; % center (new, lower-left)

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
    model.threats = [x1 y1 z1 R1; x2 y2 z2 R2; x3 y3 z3 R3; x4 y4 z4 R4];
    PlotModel(model);
end
