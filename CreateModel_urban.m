%_________________________________________________________________________%
%  Urban Canyon Terrain Model for UAV Path Planning                       %
%  Dense tall obstacles simulating city buildings                         %
%_________________________________________________________________________%

function model=CreateModel_urban()

    H = imread('ChrismasTerrain.tif');
    H (H < 0) = 0;
    MAPSIZE_X = size(H,2);
    MAPSIZE_Y = size(H,1);
    [X,Y] = meshgrid(1:MAPSIZE_X,1:MAPSIZE_Y);

    % Dense building-like threats (urban canyon)
    % [x, y, z_center, radius]
    threats = [
        300, 300, 120, 40;   % Building 1
        400, 250, 130, 35;   % Building 2
        500, 350, 140, 45;   % Building 3
        350, 500, 125, 40;   % Building 4
        600, 400, 135, 50;   % Building 5
        450, 600, 130, 35;   % Building 6
        650, 550, 145, 40;   % Building 7
        550, 700, 140, 45;   % Building 8
        700, 650, 150, 35;   % Building 9
        300, 650, 120, 30;   % Building 10
    ];

    % Map limits
    xmin=1; xmax=MAPSIZE_X;
    ymin=1; ymax=MAPSIZE_Y;
    zmin=100; zmax=200;

    % Start and end (diagonal crossing through buildings)
    start_location = [100; 100; 150];
    end_location = [750; 750; 150];

    n=10;

    model.start=start_location;
    model.end=end_location;
    model.n=n;
    model.xmin=xmin; model.xmax=xmax;
    model.zmin=zmin; model.ymin=ymin;
    model.ymax=ymax; model.zmax=zmax;
    model.MAPSIZE_X=MAPSIZE_X; model.MAPSIZE_Y=MAPSIZE_Y;
    model.X=X; model.Y=Y; model.H=H;
    model.threats = threats;
    PlotModel(model);
    title('Urban Canyon Terrain');
end
