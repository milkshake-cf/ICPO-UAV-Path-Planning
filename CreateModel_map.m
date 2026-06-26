%_________________________________________________________________________%
%  8 benchmark scenarios (inspired by SPSO paper Sec.5)                   %
%  Each modifies: obstacles, start/end, altitude limits                   %
%_________________________________________________________________________%

function model=CreateModel_map(mapID)

    H = imread('ChrismasTerrain.tif');
    H(H < 0) = 0;
    MAPSIZE_X = size(H,2);
    MAPSIZE_Y = size(H,1);
    [X,Y] = meshgrid(1:MAPSIZE_X,1:MAPSIZE_Y);
    
    switch mapID
        case 1  % Original (sparse obstacles, long range)
            threats = [400 500 100 80; 600 200 150 70; 500 350 150 80; 
                       350 200 150 70; 700 550 150 70; 650 750 150 80];
            start_loc = [200;100;150];
            end_loc   = [800;800;150];
            zmin=100; zmax=200;
            
        case 2  % Dense urban (many small obstacles)
            threats = [300 300 120 50; 400 400 130 45; 500 300 125 55;
                       350 500 120 50; 450 600 135 45; 600 400 130 50;
                       300 700 125 55; 550 150 120 45; 200 500 130 50;
                       700 650 125 55];
            start_loc = [100;50;150];
            end_loc   = [750;750;150];
            zmin=80; zmax=220;

        case 3  % Sparse threats, high altitude
            threats = [500 500 120 100; 300 300 110 90];
            start_loc = [100;700;180];
            end_loc   = [700;100;180];
            zmin=120; zmax=250;
            
        case 4  % Narrow corridor
            threats = [200 400 100 60; 600 400 100 60; 400 200 100 60;
                       400 600 100 60];
            start_loc = [100;100;150];
            end_loc   = [700;700;150];
            zmin=100; zmax=200;
            
        case 5  % Edge obstacles, center clear
            threats = [100 400 100 50; 100 200 100 50; 100 600 100 50;
                       700 400 100 50; 700 200 100 50; 700 600 100 50];
            start_loc = [400;100;150];
            end_loc   = [400;700;150];
            zmin=100; zmax=200;
            
        case 6  % Mountain pass (wide corridor)
            threats = [300 450 100 50; 550 350 90 45];
            start_loc = [100;400;160];
            end_loc   = [750;400;160];
            zmin=80; zmax=250;
            
        case 7  % Clustered threats but passable
            threats = [350 350 120 70; 500 500 120 70; 200 600 110 60;
                       600 200 110 60];
            start_loc = [100;100;150];
            end_loc   = [750;750;150];
            zmin=100; zmax=200;
            
        case 8  % Low altitude weaving
            threats = [200 400 70 35; 400 400 70 35; 600 400 70 35];
            start_loc = [50;200;100];
            end_loc   = [750;600;100];
            zmin=60; zmax=180;
            
        otherwise
            error('Map ID must be 1-8');
    end
    
    model.start=start_loc;
    model.end=end_loc;
    model.n=10;
    model.xmin=1; model.xmax=MAPSIZE_X;
    model.ymin=1; model.ymax=MAPSIZE_Y;
    model.zmin=zmin; model.zmax=zmax;
    model.MAPSIZE_X=MAPSIZE_X; model.MAPSIZE_Y=MAPSIZE_Y;
    model.X=X; model.Y=Y; model.H=H;
    model.threats=threats;
    PlotModel(model);
end
