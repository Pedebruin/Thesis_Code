clear all
close all
% This is for now the main analysis file to be made


%% Parameters
L = 5;    % mm
h = 1;     % mm

N = 10;      % Number of Horizontal Nodes

E = 210E9;
mu = 0.3;
rho = 7800;

%% Setup
% Create container
beam = createpde('structural','modal-planestress');

% Define Geometry
rects = zeros(10,N-1);
ns = [];
sf = [];
for i = 1:N-1
    rects(:,i) = [3,4,(i-1)*L/(N-1),(i-1)*L/(N-1),i*L/(N-1),i*L/(N-1),0, h, h, 0]';
    n = ['Seg',num2str(i)];
    ns = [ns;n];
    
    if ~isempty(sf)
        sf = [sf,'+',n];
    else
        sf = n;
    end
end

g = decsg(rects,sf,ns');
geometryFromEdges(beam,g);

pdegplot(beam,'EdgeLabels','on');

I = h^3/12;
analyticalOmega1 = 3.516*sqrt(E*I/(L^4*(rho*h)))/(2*pi);

% Apply material to structure
structuralProperties(beam,'YoungsModulus', E, ...
                           'PoissonsRatio', mu,...
                           'massDensity', rho);
                       
% Set up boundary conditions
structuralBC(beam,'edge',1,'Constraint','fixed');

% Generate Mesh
generateMesh(beam,'Hmax',h,...
                    'Hmin',h);
pdeplot(beam);
         
              
              
              
              
              
              
              
              
              