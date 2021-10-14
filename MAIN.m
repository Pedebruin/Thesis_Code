clear all
close all
set(0,'defaultTextInterpreter','latex','defaultAxesFontSize',12);  

%% Parameters & Settings
% Parameters
L = 100;    % mm
h = 1;     % mm

N = L/h;      % Number of Horizontal Nodes

E = 210E9;
mu = 0.3;
rho = 7800;
alpha = 10;
beta = 2;

% Settings
plotSettings.realMesh = true;
plotSettings.nodes = true;
plotSettings.nodeNumbers = false;
plotSettings.links = true;
plotSettings.linkNumbers = false;


% build FEM model of the beam
[K,M,R,nodes,links,beam] = buildBeam(L,h,N,E,mu,rho,alpha,beta,plotSettings);
numNodes = length(nodes);
numLinks = length(links);


% Find current state:
q = zeros(size(K,1),1);

[V,D] = eigs(K,M);
D = sort(real(sqrt(D)))/(2*pi);

% Update Beam
updateBeam(links,nodes,q);

% Plot beam
plotBeam(nodes,links,plotSettings);





              
              
              
              
function [links,nodes] =  updateBeam(links,nodes,q)
    numNodes = length(nodes);
    % Update beam
    % Nodes
    for i = 1:numNodes
        nodes{i}.update(q);
    end
    % Links
    for i = 1:length(links)
        links{i}.update(nodes);
    end
end             
              
              
              