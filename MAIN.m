clear all
close all
% This is for now the main analysis file to be made


%% Parameters
L = 5;    % mm
h = 1;     % mm

N = L/h;      % Number of Horizontal Nodes

E = 210E9;
mu = 0.3;
rho = 7800;
alpha = 10;
beta = 2;

[K,M,R,nodes,links,beam] = buildBeam(L,h,N,E,mu,rho,alpha,beta);

              
% q -> State vector!! (displacement vector)
%% Settings
plotSettings.nodes = true;
plotSettings.nodeNumbers = true;
plotSettings.links = true;
plotSettings.linkNumbers = false;


numNodes = length(nodes);
numLinks = length(links);
q = zeros(numNodes,1);          % State vector!

figure()
hold on
axis equal
ax = gca;
P = [];

% Plot Links
for i = 1:numLinks
    % Links
    if plotSettings.links == true
        linkPlots = links{i}.plotLink(ax,[0 0.4470 0.7410]);
        P = [P,linkPlots];
    end
    
    % Link Numbers
    if plotSettings.linkNumbers == true
        linkNumberPlots = links{i}.plotNumber(ax,'k');
        P = [P,linkNumberPlots];
    end
end

% Plot Nodes
for i = 1:numNodes
    % Nodes
    if plotSettings.nodes == true
        nodePlots = nodes{i}.plotNode(ax,'k');
        P = [P,nodePlots];
    end
    
    % Node Numbers
    if plotSettings.nodeNumbers == true
        nodeNumberPlots = nodes{i}.plotNumber(ax,'k');
        P = [P,nodeNumberPlots];
    end
end






              
              
              
              
              
              
              
              