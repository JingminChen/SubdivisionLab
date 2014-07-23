function [cm,submat]=FormTriCharacteristicMap(NameOfScheme,k,SizeInvarNhbd,NumIter,jetorder)
% FormTriSubdivisionMatrix -- Compute an approximation of a characteristic map
%  Usage
%    cm=FormTriCharacteristicMap(NameOfScheme,k,SizeInvarNhbd,NumIter,jetorder)
%  Input
%    NameOfScheme    name of subdivision shceme. e.g.: 'Loop', 'MButterfly'...
%    k               valence of extraordinary vertex
%    SizeInvarNhbd   size of the invariant neighborhood of the scheme
%    NumIter         number of subdivision steps used to approximate the
%                    characteristic map (default=4)
%    jetorder        order of jet data used by scheme (default=0)
%  Output
%    cm              approximation to characteristic map, has the same
%                    connectivity of a k-regular complex with (2^NumIter) rings
% 

if nargin<4, NumIter=4; end
if nargin<5, jetorder=0; end

m = nchoosek(jetorder+2,2); % multiplicity of an order r jss

max=getTriKRegVertex(k,SizeInvarNhbd+1,0,0)-1;
% Form Subdivision Matrix
submat=FormTriSubdivisionMatrix(NameOfScheme,k,SizeInvarNhbd,jetorder);
% Compute the 2 subdominant eigenvectors
OPTS.disp=0; [u,v]=eigs(submat,3,'LM',OPTS);
u1=u(:,2); u2=u(:,3);
% Symmetrize the 2 subdominant eigenvectors
T = inv([[u1(m+1) u2(m+1)]; ...
        [u1(m+m+1) u2(m+m+1)]]) * [1 0; cos(2*pi/k) sin(2*pi/k)];
X = T(1,1)*u1 + T(2,1)*u2;
Y = T(1,2)*u1 + T(2,2)*u2;

% Form the "eigen- control net" used to define the characteristic map
for i=1:max
    seg = ((i-1)*m+1):(i*m);
    subeig{i}=[X(seg)';Y(seg)'];
end;
cm=MakeTriKRegComplex(k,SizeInvarNhbd,subeig);

NumVertices = size(cm{1},1);
data = zeros(NumVertices,1); data(1)=1;
cm{1} = [cm{1} data];

% subdivide it a few times
for j=1:NumIter, 
    cm=TriSubdivision(cm,NameOfScheme);
    cm=ReorderTriKRegComplex(k,2^(j-1)*SizeInvarNhbd,cm);
end
% clcm=TriSubdivisionLimit(cm,NameOfScheme);
% remove the unwanted rings in cm:
cm = TrimTriKRegComplex(cm,k,2^NumIter);

bf = cm; cm{1}(:,3)=[]; % remove the data
RenderMesh(cm,0,0); % Plot the characteristic map on the x-y plane
% hold on
% RenderMesh(bf,1); % Plot the basis function in characteristic coordinates

% daspect([1,1,.5])
% pbaspect([1 1 .5])

% 
% Copyright (c) 2005. Kyle McDonald and Thomas P.-Y. Yu
%
% Copyright (c) 2006. Thomas P.-Y. Yu
%