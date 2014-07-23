function mesh = MakeTriKRegComplex(k, numRings, geometry)
% MakeTriKRegComplex -- Make a planar k-regular complex from triangles
%  Usage
%    mesh = MakeTriKRegComplex(k, numRings [, geometry])
%  Input
%    k          number of spokes
%    numRings   number of rings
%    geometry   a cell, indexed by vertex number, containing jet data
%
%  Output
%    mesh       Planar (well depends on input "geometry" data) k-regular complex made of triangles
%                mesh{1}: vertex list
%                mesh{2}: face list
%                mesh{3}: neighbor list
% 

if k<3, error('Cannot construct anything smaller than k=3.'); end;
if numRings<1, error('Must have one or more rings.'); end;

VL=[0 0]; % initialize center
FL=[];

for r=1:numRings
    for s=0:k-1
        t1=(s/k)*2*pi;
        curpos=[cos(t1) sin(t1)]*r;
        t2=((s+1)/k)*2*pi;
        nextpos=[cos(t2) sin(t2)]*r;
        diff=nextpos-curpos;
        
        for o=0:r-1
            VL=[VL; curpos+diff*(o/r)];
            cv=getTriKRegVertex(k,r,s,o);
            if(o==0) % spoke
                FL=[FL;...
                    cv,...
                    getTriKRegVertex(k, r, s, 1),... % next (counter clockwise)
                    getTriKRegVertex(k, r-1, s, 0)]; % inside
            else % bead
                FL=[FL;...
                    cv,...
                    getTriKRegVertex(k, r-1, s, o),...
                    getTriKRegVertex(k, r-1, s, o-1)];
                FL=[FL;...
                    cv,...
                    getTriKRegVertex(k, r, s, o+1),...
                    getTriKRegVertex(k, r-1, s, o)];
            end
        end
    end
end

% replace VL if geometry is specified
if nargin>2
    rep=find(cellfun('isempty',geometry)==0); % which vertices are specified
    rows=size(geometry{rep(1)},1);
    cols=size(geometry{rep(1)},2);
    VL=zeros(size(VL,1),rows,cols);
    for cur=rep; VL(cur,:,:)=geometry{cur}; end;
end

mesh{1}=VL;
mesh{2}=FL;
mesh{3}=MakeTriKRegComplexNL(k, numRings);

% 
% Copyright (c) 2005. Kyle McDonald
%