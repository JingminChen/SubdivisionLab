function mat = FormTriSubdivisionMatrix(scheme, k, neighborhood, jetorder)
% FormTriSubdivisionMatrix -- Form the subdivision matrix for a given subdivision
%                             scheme around a vertex of specified valence
%  Usage
%    mat = FormTriSubdivisionMatrix(scheme, k, neighborhood, jetorder)
%  Input
%    scheme          name of subdivision shceme. e.g.: 'Loop', 'MButterfly'...
%    k               valence of central vertex
%    neighborhood    invariant neighborhood of the scheme
%    jetorder        order of jet data used by scheme (default=0)
%  Output
%    n               integer, vertex number
% 

if nargin<4, jetorder=0; end
m = nchoosek(jetorder+2,2); % multiplicity of an order r jss
mat=[];
max=getTriKRegVertex(k,neighborhood+1,0,0)-1;

for i=1:max, geom{i}=zeros(1,m); end
coarse=MakeTriKRegComplex(k,neighborhood,geom);
for cv=1:max,
    for nu=1:m,
        coarse{1}(cv,1,nu)=1;     
        fine=ReorderTriKRegComplex(k,neighborhood,TriSubdivision(coarse,scheme));        
        column = [];
        for i=1:max,
            segment = shiftdim(fine{1}(i,:,:),1);
            column = [column; segment'];
        end
        mat = [mat column];
        coarse{1}(cv,1,nu)=0;
    end
end

% 
% Copyright (c) 2005. Kyle McDonald
%