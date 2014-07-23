function n=getTriKRegVertex(k, ring, sector, offset)
% ReorderTriKRegComplex -- recursively retrieve k-regular vertex number
%  Usage
%    finen = ReorderTriKRegComplex(k, numRings, infineo)
%  Input
%    k           number of spokes in k-regular complex
%    ring        ring on which the vertex lies (center=0)
%    sector      sector in which the vertex lies (first sector=0)
%    offset      distance from start of sector (first position offset=0)
%  Output
%    n           integer, vertex number
% 

if (ring==0)
    n=1;
elseif (offset==0 && sector==0)
    n=0;
    if ring==1; n=1; end; % first building on center
    n=n+getTriKRegVertex(k, ring-1,0,0)+((ring-1)*k);
else
    n=getTriKRegVertex(k, ring,0,0)+mod(((sector*ring)+offset),k*ring);
end

% 
% Copyright (c) 2005. Kyle McDonald
% 