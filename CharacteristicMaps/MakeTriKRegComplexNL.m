function NL=MakeTriKRegComplexNL(ink,numRings)
% MakeTriKRegComplexNL -- Make the neighbor lists for a k-regular complex
%  Usage
%    NL = MakeTriKRegComplexNL(ink, numRings)
%  Input
%    ink        number of spokes
%    numRings   number of rings
%
%  Output
%    NL   A cell of matrices, index by vertex number,
%         each matrix corresponding to the neighbor list of that vertex
% 

global k max
k=ink;
max=getTriKRegVertex(k,numRings+1,0,0)-1;

NL{1}=(1:k)+1; % center
for r=1:numRings
    for s=0:k-1
        for o=0:r-1
            cv=getTriKRegVertex(k,r,s,o);
            if(o==0) % spoke
                if r==numRings  % boundary
                    NL{cv}=[v(r,s,1), v(r-1,s,0), v(r,s,-1)];
                else
                    NL{cv}=[
                        v(r-1, s,  0),...
                            v(r+0, s, -1),...
                            v(r+1, s, -1),...
                            v(r+1, s,  0),...
                            v(r+1, s, +1),...
                            v(r+0, s, +1)];
                end
            else % bead
                if r==numRings % boundary
                    NL{cv}=[v(r,s,o+1), v(r-1,s,o), v(r-1,s,o-1), v(r,s,o-1)];
                else
                    NL{cv}=[
                        v(r+0, s, o-1),...
                            v(r+1, s, o+0),...
                            v(r+1, s, o+1),...
                            v(r+0, s, o+1),...
                            v(r-1, s, o+0),...
                            v(r-1, s, o-1)];
                end
            end
        end
    end
end



function vertex=v(r,s,o)
global k max
vertex=getTriKRegVertex(k,r,s,o);
if(vertex>max) vertex=[]; end; % out of range vertex

% 
% Copyright (c) 2005. Kyle McDonald
%
%  Modified by Gang Xie, May 2007