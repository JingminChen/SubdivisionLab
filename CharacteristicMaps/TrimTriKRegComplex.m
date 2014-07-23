function trimmedkmesh = TrimTriKRegComplex(kmesh,k,numRings,ChangeVecField)
% TrimTriKRegComplex -- Trim a k-regular complex
%  Usage
%    trimmedkmesh = TrimTriKRegComplex(kmesh,k,numRings)
%  Input
%    k             number of spokes
%    numRings      number of rings to retain
%    kmesh         k-regular complex mesh, with at least (numRings) rings
%  Output
%    trimmedkmesh  the trimmed k-regular mesh
% 
if nargin==3
    ChangeVecField = 1;
end
n=size(kmesh,2);
max=getTriKRegVertex(k,numRings+1,0,0)-1;

trimmedkmesh{1} = kmesh{1}(1:max,:,:);
if n>3
    trimmedkmesh{4} = kmesh{4}(1:max,:,:);
    trimmedkmesh{5} = kmesh{5};
end


trimmedkmesh{2} = kmesh{2};
[I,J]=find(kmesh{2}>max);
trimmedkmesh{2}(I,:) = [];

trimmedkmesh{3} = kmesh{3}(1:max);
for s = 0:k-1
    for o = 0:numRings-1
        if o==0
            NLvn = [getTriKRegVertex(k,numRings,s,1), getTriKRegVertex(k,numRings-1,s,0), getTriKRegVertex(k,numRings,s,-1)];
        else
            NLvn = [getTriKRegVertex(k,numRings,s,o+1), getTriKRegVertex(k,numRings-1,s,o), getTriKRegVertex(k,numRings-1,s,o-1), getTriKRegVertex(k,numRings,s,o-1)];
        end
        cv = getTriKRegVertex(k,numRings,s,o);
        trimmedkmesh{3}{cv} = NLvn;
        if n>3 & ChangeVecField==1
            NLvo = kmesh{3}{cv};
            i = find(NLvo==NLvn(1))-1;
            theta = i*pi/3;
            NLvo = [NLvo NLvo(1)];
            if NLvn(2)==NLvo(i+2)
                %same order
                trimmedkmesh{4}(cv,:) = kmesh{4}(cv,:)*[cos(theta) sin(theta); -sin(theta) cos(theta)]';
            else
                %opposite order
                trimmedkmesh{4}(cv,:) = kmesh{4}(cv,:)*[cos(theta) sin(theta); sin(theta) -cos(theta)]';
            end
        end
    end
end
% for i=1:max,
%     trimmedkmesh{3}{i}(find(trimmedkmesh{3}{i}>max))=[]; % erase out of bound vertices
% end


% 
% Copyright (c) 2005. Thomas P.-Y. Yu
%
% Modified by Gang Xie, June 2007