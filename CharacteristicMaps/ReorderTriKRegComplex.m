function rekmesh = ReorderTriKRegComplex(k,numRings,kmesh)
% ReorderTriKRegComplex -- Reorder vertices of a subdivided k-regular complex
%  Usage
%    rekmesh = ReorderTriKRegComplex(k, numRings, kmesh)
%  Input
%    k           number of spokes
%    numRings    number of rings in the coarse mesh
%    kmesh       subdivided k-regular complex mesh
%  Output
%    rekmesh     the reordered k-regular mesh
% 

global rekmesh kmesho OtoN NtoO n
n = size(kmesh, 2);
nv=size(kmesh{1},1);
OtoN=sparse(nv,nv);
NtoO=sparse(nv,nv);
kmesho=kmesh;
rekmesh=kmesh;

%% fix positions and faces
renameVertex(1,1); % set center
for r=1:numRings % then handle the rest
    for s=0:k-1
        for o=0:r-1
            % rename current vertex
            renameVertex(getTriKRegVertex(k,r,s,o),getTriKRegVertex(k,2*r,s,2*o));
            % rename associated triangle (some repetition)
            fixBetween(k, r,s,o, r,s,o+1);   % cur-next edge
            fixBetween(k, r,s,o, r-1,s,o);   % cur-inside edge
            fixBetween(k, r,s,o+1, r-1,s,o); % next-inside edge
        end
    end
end

%% fix NLs
rekmesh{3}=MakeTriKRegComplexNL(k, numRings*2); % generate/set new NLs
% if size(kmesh{1},3)>1 % if order > 0, rotate jets
%     if size(kmesh{1},3)==3, 
%         jetorder=1;
%     elseif size(kmesh{1},3)==6,
%         jetorder=2;
%     else
%         error('This code can be easily generalized to JSS of order > 2');
%     end
%     Flip = [0 1; 1 0];
    for i=1:nv
        rekmesh{3}{i} = OtoN([kmesho{3}{NtoO(i)}]);
%         o=OtoN([kmesho{3}{NtoO(i)}]); % reinterpret old NL
%         n=rekmesh{3}{i}; % new NL
%         k=length(o);
%         i1=find(n(1)==o)-1; i2=find(n(2)==o)-1;
%         cv=shiftdim(rekmesh{1}(i,:,:),1);
%         rekmesh{1}(i,:,:) = cv*S(inv([omega(0,k),omega(1,k)])*...
%                                  [omega(i1,k),omega(i2,k)],jetorder);
    end
% end

% fix the vertex between two old (coarse mesh) vertices
function fixBetween(k, r1,s1,o1, r2,s2,o2)
global kmesho
b=intersect( kmesho{3}{getTriKRegVertex(k,r1,s1,o1)},...
    kmesho{3}{getTriKRegVertex(k,r2,s2,o2)});
bn=getTriKRegVertex(k,r1+r2,s1,o1+o2);
renameVertex(b,bn);

% copy VL and fix FL
function renameVertex(old,new)
global rekmesh kmesho OtoN NtoO n;
OtoN(old)=new;
NtoO(new)=old;
rekmesh{1}(new,:)=kmesho{1}(old,:); % set vertex data
if n>3
    rekmesh{4}(new,:)=kmesho{4}(old,:); % set vector field data
end

rekmesh{2}(find(kmesho{2}==old))=new; % replace in FLs

% distance from v in NL1 to NL2
% 0 < off < length(NL)
function off=voffset(NL1,NL2,v)
off=mod(find(NL2==v)-find(NL1==v),length(NL1));

% 
% Copyright (c) 2005. Kyle McDonald
%
% Modified by Gang Xie, June 2007