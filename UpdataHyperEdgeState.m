function newHyperEdgeState = UpdataHyperEdgeState(H, NodeState)
    
M = size(H,2);
newHyperEdgeState = ones(M,1);
for i=1:M
    if(isempty( find(NodeState(H(:,i)~=0)==1, 1) ))
        newHyperEdgeState(i) = 0;
    end
end
end
