%找到i的邻居节点
function NeighborNode = find_NeighborNode(A,i)
NeighborNode = find(A(i,:)~=0);
end
