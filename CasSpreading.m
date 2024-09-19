function Result_G = CasSpreading(A, H, alpha, beta, T, attack_num, DistributeMethod_Flag, AttachMethod_Flag)
 
degree = sum(A,2);
hyperdegree = sum(H,2);
% 定义负载
% alpha = 10;
% beta = 1;
L_init = alpha.*((degree + hyperdegree).^beta);
% L_init = alpha.*hyperdegree + degree.*beta;
L_real = L_init;
% 定义容量
gamma = T;
C_max = (1 + gamma) .* L_init;
% -----------方式定义-------------
% DistributeMethod_Flag = 1; % 1--均匀分配  2--度优先分配
% AttachMethod_Flag = 1; % 1--最大负载攻击 2--最小负载攻击 3--随机攻击（重复X次，取平均）
% 攻击S个节点
% attack_num = 5;
N = length(A);
M = size(H,2);
NodeState = ones(N,1); % 节点状态初始化 正常为1  正处于危机为0  危机后移除为-1
HyperEdgeState = ones(M,1); % 超边状态初始化
if( 1 == AttachMethod_Flag )
    % --------最大负载攻击--------
    [~,ld] = sort(L_init,"descend");
    attack_node = ld(1:attack_num);
elseif( 2 == AttachMethod_Flag )
    % --------最小负载攻击--------
    [~,ld] = sort(L_init,"ascend");
    attack_node = ld(1:attack_num);
elseif( 3 == AttachMethod_Flag )
    % ----------随机攻击----------
    ld = randperm(N);
    attack_node = ld(1:attack_num);
end
if( AttachMethod_Flag == 1 || AttachMethod_Flag == 2)
    NodeState(attack_node) = 0;
    AA = A;
    TimeStep = 1;
    while(sum(NodeState==0)~=0)
%         all_the_indexs = figure_NetEff( double(AA) ,length(AA));
%         Result_Max_Sub(TimeStep) = all_the_indexs.Max_sub_num/length(AA); %最大联通子图
%         Result_Net_Eff(TimeStep)= all_the_indexs.Eglob; % 网络效率
%         Result_G(TimeStep) = sum(UpdataHyperEdgeState(H,NodeState)) / size(H,2);
        Disuse_Node = find(NodeState==0); % 找出处于失效的节点
        Num_Disuse_Node = length(Disuse_Node); 
        % 更改处于失效的节点的状态、超边状态
        NodeState(Disuse_Node)=-1; 
        % 更新超边状态
        HyperEdgeState = UpdataHyperEdgeState(H,NodeState);
        % 找到失效节点的未失效超边，并将负载均分，再进一步分配给超边内的邻居节点
        for i = 1:Num_Disuse_Node
            % 负载分配策略->未失效超边均分
            % 找到失效节点的未失效超边
            HyEdges = find(H(Disuse_Node(i),:)~=0);
            undis_hyEdges = HyEdges(HyperEdgeState(HyEdges')~=0);
            % 找到处于失效的节点的邻居节点
            Neighbor_i = find_NeighborNode(AA,Disuse_Node(i)); % 找到处于失效的节点的邻居节点
            % 给每条超边分配的负载
            DetaL_To_EveHyEdge = L_real(Disuse_Node(i)) / length(undis_hyEdges);
            % 找到一条超边里属于邻居节点的那些
            for j = 1:length(undis_hyEdges)
                NodesInHyEdge = find(H(:,undis_hyEdges(j))~=0);
                % 找交集
                NeighNodesInHyEdge = intersect(NodesInHyEdge, Neighbor_i);
                % 剔除状态失效节点
                NeighNodesWithoutDisNode = [];
                for k = 1:length(NeighNodesInHyEdge)
                    if(NodeState(NeighNodesInHyEdge(k))~= -1)
                        NeighNodesWithoutDisNode = [NeighNodesWithoutDisNode;NeighNodesInHyEdge(k)];
                    end
                end
                % 分配负载给节点
                if( 1 == DistributeMethod_Flag)
                % ---------均匀分配------------
                    for k = 1:length(NeighNodesWithoutDisNode)
                        L_real(NeighNodesWithoutDisNode(k)) = L_real(NeighNodesWithoutDisNode(k)) + DetaL_To_EveHyEdge / length(NeighNodesWithoutDisNode);
                    end
                elseif( 2 == DistributeMethod_Flag)
                % ---------度优先分配------------
                    % 度之和
                    sum_deg = 0;
                    for k = 1:length(NeighNodesWithoutDisNode)
                        sum_deg = sum_deg + degree(NeighNodesWithoutDisNode(k));
                    end
                    % 按度分配
                    for k = 1:length(NeighNodesWithoutDisNode)
                        L_real(NeighNodesWithoutDisNode(k)) = L_real(NeighNodesWithoutDisNode(k)) + DetaL_To_EveHyEdge * degree(NeighNodesWithoutDisNode(k)) / sum_deg;
                    end
                end
            end
            
            % 更新当前节点状态 进行状态标记
            for k =1:length(AA)
                if NodeState(k) ~= -1
                    if L_real(k)> C_max(k) 
                        NodeState(k) = 0; 
                    end
                end
            end
        end
        % 将目标节点移除 置零
        AA(NodeState==-1,:)=0;
        AA(:,NodeState==-1)=0;
        % 更新度值
        degree = sum(AA,2);
        TimeStep = TimeStep + 1;
end
%     all_the_indexs = figure_NetEff( double(AA) ,length(AA));
%     Result_Max_Sub(TimeStep) = all_the_indexs.Max_sub_num/length(AA); %最大联通子图
%     Result_Net_Eff(TimeStep)= all_the_indexs.Eglob; % 网络效率
    Result_G = sum(UpdataHyperEdgeState(H,NodeState)) / size(H,2);
elseif(AttachMethod_Flag == 3)
    exec_time = 0;
    loop_time = 0.5;
    sum_G = 0;
    % ----------随机攻击----------
    while( exec_time <= loop_time)
        ld = randperm(N);
        attack_node = ld(1:attack_num);
        NodeState = ones(N,1); % 节点状态初始化 正常为1  正处于危机为0  危机后移除为-1
        HyperEdgeState = ones(M,1); % 超边状态初始化
        NodeState(attack_node) = 0;
        AA = A;
        TimeStep = 0;
        while(sum(NodeState==0)~=0)
%             all_the_indexs = figure_NetEff( double(AA) ,length(AA));
%             Result_Max_Sub(TimeStep) = all_the_indexs.Max_sub_num/length(AA); %最大联通子图
%             Result_Net_Eff(TimeStep)= all_the_indexs.Eglob; % 网络效率
%             Result_G(TimeStep) = sum(UpdataHyperEdgeState(H,NodeState)) / size(H,2);
            Disuse_Node = find(NodeState==0); % 找出处于失效的节点
            Num_Disuse_Node = length(Disuse_Node); 
            % 更改处于失效的节点的状态、超边状态
            NodeState(Disuse_Node)=-1; 
            % 更新超边状态
            HyperEdgeState = UpdataHyperEdgeState(H,NodeState);
            % 找到失效节点的未失效超边，并将负载均分，再进一步分配给超边内的邻居节点
            for i = 1:Num_Disuse_Node
                % 负载分配策略->未失效超边均分
                % 找到失效节点的未失效超边
                HyEdges = find(H(Disuse_Node(i),:)~=0);
                undis_hyEdges = HyEdges(HyperEdgeState(HyEdges')~=0);
                % 找到处于失效的节点的邻居节点
                Neighbor_i = find_NeighborNode(AA,Disuse_Node(i)); % 找到处于失效的节点的邻居节点
                % 给每条超边分配的负载
                DetaL_To_EveHyEdge = L_real(Disuse_Node(i)) / length(undis_hyEdges);
                % 找到一条超边里属于邻居节点的那些
                for j = 1:length(undis_hyEdges)
                    NodesInHyEdge = find(H(:,undis_hyEdges(j))~=0);
                    % 找交集
                    NeighNodesInHyEdge = intersect(NodesInHyEdge, Neighbor_i);
                    % 剔除状态失效节点
                    NeighNodesWithoutDisNode = [];
                    for k = 1:length(NeighNodesInHyEdge)
                        if(NodeState(NeighNodesInHyEdge(k))~= -1)
                            NeighNodesWithoutDisNode = [NeighNodesWithoutDisNode;NeighNodesInHyEdge(k)];
                        end
                    end
                    % 分配负载给节点
                    if( 1 == DistributeMethod_Flag)
                    % ---------均匀分配------------
                        for k = 1:length(NeighNodesWithoutDisNode)
                            L_real(NeighNodesWithoutDisNode(k)) = L_real(NeighNodesWithoutDisNode(k)) + DetaL_To_EveHyEdge / length(NeighNodesWithoutDisNode);
                        end
                    elseif( 2 == DistributeMethod_Flag)
                    % ---------度优先分配------------
                        % 度之和
                        sum_deg = 0;
                        for k = 1:length(NeighNodesWithoutDisNode)
                            sum_deg = sum_deg + degree(NeighNodesWithoutDisNode(k));
                        end
                        % 按度分配
                        for k = 1:length(NeighNodesWithoutDisNode)
                            L_real(NeighNodesWithoutDisNode(k)) = L_real(NeighNodesWithoutDisNode(k)) + DetaL_To_EveHyEdge * degree(NeighNodesWithoutDisNode(k)) / sum_deg;
                        end
                    end
                end
                
                % 更新当前节点状态 进行状态标记
                for k =1:length(AA)
                    if NodeState(k) ~= -1
                        if L_real(k)> C_max(k) 
                            NodeState(k) = 0; 
                        end
                    end
                end
            end
            % 将目标节点移除 置零
            AA(NodeState==-1,:)=0;
            AA(:,NodeState==-1)=0;
            % 更新度值
            degree = sum(AA,2);
            TimeStep = TimeStep + 1;
        end
%         all_the_indexs = figure_NetEff( double(AA) ,length(AA));
%         Result_Max_Sub(TimeStep) = all_the_indexs.Max_sub_num/length(AA); %最大联通子图
%         Result_Net_Eff(TimeStep)= all_the_indexs.Eglob; % 网络效率
        Result_G = sum(UpdataHyperEdgeState(H,NodeState)) / size(H,2);
        sum_G = sum_G + Result_G;
        exec_time = exec_time + 1;
    end
    Result_G = sum_G/ exec_time;
end