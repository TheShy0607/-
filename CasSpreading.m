function Result_G = CasSpreading(A, H, alpha, beta, T, attack_num, DistributeMethod_Flag, AttachMethod_Flag)
 
degree = sum(A,2);
hyperdegree = sum(H,2);
% ���帺��
% alpha = 10;
% beta = 1;
L_init = alpha.*((degree + hyperdegree).^beta);
% L_init = alpha.*hyperdegree + degree.*beta;
L_real = L_init;
% ��������
gamma = T;
C_max = (1 + gamma) .* L_init;
% -----------��ʽ����-------------
% DistributeMethod_Flag = 1; % 1--���ȷ���  2--�����ȷ���
% AttachMethod_Flag = 1; % 1--����ع��� 2--��С���ع��� 3--����������ظ�X�Σ�ȡƽ����
% ����S���ڵ�
% attack_num = 5;
N = length(A);
M = size(H,2);
NodeState = ones(N,1); % �ڵ�״̬��ʼ�� ����Ϊ1  ������Σ��Ϊ0  Σ�����Ƴ�Ϊ-1
HyperEdgeState = ones(M,1); % ����״̬��ʼ��
if( 1 == AttachMethod_Flag )
    % --------����ع���--------
    [~,ld] = sort(L_init,"descend");
    attack_node = ld(1:attack_num);
elseif( 2 == AttachMethod_Flag )
    % --------��С���ع���--------
    [~,ld] = sort(L_init,"ascend");
    attack_node = ld(1:attack_num);
elseif( 3 == AttachMethod_Flag )
    % ----------�������----------
    ld = randperm(N);
    attack_node = ld(1:attack_num);
end
if( AttachMethod_Flag == 1 || AttachMethod_Flag == 2)
    NodeState(attack_node) = 0;
    AA = A;
    TimeStep = 1;
    while(sum(NodeState==0)~=0)
%         all_the_indexs = figure_NetEff( double(AA) ,length(AA));
%         Result_Max_Sub(TimeStep) = all_the_indexs.Max_sub_num/length(AA); %�����ͨ��ͼ
%         Result_Net_Eff(TimeStep)= all_the_indexs.Eglob; % ����Ч��
%         Result_G(TimeStep) = sum(UpdataHyperEdgeState(H,NodeState)) / size(H,2);
        Disuse_Node = find(NodeState==0); % �ҳ�����ʧЧ�Ľڵ�
        Num_Disuse_Node = length(Disuse_Node); 
        % ���Ĵ���ʧЧ�Ľڵ��״̬������״̬
        NodeState(Disuse_Node)=-1; 
        % ���³���״̬
        HyperEdgeState = UpdataHyperEdgeState(H,NodeState);
        % �ҵ�ʧЧ�ڵ��δʧЧ���ߣ��������ؾ��֣��ٽ�һ������������ڵ��ھӽڵ�
        for i = 1:Num_Disuse_Node
            % ���ط������->δʧЧ���߾���
            % �ҵ�ʧЧ�ڵ��δʧЧ����
            HyEdges = find(H(Disuse_Node(i),:)~=0);
            undis_hyEdges = HyEdges(HyperEdgeState(HyEdges')~=0);
            % �ҵ�����ʧЧ�Ľڵ���ھӽڵ�
            Neighbor_i = find_NeighborNode(AA,Disuse_Node(i)); % �ҵ�����ʧЧ�Ľڵ���ھӽڵ�
            % ��ÿ�����߷���ĸ���
            DetaL_To_EveHyEdge = L_real(Disuse_Node(i)) / length(undis_hyEdges);
            % �ҵ�һ�������������ھӽڵ����Щ
            for j = 1:length(undis_hyEdges)
                NodesInHyEdge = find(H(:,undis_hyEdges(j))~=0);
                % �ҽ���
                NeighNodesInHyEdge = intersect(NodesInHyEdge, Neighbor_i);
                % �޳�״̬ʧЧ�ڵ�
                NeighNodesWithoutDisNode = [];
                for k = 1:length(NeighNodesInHyEdge)
                    if(NodeState(NeighNodesInHyEdge(k))~= -1)
                        NeighNodesWithoutDisNode = [NeighNodesWithoutDisNode;NeighNodesInHyEdge(k)];
                    end
                end
                % ���为�ظ��ڵ�
                if( 1 == DistributeMethod_Flag)
                % ---------���ȷ���------------
                    for k = 1:length(NeighNodesWithoutDisNode)
                        L_real(NeighNodesWithoutDisNode(k)) = L_real(NeighNodesWithoutDisNode(k)) + DetaL_To_EveHyEdge / length(NeighNodesWithoutDisNode);
                    end
                elseif( 2 == DistributeMethod_Flag)
                % ---------�����ȷ���------------
                    % ��֮��
                    sum_deg = 0;
                    for k = 1:length(NeighNodesWithoutDisNode)
                        sum_deg = sum_deg + degree(NeighNodesWithoutDisNode(k));
                    end
                    % ���ȷ���
                    for k = 1:length(NeighNodesWithoutDisNode)
                        L_real(NeighNodesWithoutDisNode(k)) = L_real(NeighNodesWithoutDisNode(k)) + DetaL_To_EveHyEdge * degree(NeighNodesWithoutDisNode(k)) / sum_deg;
                    end
                end
            end
            
            % ���µ�ǰ�ڵ�״̬ ����״̬���
            for k =1:length(AA)
                if NodeState(k) ~= -1
                    if L_real(k)> C_max(k) 
                        NodeState(k) = 0; 
                    end
                end
            end
        end
        % ��Ŀ��ڵ��Ƴ� ����
        AA(NodeState==-1,:)=0;
        AA(:,NodeState==-1)=0;
        % ���¶�ֵ
        degree = sum(AA,2);
        TimeStep = TimeStep + 1;
end
%     all_the_indexs = figure_NetEff( double(AA) ,length(AA));
%     Result_Max_Sub(TimeStep) = all_the_indexs.Max_sub_num/length(AA); %�����ͨ��ͼ
%     Result_Net_Eff(TimeStep)= all_the_indexs.Eglob; % ����Ч��
    Result_G = sum(UpdataHyperEdgeState(H,NodeState)) / size(H,2);
elseif(AttachMethod_Flag == 3)
    exec_time = 0;
    loop_time = 0.5;
    sum_G = 0;
    % ----------�������----------
    while( exec_time <= loop_time)
        ld = randperm(N);
        attack_node = ld(1:attack_num);
        NodeState = ones(N,1); % �ڵ�״̬��ʼ�� ����Ϊ1  ������Σ��Ϊ0  Σ�����Ƴ�Ϊ-1
        HyperEdgeState = ones(M,1); % ����״̬��ʼ��
        NodeState(attack_node) = 0;
        AA = A;
        TimeStep = 0;
        while(sum(NodeState==0)~=0)
%             all_the_indexs = figure_NetEff( double(AA) ,length(AA));
%             Result_Max_Sub(TimeStep) = all_the_indexs.Max_sub_num/length(AA); %�����ͨ��ͼ
%             Result_Net_Eff(TimeStep)= all_the_indexs.Eglob; % ����Ч��
%             Result_G(TimeStep) = sum(UpdataHyperEdgeState(H,NodeState)) / size(H,2);
            Disuse_Node = find(NodeState==0); % �ҳ�����ʧЧ�Ľڵ�
            Num_Disuse_Node = length(Disuse_Node); 
            % ���Ĵ���ʧЧ�Ľڵ��״̬������״̬
            NodeState(Disuse_Node)=-1; 
            % ���³���״̬
            HyperEdgeState = UpdataHyperEdgeState(H,NodeState);
            % �ҵ�ʧЧ�ڵ��δʧЧ���ߣ��������ؾ��֣��ٽ�һ������������ڵ��ھӽڵ�
            for i = 1:Num_Disuse_Node
                % ���ط������->δʧЧ���߾���
                % �ҵ�ʧЧ�ڵ��δʧЧ����
                HyEdges = find(H(Disuse_Node(i),:)~=0);
                undis_hyEdges = HyEdges(HyperEdgeState(HyEdges')~=0);
                % �ҵ�����ʧЧ�Ľڵ���ھӽڵ�
                Neighbor_i = find_NeighborNode(AA,Disuse_Node(i)); % �ҵ�����ʧЧ�Ľڵ���ھӽڵ�
                % ��ÿ�����߷���ĸ���
                DetaL_To_EveHyEdge = L_real(Disuse_Node(i)) / length(undis_hyEdges);
                % �ҵ�һ�������������ھӽڵ����Щ
                for j = 1:length(undis_hyEdges)
                    NodesInHyEdge = find(H(:,undis_hyEdges(j))~=0);
                    % �ҽ���
                    NeighNodesInHyEdge = intersect(NodesInHyEdge, Neighbor_i);
                    % �޳�״̬ʧЧ�ڵ�
                    NeighNodesWithoutDisNode = [];
                    for k = 1:length(NeighNodesInHyEdge)
                        if(NodeState(NeighNodesInHyEdge(k))~= -1)
                            NeighNodesWithoutDisNode = [NeighNodesWithoutDisNode;NeighNodesInHyEdge(k)];
                        end
                    end
                    % ���为�ظ��ڵ�
                    if( 1 == DistributeMethod_Flag)
                    % ---------���ȷ���------------
                        for k = 1:length(NeighNodesWithoutDisNode)
                            L_real(NeighNodesWithoutDisNode(k)) = L_real(NeighNodesWithoutDisNode(k)) + DetaL_To_EveHyEdge / length(NeighNodesWithoutDisNode);
                        end
                    elseif( 2 == DistributeMethod_Flag)
                    % ---------�����ȷ���------------
                        % ��֮��
                        sum_deg = 0;
                        for k = 1:length(NeighNodesWithoutDisNode)
                            sum_deg = sum_deg + degree(NeighNodesWithoutDisNode(k));
                        end
                        % ���ȷ���
                        for k = 1:length(NeighNodesWithoutDisNode)
                            L_real(NeighNodesWithoutDisNode(k)) = L_real(NeighNodesWithoutDisNode(k)) + DetaL_To_EveHyEdge * degree(NeighNodesWithoutDisNode(k)) / sum_deg;
                        end
                    end
                end
                
                % ���µ�ǰ�ڵ�״̬ ����״̬���
                for k =1:length(AA)
                    if NodeState(k) ~= -1
                        if L_real(k)> C_max(k) 
                            NodeState(k) = 0; 
                        end
                    end
                end
            end
            % ��Ŀ��ڵ��Ƴ� ����
            AA(NodeState==-1,:)=0;
            AA(:,NodeState==-1)=0;
            % ���¶�ֵ
            degree = sum(AA,2);
            TimeStep = TimeStep + 1;
        end
%         all_the_indexs = figure_NetEff( double(AA) ,length(AA));
%         Result_Max_Sub(TimeStep) = all_the_indexs.Max_sub_num/length(AA); %�����ͨ��ͼ
%         Result_Net_Eff(TimeStep)= all_the_indexs.Eglob; % ����Ч��
        Result_G = sum(UpdataHyperEdgeState(H,NodeState)) / size(H,2);
        sum_G = sum_G + Result_G;
        exec_time = exec_time + 1;
    end
    Result_G = sum_G/ exec_time;
end