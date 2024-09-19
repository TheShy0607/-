clear;clc;tic;
% 初始节点数
m0 = 5;
% 新增节点数
m = m0 - 1;
% 限制参数
M = 100; % 超边数
N = m0 +  (m0-1)*(M-1); % 节点数

%% 开始连接超图
% 需要一个关联矩阵表示节点和超边的关系，即能表示节点属于哪个超边
%%%%%%%%%%% 初始化
for i=1:m0
    edge{i} = i;
end
hyperedges{1,1} = edge;
hyperedges_num = 2;
node_num = m0;
%%%%%%%%%%% 1--随机连接 or 2--超度优先连接
connect_flag = 1;
%%%%%%%%%%%
if(1 == connect_flag)
    while(hyperedges_num <= M && node_num <= N)
        % 随机连接
        % 等概率连接 1/num 随机取一个
        random_rank = randperm(node_num);
        for i=1:m
            edge{i} = node_num + i;
        end
        edge{m0} = random_rank(1);
        
        hyperedges{hyperedges_num,1} = edge;
        hyperedges_num = hyperedges_num + 1;
        node_num = node_num + m;
    end
    
elseif(2 == connect_flag)
    % 超度优先连接
    while(hyperedges_num <= M && node_num <= N)
        for i=1:m
            edge{i} = node_num + i;
        end
        % 需要计算每个节点的超度
        for i=1:length(hyperedges)
            hyperedge = hyperedges{i,1};
            for j=1:length(hyperedge)
                H(hyperedge{j},i) = 1;
            end
        end
        HyperDegree = sum(H,2);
        choice_seq = [];
        for i=1:length(HyperDegree)
            for j=1:HyperDegree(i)
                choice_seq = [choice_seq i];
            end
        end
        random_rank = randperm(length(choice_seq));
        edge{m0} = choice_seq(random_rank(1));

        hyperedges{hyperedges_num,1} = edge;
        hyperedges_num = hyperedges_num + 1;
        node_num = node_num + m;
    end
end

%%%%%%%%%%%% 转换为关联矩阵
% 需要计算每个节点的超度
for i=1:length(hyperedges)
    hyperedge = hyperedges{i,1};
    for j=1:length(hyperedge)
        H(hyperedge{j},i) = 1;
    end
end

%% 超边内连接
% 1--随机 2--完全 3--超度优先
node_num = size(H,1);
hyperedges_num = size(H,2);
A = zeros(node_num,node_num);
%%%%%%%%%%%%%%%%%
connect_flag = 1;
% 随机连接
if(1 == connect_flag)
    prob = 0.5;
    while(node_num ~= length(bfsearch(graph(A),1)))
        for i=1:hyperedges_num
            for j=1:length(hyperedges{i,1})
                degree = sum(A);
                    for k=1:length(hyperedges{i,1})
                        if( rand() < prob && hyperedges{i,1}{j}~=hyperedges{i,1}{k})
                            A(hyperedges{i,1}{j}, hyperedges{i,1}{k}) = 1;
                            A(hyperedges{i,1}{k}, hyperedges{i,1}{j}) = 1;
                        end
                    end
            end
        end
    end
% 完全连接
elseif(2 == connect_flag)
    for i=1:hyperedges_num
        for j=1:length(hyperedges{i,1})
            for k=1+j:length(hyperedges{i,1})
                A(hyperedges{i,1}{j}, hyperedges{i,1}{k}) = 1;
                A(hyperedges{i,1}{k}, hyperedges{i,1}{j}) = 1;
            end
        end
    end
% 超度优先连接
elseif(3 == connect_flag)
    while(node_num ~= length(bfsearch(graph(A),1)))
        for i=1:hyperedges_num
            for j=1:length(hyperedges{i,1})
                hyperdegree = sum(H,2);
                sum_hd = 0;
                prob = [];
                for jj=1:length(hyperedges{i,1})
                    sum_hd = sum_hd + hyperdegree(hyperedges{i,1}{jj});
                    prob = [prob; hyperdegree(hyperedges{i,1}{jj})];
                end
                prob = prob ./ sum_hd;
                for k=1:length(hyperedges{i,1})
                    % 在此处根据节点度信息分配连接
                    if rand() < prob(k) && hyperedges{i,1}{j} ~= hyperedges{i,1}{k}
                        A(hyperedges{i,1}{j}, hyperedges{i,1}{k}) = 1;
                        A(hyperedges{i,1}{k}, hyperedges{i,1}{j}) = 1;
                    end
                end
            end
        end
    end

end

h = plot(graph(A));
layout(h,'force','UseGravity',true)
title_text = ['K = ' num2str(m0)];
title(title_text);
degree = sum(A,2);
hyperdegree = sum(H,2);
PlotDegDistribution(degree, hyperdegree);
title(title_text);
% 
% 绘制概率分布
figure;  % 创建一个新的图形窗口
histogram(degree, 'Normalization', 'probability');  % 创建直方图，并将直方图归一化为概率分布
title('度概率分布','FontName',"宋体");  % 添加标题
xlabel('度','FontName',"宋体");  % 添加x轴标签
ylabel('概率','FontName',"宋体");  % 添加y轴标签
figure;  % 创建一个新的图形窗口
histogram(hyperdegree, 'Normalization', 'probability');  % 创建直方图，并将直方图归一化为概率分布
title('超度概率分布','FontName',"宋体");  % 添加标题
xlabel('超度','FontName',"宋体");  % 添加x轴标签
ylabel('概率','FontName',"宋体");  % 添加y轴标签