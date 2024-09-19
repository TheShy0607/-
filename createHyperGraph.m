clear;clc;tic;
% ��ʼ�ڵ���
m0 = 5;
% �����ڵ���
m = m0 - 1;
% ���Ʋ���
M = 100; % ������
N = m0 +  (m0-1)*(M-1); % �ڵ���

%% ��ʼ���ӳ�ͼ
% ��Ҫһ�����������ʾ�ڵ�ͳ��ߵĹ�ϵ�����ܱ�ʾ�ڵ������ĸ�����
%%%%%%%%%%% ��ʼ��
for i=1:m0
    edge{i} = i;
end
hyperedges{1,1} = edge;
hyperedges_num = 2;
node_num = m0;
%%%%%%%%%%% 1--������� or 2--������������
connect_flag = 1;
%%%%%%%%%%%
if(1 == connect_flag)
    while(hyperedges_num <= M && node_num <= N)
        % �������
        % �ȸ������� 1/num ���ȡһ��
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
    % ������������
    while(hyperedges_num <= M && node_num <= N)
        for i=1:m
            edge{i} = node_num + i;
        end
        % ��Ҫ����ÿ���ڵ�ĳ���
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

%%%%%%%%%%%% ת��Ϊ��������
% ��Ҫ����ÿ���ڵ�ĳ���
for i=1:length(hyperedges)
    hyperedge = hyperedges{i,1};
    for j=1:length(hyperedge)
        H(hyperedge{j},i) = 1;
    end
end

%% ����������
% 1--��� 2--��ȫ 3--��������
node_num = size(H,1);
hyperedges_num = size(H,2);
A = zeros(node_num,node_num);
%%%%%%%%%%%%%%%%%
connect_flag = 1;
% �������
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
% ��ȫ����
elseif(2 == connect_flag)
    for i=1:hyperedges_num
        for j=1:length(hyperedges{i,1})
            for k=1+j:length(hyperedges{i,1})
                A(hyperedges{i,1}{j}, hyperedges{i,1}{k}) = 1;
                A(hyperedges{i,1}{k}, hyperedges{i,1}{j}) = 1;
            end
        end
    end
% ������������
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
                    % �ڴ˴����ݽڵ����Ϣ��������
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
% ���Ƹ��ʷֲ�
figure;  % ����һ���µ�ͼ�δ���
histogram(degree, 'Normalization', 'probability');  % ����ֱ��ͼ������ֱ��ͼ��һ��Ϊ���ʷֲ�
title('�ȸ��ʷֲ�','FontName',"����");  % ��ӱ���
xlabel('��','FontName',"����");  % ���x���ǩ
ylabel('����','FontName',"����");  % ���y���ǩ
figure;  % ����һ���µ�ͼ�δ���
histogram(hyperdegree, 'Normalization', 'probability');  % ����ֱ��ͼ������ֱ��ͼ��һ��Ϊ���ʷֲ�
title('���ȸ��ʷֲ�','FontName',"����");  % ��ӱ���
xlabel('����','FontName',"����");  % ���x���ǩ
ylabel('����','FontName',"����");  % ���y���ǩ