%% 
% 3 Bucket model and Zipf distribution with three types of users and
% Freshness requirement as follows:
% Freshness=[F_a F_b F_c] where F_b=F_a*10^2 or some higher value.
% Setting is as follows:
% 
% 3 Bucket Uniform Distribution: 
% Prob_a, Prob_b=0.99*Prob_a, Prob_c=1-Prob_a-Prob_b
% Freshness requirement as above
% Number of Producers in each bucket are N_a=C,N_b=C,N_c=N-N_a-N_b where C
% is cache size and N is total number of producers.
% 
% Zipf Distribution:
% Use moderate Zipf parameter beta for probability distribution and use
% remaining setting as 3 Bucket Uniform distribution.
%% Probabilistic Save Implementation
clear all;
close all;
% clc;
%%
count=10^6;
Producers=100; % Number of Producers
global Pop_producers

global Freshness_requirment
const=10^2;
F_a=5;
F_b=const*F_a;
F_c=F_a;
Freshness_requirment=[F_a F_b F_c];

global Router1_hit_count

ProbForSavingVectorR1=1;%0.2:0.2:1.0;%1.0;
CacheSize=10:5:40;

Prob_a=0.25:0.05:0.45;
%% Least Expected Variables :::::::::::::::::::::::::::::::::::::::::::::::
R1_hit_count_Uni_LeastExpe=zeros(Producers,length(CacheSize),length(Prob_a));

R1_hit_count_Zipf_LeastExpe=zeros(Producers,length(CacheSize),length(beta));

N_min_3Bucket_LeastExpe_LeastExpe=zeros(length(CacheSize),length(Prob_a));
N_max_3Bucket_LeastExpe_LeastExpe=zeros(length(CacheSize),length(Prob_a));

N_min_Zipf_LeastExpe_LeastExpe=zeros(length(CacheSize),length(beta));
N_max_Zipf_LeastExpe_LeastExpe=zeros(length(CacheSize),length(beta));

%% LRU Variables ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
R1_hit_count_Uni_LRU=zeros(Producers,length(CacheSize),length(Prob_a));

R1_hit_count_Zipf_LRU=zeros(Producers,length(CacheSize),length(beta));

N_min_3Bucket_LRU_LRU=zeros(length(CacheSize),length(Prob_a));
N_max_3Bucket_LRU_LRU=zeros(length(CacheSize),length(Prob_a));

N_min_Zipf_LRU_LRU=zeros(length(CacheSize),length(beta));
N_max_Zipf_LRU_LRU=zeros(length(CacheSize),length(beta));
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

% Characteristic_time_Uni=zeros(Producers,length(CacheSize));
% Characteristic_time_Zipf=zeros(Producers,length(CacheSize));

global memoryR1_LeastExpe memoryR1_LRU Probability_producers

global count1 count2 % Checks cache is empty or not.

% Exponential inter-arrival time
time=cumsum(exprnd(1,count,1));

%% 3 Bucket Uniform Distribution

tic;
Prob_b=0.99*Prob_a;
Prob_c=ones(1,length(Prob_a))-Prob_a-Prob_b;
ProducersProbability_Uni=zeros(Producers,length(CacheSize),length(Prob_a));
for cache=1:length(CacheSize)
    N_a=CacheSize(cache);
    N_b=CacheSize(cache);
    N_c=Producers-N_a-N_b;
    ProducersProbability_Uni(:,cache,:)=[repmat(Prob_a./N_a,N_a,1);repmat(Prob_b./N_b,N_b,1);repmat(Prob_c./N_c,N_c,1)];
end
producersRequest_Uni=zeros(count,length(CacheSize),length(Prob_a));

for nn=1:length(Prob_a)
    for cache=1:length(CacheSize)
        r=rand(1,count);
        temp1(1,:)=ProducersProbability_Uni(:,cache,nn);
        for i=1:count
            producersRequest_Uni(i,cache,nn) = sum(r(1,i) >= cumsum([0, temp1]));
        end
    end
end
clear r;


%% 3 Bucket Uniform Distribution Least Expected

for nn=1:length(Prob_a)
    N_min=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    N_max=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    for cache=1:length(CacheSize)
        Probability_producers(1,:)=ProducersProbability_Uni(:,cache,nn);
        N_a=CacheSize(cache);
        N_b=CacheSize(cache);
        N_c=Producers-N_a-N_b;
        Pop_producers=[N_a N_b N_c];
        
        for jj=1:length(ProbForSavingVectorR1)
            ProbForSavingR1=ProbForSavingVectorR1(jj);
% memoryR1_LeastExpe and memoryR2_LeastExpe have following structure.
% First Column: Producer number
% Second Column: t_inst at which it is being fetched from producer

            memoryR1_LeastExpe=zeros(CacheSize(cache),2);

            Router1_hit_count=zeros(Producers,1);

            count1=0;
            count2=0;

            message=sprintf('Running for Cache Size=%d and ProbForSavingR1=%f Probability_a=%f'...
                            ,CacheSize(cache),ProbForSavingR1,Prob_a(nn));
            h=msgbox(message);
            clear message

            N_min_temp=0;
            N_max_temp=0;
            display('3 Bucket Uniform Distribution');
            for ii=1:length(time)
%                 display(t_inst);
                produ=producersRequest_Uni(ii,cache,nn);
                t_inst=time(ii);
%                 Frshness=Freshness(t_inst,1);
                [N_min_temp,N_max_temp]=router1_LeastExpe_plain_3Bucket(produ,t_inst,...
                                                            ProbForSavingR1,N_min_temp,N_max_temp);

            end
            N_min(jj,cache)=N_min_temp;
            N_max(jj,cache)=N_max_temp;

            delete(h);
            clear('h');
        end        

        R1_hit_count_Uni_LeastExpe(:,cache,nn)=Router1_hit_count;
    end
    N_min_3Bucket_LeastExpe_LeastExpe(:,nn)=N_min;
    N_max_3Bucket_LeastExpe_LeastExpe(:,nn)=N_max;
end
toc

% delete(h1);
% clear('h1');
display('Done Uniform');


requests_Uni=zeros(Producers,length(CacheSize),length(Prob_a));
for nn=1:length(Prob_a)
    for cache=1:length(CacheSize)
        for ii=1:count
            requests_Uni(producersRequest_Uni(ii,cache,nn),cache,nn)=requests_Uni(producersRequest_Uni(ii,cache,nn),cache,nn)+1;
        end
    end
end

% requests_Uni(1:3,1);
% R1_hit_count_Uni(1:3,1)
% temp1(:,1)=R1_hit_count(1:3,1);
% temp2(:,1)=requests(1:3,1);
% display('Hit rate Simulation for  2 Bucket Uniform Distribution')
% temp1./temp2
% temp2=repmat(requests_Uni,1,length(CacheSize));
hit_rate_Simul_Uni_LeastExpe=R1_hit_count_Uni_LeastExpe./requests_Uni;
clear temp1 temp2
hit_rate_total_Sim_Uni_LeastExpe=zeros(length(Prob_a),length(CacheSize));
for nn=1:length(Prob_a)
    temp1(:,:)=R1_hit_count_Uni_LeastExpe(:,:,nn);
    hit_rate_total_Sim_Uni_LeastExpe(nn,:)=sum(temp1)/count;
end

%% 3 Bucket Uniform Distribution LRU (Least Recently Used)

for nn=1:length(Prob_a)
    N_min=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    N_max=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    for cache=1:length(CacheSize)
        Probability_producers(1,:)=ProducersProbability_Uni(:,cache,nn);
        N_a=CacheSize(cache);
        N_b=CacheSize(cache);
        N_c=Producers-N_a-N_b;
        Pop_producers=[N_a N_b N_c];
        
        for jj=1:length(ProbForSavingVectorR1)
            ProbForSavingR1=ProbForSavingVectorR1(jj);
% memoryR1_LRU have following structure.
% First Column: t_stamp
% Second Column: Producer number
% Third Column: t_inst at which it is being fetched from producer
            memoryR1_LRU=zeros(CacheSize(cache),3);

            Router1_hit_count=zeros(Producers,1);

            message=sprintf('Running for Cache Size=%d and ProbForSavingR1=%f and Probability_a=%f'...
                            ,CacheSize(cache),ProbForSavingR1,Prob_a(nn));
            h=msgbox(message);
            clear message

            N_min_temp=0;
            N_max_temp=0;
            display('3 Bucket Uniform Distribution');
            for ii=1:length(time)
%                 display(t_inst);
                produ=producersRequest_Uni(ii,1);
                t_inst=time(ii);
%                 Frshness=Freshness(t_inst,1);
                [N_min_temp,N_max_temp]=router1_LRU_plain_3Bucket(produ,t_inst,...
                                                            ProbForSavingR1,N_min_temp,N_max_temp);%,Sele_1,Sele_2);

            end
            N_min(jj,cache)=N_min_temp;
            N_max(jj,cache)=N_max_temp;

            delete(h);
            clear('h');
        end

        R1_hit_count_Uni_LRU(:,cache,nn)=Router1_hit_count;
    end
    
    N_min_3Bucket_LRU_LRU(:,nn)=N_min;
    N_max_3Bucket_LRU_LRU(:,nn)=N_max;
end
toc


display('Done Uniform');

% requests_Uni=zeros(Producers,1);
% for ii=1:count
%     requests_Uni(producersRequest_Uni(ii))=requests_Uni(producersRequest_Uni(ii))+1;
% end

% requests_Uni(1:3,1);
% R1_hit_count_Uni(1:3,1)
% temp1(:,1)=R1_hit_count(1:3,1);
% temp2(:,1)=requests(1:3,1);
% display('Hit rate Simulation for  2 Bucket Uniform Distribution')
% temp1./temp2
% temp2=repmat(requests_Uni,1,length(CacheSize));
hit_rate_Simul_Uni_LRU=R1_hit_count_Uni_LRU./request_Uni;
clear temp1 temp2
hit_rate_total_Sim_Uni_LRU=zeros(length(Prob_a),length(CacheSize));
for nn=1:length(Prob_a)
    temp1(:,:)=R1_hit_count_Uni_LRU(:,:,nn);
    hit_rate_total_Sim_Uni_LRU(nn,:)=sum(temp1)/count;
end

%% Theorotical 2 Bucket Model
% for cache=1:length(CacheSize)
%     for jj=1:Producers
%         memory1=CacheSize(cache);
%         Prob(:,1)=ProducersProbability_Uni([1:jj-1 jj+1:Producers]);
%         f = @(t_c)parameterfunRAND(t_c,Prob,memory1);
%         x=fzero(f,1);
%         Characteristic_time_Uni(jj,cache)=x;
% %         clear Prob
%     end
% end
% 
% hit_rate_theor_Uni=zeros(Producers,length(CacheSize));
% for ii=1:3
%     hit_rate_theor_Uni(ii,:)=ones(1,Pop_producers)*(ProducersProbability_Uni(ii)*FreshnessMin)/(1+(ProducersProbability_Uni(ii)*FreshnessMin));
% end
% clear temp1 temp2;
% temp1=repmat(ProducersProbability_Uni',1,length(CacheSize));
% hit_rate_theor_Uni=(temp1.*Characteristic_time_Uni)./(ones(Producers,length(CacheSize))+(temp1.*Characteristic_time_Uni));
% 
% clear temp1;
% clear N_min N_inter N_max N_min_temp N_inter_temp N_max_temp jj kk t_inst cache produ i ii
% clear ProbForSavingR1 ProbForSavingR2
% 
% temp1=repmat(ProducersProbability_Uni',1,length(CacheSize));
% hit_rate_total_theor_Uni=sum(hit_rate_theor_Uni.*temp1)
tic
%% Zipf distribution with parameter beta
beta=0.5:.3:1.7;
nn=1:Producers;
ProducersProbability_Zipf=zeros(length(beta),Producers);
for ii=1:length(beta)
    ProducersProbability_Zipf(ii,:)=(nn.^-beta(ii))/sum((nn.^-beta(ii)));
end
producersRequest_Zipf=zeros(count,length(beta));

for jj=1:length(beta)
    r=rand(1,count);
    for ii=1:count
        temp1(1,:)=ProducersProbability_Zipf(jj,:);
        producersRequest_Zipf(ii,jj) = sum(r(1,i) >= cumsum([0, temp1]));
    end
end
clear r;

%% Zipf distribution with parameter beta Least Expected
for nn=1:length(beta)
    N_min=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    N_max=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    Probability_producers(1,:)=ProducersProbability_Zipf(nn,:);

    for cache=1:length(CacheSize)
        N_a=CacheSize(cache);
        N_b=CacheSize(cache);
        N_c=Producers-N_a-N_b;
        Pop_producers=[N_a N_b N_c];

        for jj=1:length(ProbForSavingVectorR1)
            ProbForSavingR1=ProbForSavingVectorR1(jj);
    % memoryR1_LeastExpe and memoryR2_LeastExpe have following structure.
    % First Column: latest time_instant when data was being used under
    % condition it was fresh.
    % Second Column: Producer number
    % Third Column: time_stamp at ehich data for corresponding producer was
    % being fetched and stored.
            memoryR1_LeastExpe=zeros(CacheSize(cache),2);

            Router1_hit_count=zeros(Producers,1);

            count1=0;
            count2=0;

            message=sprintf('Running for Cache Size=%d and ProbForSavingR1=%f and beta=%f'...
                            ,CacheSize(cache),ProbForSavingR1,beta(nn));
            h=msgbox(message);
            clear message

            N_min_temp=0;
            N_max_temp=0;
            display('Zipf Distribution');
            for ii=1:length(time)
    %                 display(t_inst);
                produ=producersRequest_Zipf(ii,nn);
                t_inst=time(ii);
    %                 Frshness=Freshness(t_inst,1);
                [N_min_temp,N_max_temp]=router1_LeastExpe_plain_3Bucket(produ,t_inst,...
                                                            ProbForSavingR1,N_min_temp,N_max_temp);%,Sele_1,Sele_2);

            end
            N_min(jj,cache)=N_min_temp; % jj-> row number; kk-> column number
            N_max(jj,cache)=N_max_temp;

            delete(h);
            clear('h');
        end

        R1_hit_count_Zipf_LeastExpe(:,cache,nn)=Router1_hit_count;
    end
    N_min_Zipf_LeastExpe_LeastExpe(:,nn)=N_min;
    N_max_Zipf_LeastExpe_LeastExpe(:,nn)=N_max;
end
toc
display('Done Zipf');

requests_Zipf=zeros(Producers,length(beta));
for nn=1:length(beta)
    for ii=1:count
        requests_Zipf(producersRequest_Zipf(ii,nn),nn)=requests_Zipf(producersRequest_Zipf(ii,nn))+1;
    end
end

% requests_Uni(1:3,1)
% R1_hit_count_Uni(1:3,1)
% temp1(:,1)=R1_hit_count(1:3,1);
% temp2(:,1)=requests(1:3,1);
% display('Hit rate Simulation for Zipf Distribution')
% temp1./temp2
% R1_hit_count_Zipf./requests_Zipf
temp2=repmat(requests_Zipf,1,1,length(CacheSize));
hit_rate_Simul_Zipf_LeastExpe=R1_hit_count_Zipf_LeastExpe./temp2;
clear temp1 temp2
hit_rate_total_Sim_Zipf_LeastExpe=zeros(length(beta),length(CacheSize));
for nn=1:length(beta)
    temp1(:,:)=R1_hit_count_Zipf_LeastExpe(:,:,nn);
    hit_rate_total_Sim_Zipf_LeastExpe(nn,:)=sum(temp1)/count;
end
%% Zipf distribution with parameter beta LRU
for cache=1:length(CacheSize)
    N_min=zeros(length(ProbForSavingVectorR1),length(ProbForSavingVectorR2));
    N_inter=zeros(length(ProbForSavingVectorR1),length(ProbForSavingVectorR2));
    N_max=zeros(length(ProbForSavingVectorR1),length(ProbForSavingVectorR2));
    for kk=1:length(ProbForSavingVectorR2)
        ProbForSavingR2=ProbForSavingVectorR2(kk);
        for jj=1:length(ProbForSavingVectorR1)
            ProbForSavingR1=ProbForSavingVectorR1(jj);
% memoryR1_LeastExpe and memoryR2_LeastExpe have following structure.
% First Column: latest time_instant when data was being used under
% condition it was fresh.
% Second Column: Producer number
% Third Column: time_stamp at ehich data for corresponding producer was
% being fetched and stored.
            memoryR1_LRU=zeros(CacheSize(cache),3);
            memoryR2_LRU=zeros(CacheSize(cache),3);

            Router1_hit_count=zeros(Producers,1);
            Router2_hit_count=zeros(Producers,1);

            message=sprintf('Running for Cache Size=%d and ProbForSavingR1=%f and ProbForSavingR2=%f'...
                            ,CacheSize(cache),ProbForSavingR1,ProbForSavingR2);
            h=msgbox(message);
            clear message

            N_min_temp=0;
            N_inter_temp=0;
            N_max_temp=0;
            display('Zipf Distribution');
            for ii=1:length(time)
%                 display(t_inst);
                produ=producersRequest_Zipf(ii,1);
                t_inst=time(ii);
%                 Frshness=Freshness(t_inst,1);
                [N_min_temp,N_inter_temp,N_max_temp]=router1_LRU_plain(produ,t_inst,...
                                                            ProbForSavingR1,ProbForSavingR2,...
                                                            N_min_temp,N_inter_temp,N_max_temp);%,Sele_1,Sele_2);

            end
            N_min(jj,kk)=N_min_temp; % jj-> row number; kk-> column number
            N_inter(jj,kk)=N_inter_temp;
            N_max(jj,kk)=N_max_temp;
            
       

% in two dimensional data N_min, N_inter and N_max, rows corresponding to
% data with respect ProbForSavingVectorR1 and columns corresponding to data
% with respect to ProbForSavingVectorR2.
            delete(h);
            clear('h');
        end
    end
        
    N_min_Zipf_LRU_LRU(cache,:,:)=N_min;
    N_inter_Zipf_LRU_LRU(cache,:,:)=N_inter;
    N_max_Zipf_LRU_LRU(cache,:,:)=N_max;
    
    R1_hit_count_Zipf_LRU(:,cache)=Router1_hit_count;
    R2_hit_count_Zipf_LRU(:,cache)=Router2_hit_count;
end
toc
display('Done Zipf');

% requests_Zipf=zeros(Producers,1);
% for ii=1:count
%     requests_Zipf(producersRequest_Zipf(ii))=requests_Zipf(producersRequest_Zipf(ii))+1;
% end

% requests_Uni(1:3,1)
% R1_hit_count_Uni(1:3,1)
% temp1(:,1)=R1_hit_count(1:3,1);
% temp2(:,1)=requests(1:3,1);
% display('Hit rate Simulation for Zipf Distribution')
% temp1./temp2
% R1_hit_count_Zipf./requests_Zipf
temp2=repmat(requests_Zipf,1,length(CacheSize));
hit_rate_Simul_Zipf_LRU=R1_hit_count_Zipf_LRU./temp2;
% clear temp1 temp2

hit_rate_total_Sim_Zipf_LRU=sum(R1_hit_count_Zipf_LRU)/count
%% Theorotical calculations
% for cache=1:length(CacheSize)
%     ProducersProbability_Zipf_sele(1,:)=ProducersProbability_Zipf(1:Producers);
%     for jj=1:Producers
%         memory1=CacheSize(cache);        
%         Prob(:,1)=ProducersProbability_Zipf_sele([1:jj-1 jj+1:Producers]);
%         f = @(t_c)parameterfunLeastExpe(t_c,Prob,memory1);
%         x=fzero(f,1);
%         Characteristic_time_Zipf(jj,cache)=x;
%         clear Prob
%     end
% end    
% temp1=repmat(ProducersProbability_Zipf',1,length(CacheSize));
% hit_rate_theor_Zipf=(temp1.*Characteristic_time_Zipf)./(ones(Producers,length(CacheSize))+(temp1.*Characteristic_time_Zipf));
% 
% clear temp1;
% clear memoryR1_LeastExpe memoryR2_LeastExpe N_min N_inter N_max N_min_temp N_inter_temp N_max_temp jj kk t_inst cache produ i ii
% clear ProbForSavingR1 ProbForSavingR2
% 
% temp1=repmat(ProducersProbability_Zipf',1,length(CacheSize));
% hit_rate_total_theor_Zipf=sum(hit_rate_theor_Zipf.*temp1)
%% always change the dataname for saving. Keep it simple and discriptive.
% temp1=cd;
cd('D:\IoT\IoT\31Jan\Least Expected\Data')
save('cmp_least_expected_LRU');
