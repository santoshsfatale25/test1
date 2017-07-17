%% 
% Program to simulate LRU-LRU policy at both routers. Producers are divide 
% into two bucket 1 and 2 such that bucket 1 producers are uniformlly high 
% popular and bucket 2 producers are uniformlly low popular. Simulation is
% done for various probabiliy for saving, cache size. Two producer type
% have two diffrent freshness parameters.

%% Probabilistic Save Implementation
clear all;
close all;
% clc;
%%
count=10^6;
Producers=50; % Number of Producers
Pop_producers=3;

% Freshness= floor(rand(count,1)*(Frshness_max-Frshness_min))+1;
% Freshness considered to be more than cache size. So this simulation is
% with Freshness 30 and cache size as 10 and 15 for Zipf distribution.

global FreshnessMin FreshnessMax
FreshnessMin= 10000;
FreshnessMax=10000;

global Router1_hit_count Router2_hit_count

ProbForSavingVectorR2=1;%0.2:0.2:1.0;%0.0:0.02:0.08;
ProbForSavingVectorR1=1;%0.2:0.2:1.0;%1.0;
CacheSize=10:5:20;

R1_hit_count_Uni=zeros(Producers,length(CacheSize));
R2_hit_count_Uni=zeros(Producers,length(CacheSize));

R1_hit_count_Zipf=zeros(Producers,length(CacheSize));
R2_hit_count_Zipf=zeros(Producers,length(CacheSize));

N_min_2Bucket_FIFO_FIFO=zeros(length(CacheSize),length(ProbForSavingVectorR1),length(ProbForSavingVectorR2));
N_inter_2Bucket_FIFO_FIFO=zeros(length(CacheSize),length(ProbForSavingVectorR1),length(ProbForSavingVectorR2));
N_max_2Bucket_FIFO_FIFO=zeros(length(CacheSize),length(ProbForSavingVectorR1),length(ProbForSavingVectorR2));

N_min_Zipf_FIFO_FIFO=zeros(length(CacheSize),length(ProbForSavingVectorR1),length(ProbForSavingVectorR2));
N_inter_Zipf_FIFO_FIFO=zeros(length(CacheSize),length(ProbForSavingVectorR1),length(ProbForSavingVectorR2));
N_max_Zipf_FIFO_FIFO=zeros(length(CacheSize),length(ProbForSavingVectorR1),length(ProbForSavingVectorR2));

Characteristic_time_Uni=zeros(Producers-Pop_producers,length(CacheSize));
Characteristic_time_Zipf=zeros(Producers-Pop_producers,length(CacheSize));

global memoryR1_FIFO %Frequency_R1 % memory for cache: C_1-> producer;C3->t_stamp

global memoryR2_FIFO %Frequency_R2 % memory for cache: C_1-> producer;C3->t_stamp

global count1 count2 % Checks cache is empty or not.
count1=0;
count2=0;

%% 2 Bucket Uniform Distribution

tic

ProducersProbability_Uni=[ones(1,3)*0.25 ones(1,Producers-3)*0.25/(Producers-3)];
producersRequest_Uni=zeros(count,1);
r=rand(1,count);
for i=1:count
    producersRequest_Uni(i,1) = sum(r(1,i) >= cumsum([0, ProducersProbability_Uni]));
end
clear r;


% Exponential inter-arrival time
time=cumsum(exprnd(1,count,1));

for cache=1:length(CacheSize)
    N_min=zeros(length(ProbForSavingVectorR1),length(ProbForSavingVectorR2));
    N_inter=zeros(length(ProbForSavingVectorR1),length(ProbForSavingVectorR2));
    N_max=zeros(length(ProbForSavingVectorR1),length(ProbForSavingVectorR2));
    for kk=1:length(ProbForSavingVectorR2)
        ProbForSavingR2=ProbForSavingVectorR2(kk);
        for jj=1:length(ProbForSavingVectorR1)
            ProbForSavingR1=ProbForSavingVectorR1(jj);
% memoryR1_FIFO and memoryR2_FIFO have following structure.
% First Column: Producer number
% Second Column: t_inst at which it is being fetched from producer

            memoryR1_FIFO=zeros(CacheSize(cache),2);
            memoryR2_FIFO=zeros(CacheSize(cache),2);

            Router1_hit_count=zeros(Producers,1);
            Router2_hit_count=zeros(Producers,1);

            message=sprintf('Running for Cache Size=%d and ProbForSavingR1=%f and ProbForSavingR2=%f'...
                            ,CacheSize(cache),ProbForSavingR1,ProbForSavingR2);
            h=msgbox(message);
            clear message

            N_min_temp=0;
            N_inter_temp=0;
            N_max_temp=0;
            display('2 Bucket Uniform Distribution');
            for ii=1:length(time)
%                 display(t_inst);
                produ=producersRequest_Uni(ii,1);
                t_inst=time(ii);
%                 Frshness=Freshness(t_inst,1);
                [N_min_temp,N_inter_temp,N_max_temp]=router1_FIFO_plain(produ,t_inst,...
                                                            ProbForSavingR1,ProbForSavingR2,Pop_producers,...
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
    % To calculate Characteristic time T_C

    N_min_2Bucket_FIFO_FIFO(cache,:,:)=N_min;
    N_inter_2Bucket_FIFO_FIFO(cache,:,:)=N_inter;
    N_max_2Bucket_FIFO_FIFO(cache,:,:)=N_max;
    
    R1_hit_count_Uni(:,cache)=Router1_hit_count;
    R2_hit_count_Uni(:,cache)=Router2_hit_count;
end
toc

% delete(h1);
% clear('h1');
display('Done Uniform');

requests_Uni=zeros(Producers,1);
for ii=1:count
    requests_Uni(producersRequest_Uni(ii))=requests_Uni(producersRequest_Uni(ii))+1;
end

% requests_Uni(1:3,1);
% R1_hit_count_Uni(1:3,1)
% temp1(:,1)=R1_hit_count(1:3,1);
% temp2(:,1)=requests(1:3,1);
% display('Hit rate Simulation for  2 Bucket Uniform Distribution')
% temp1./temp2
temp2=repmat(requests_Uni,1,length(CacheSize));
hit_rate_Simul_Uni=R1_hit_count_Uni./temp2;
% clear temp1 temp2

hit_rate_total_Sim_Uni=sum(R1_hit_count_Uni)/count
%% Theorotical 2 Bucket Model
for cache=1:length(CacheSize)
    for jj=1:Producers
        memory1=CacheSize(cache);
        Prob(:,1)=ProducersProbability_Uni([1:jj-1 jj+1:Producers]);
        f = @(t_c)parameterfunFIFO(t_c,Prob,memory1);
        x=fzero(f,1);
        Characteristic_time_Uni(jj,cache)=x;
%         clear Prob
    end
end
%%
% hit_rate_theor_Uni=zeros(Producers,length(CacheSize));
% for ii=1:3
%     hit_rate_theor_Uni(ii,:)=ones(1,Pop_producers)*(ProducersProbability_Uni(ii)*FreshnessMin)/(1+(ProducersProbability_Uni(ii)*FreshnessMin));
% end
% clear temp1 temp2;
temp1=repmat(ProducersProbability_Uni',1,length(CacheSize));
hit_rate_theor_Uni=(temp1.*Characteristic_time_Uni)./(ones(Producers,length(CacheSize))+(temp1.*Characteristic_time_Uni));

clear temp1;
clear N_min N_inter N_max N_min_temp N_inter_temp N_max_temp jj kk t_inst cache produ i ii
clear ProbForSavingR1 ProbForSavingR2

temp1=repmat(ProducersProbability_Uni',1,length(CacheSize));
hit_rate_total_theor_Uni=sum(hit_rate_theor_Uni.*temp1)
tic
%% Zipf distribution with parameter beta
beta=0.5;
nn=1:Producers;
ProducersProbability_Zipf(1,:)=(nn.^-beta)/sum((nn.^-beta));
producersRequest_Zipf=zeros(count,1);
r=rand(1,count);
for i=1:count
    producersRequest_Zipf(i,1) = sum(r(1,i) >= cumsum([0, ProducersProbability_Zipf]));
end
clear r;

for cache=1:length(CacheSize)
    N_min=zeros(length(ProbForSavingVectorR1),length(ProbForSavingVectorR2));
    N_inter=zeros(length(ProbForSavingVectorR1),length(ProbForSavingVectorR2));
    N_max=zeros(length(ProbForSavingVectorR1),length(ProbForSavingVectorR2));
    for kk=1:length(ProbForSavingVectorR2)
        ProbForSavingR2=ProbForSavingVectorR2(kk);
        for jj=1:length(ProbForSavingVectorR1)
            ProbForSavingR1=ProbForSavingVectorR1(jj);
% memoryR1_FIFO and memoryR2_FIFO have following structure.
% First Column: latest time_instant when data was being used under
% condition it was fresh.
% Second Column: Producer number
% Third Column: time_stamp at ehich data for corresponding producer was
% being fetched and stored.
            memoryR1_FIFO=zeros(CacheSize(cache),2);
            memoryR2_FIFO=zeros(CacheSize(cache),2);

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
                [N_min_temp,N_inter_temp,N_max_temp]=router1_FIFO_plain(produ,t_inst,...
                                                            ProbForSavingR1,ProbForSavingR2,Pop_producers,...
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
        
    N_min_Zipf_FIFO_FIFO(cache,:,:)=N_min;
    N_inter_Zipf_FIFO_FIFO(cache,:,:)=N_inter;
    N_max_Zipf_FIFO_FIFO(cache,:,:)=N_max;
    
    R1_hit_count_Zipf(:,cache)=Router1_hit_count;
    R2_hit_count_Zipf(:,cache)=Router2_hit_count;
end
toc
display('Done Zipf');

requests_Zipf=zeros(Producers,1);
for ii=1:count
    requests_Zipf(producersRequest_Zipf(ii))=requests_Zipf(producersRequest_Zipf(ii))+1;
end

% requests_Uni(1:3,1)
% R1_hit_count_Uni(1:3,1)
% temp1(:,1)=R1_hit_count(1:3,1);
% temp2(:,1)=requests(1:3,1);
% display('Hit rate Simulation for Zipf Distribution')
% temp1./temp2
% R1_hit_count_Zipf./requests_Zipf
temp2=repmat(requests_Zipf,1,length(CacheSize));
hit_rate_Simul_Zipf=R1_hit_count_Zipf./temp2;
% clear temp1 temp2

hit_rate_total_Sim_Zipf=sum(R1_hit_count_Zipf)/count
%% Theorotical calculations
for cache=1:length(CacheSize)
    ProducersProbability_Zipf_sele(1,:)=ProducersProbability_Zipf(1:Producers);
    for jj=1:Producers
        memory1=CacheSize(cache);        
        Prob(:,1)=ProducersProbability_Zipf_sele([1:jj-1 jj+1:Producers]);
        f = @(t_c)parameterfunFIFO(t_c,Prob,memory1);
        x=fzero(f,1);
        Characteristic_time_Zipf(jj,cache)=x;
        clear Prob
    end
end    
temp1=repmat(ProducersProbability_Zipf',1,length(CacheSize));
hit_rate_theor_Zipf=(temp1.*Characteristic_time_Zipf)./(ones(Producers,length(CacheSize))+(temp1.*Characteristic_time_Zipf));

clear temp1;
clear memoryR1_FIFO memoryR2_FIFO N_min N_inter N_max N_min_temp N_inter_temp N_max_temp jj kk t_inst cache produ i ii
clear ProbForSavingR1 ProbForSavingR2

temp1=repmat(ProducersProbability_Zipf',1,length(CacheSize));
hit_rate_total_theor_Zipf=sum(hit_rate_theor_Zipf.*temp1)
%% always change the dataname for saving. Keep it simple and discriptive.
temp1=cd;
cd('D:\IoT\IoT\31Jan\FIFO\Data')
save('Che_Approximation_check__FIFO');
