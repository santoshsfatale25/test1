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
count=10^5;
Producers=[50 75 100 125 150]; % Number of Producers
global Pop_producers

global Freshness_requirment
const=20;
% F_a=5;
% F_b=const*F_a;
% F_c=F_a;
% Freshness_requirment=[F_a F_b F_c];

global Router1_hit_count

ProbForSavingVectorR1=1;%0.2:0.2:1.0;%1.0;
CacheSize=20;%10:5:40;

Prob_a=0.4;%0.25:0.05:0.45;
beta=0.8;%0.5:0.3:1.7;

%% Least Expected Variables :::::::::::::::::::::::::::::::::::::::::::::::
R1_hit_count_Uni_LeastExpe=zeros(max(Producers),length(Producers));

R1_hit_count_Zipf_LeastExpe=zeros(max(Producers),length(Producers));

N_min_3Bucket_LeastExpe=zeros(length(CacheSize),length(Producers));
N_max_3Bucket_LeastExpe=zeros(length(CacheSize),length(Producers));

N_min_Zipf_LeastExpe=zeros(length(CacheSize),length(Producers));
N_max_Zipf_LeastExpe=zeros(length(CacheSize),length(Producers));

%% LRU Variables ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
R1_hit_count_Uni_LRU=zeros(max(Producers),length(Producers));

R1_hit_count_Zipf_LRU=zeros(max(Producers),length(Producers));

N_min_3Bucket_LRU=zeros(length(CacheSize),length(Producers));
N_max_3Bucket_LRU=zeros(length(CacheSize),length(Producers));

N_min_Zipf_LRU=zeros(length(CacheSize),length(Producers));
N_max_Zipf_LRU=zeros(length(CacheSize),length(Producers));

%% LFU Variables ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
R1_hit_count_Uni_LFU=zeros(max(Producers),length(Producers));

R1_hit_count_Zipf_LFU=zeros(max(Producers),length(Producers));

N_min_3Bucket_LFU=zeros(length(CacheSize),length(Producers));
N_max_3Bucket_LFU=zeros(length(CacheSize),length(Producers));

N_min_Zipf_LFU=zeros(length(CacheSize),length(Producers));
N_max_Zipf_LFU=zeros(length(CacheSize),length(Producers));

%% RAND Variables ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
R1_hit_count_Uni_RAND=zeros(max(Producers),length(Producers));

R1_hit_count_Zipf_RAND=zeros(max(Producers),length(Producers));

N_min_3Bucket_RAND=zeros(length(CacheSize),length(Producers));
N_max_3Bucket_RAND=zeros(length(CacheSize),length(Producers));

N_min_Zipf_RAND=zeros(length(CacheSize),length(Producers));
N_max_Zipf_RAND=zeros(length(CacheSize),length(Producers));
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

% Characteristic_time_Uni=zeros(Producers,length(CacheSize));
% Characteristic_time_Zipf=zeros(Producers,length(CacheSize));

global memoryR1_LeastExpe memoryR1_LRU memoryR1_LFU memoryR1_RAND Probability_producers
Probability_producers=zeros(1,max(Producers));

global count1 count2 % Checks cache is empty or not.

% Exponential inter-arrival time
time=cumsum(exprnd(1,count,1));
%% ######################################### 3 Bucket Uniform Distribution ################################################

tic;
Prob_b=0.99*Prob_a;
Prob_c=ones(1,length(Prob_a))-Prob_a-Prob_b;
ProducersProbability_Uni=zeros(max(Producers),length(Producers));
% Freshness_Uni=zeros(max(Producers),length(Producers));
for pp=1:length(Producers)
    N_a=CacheSize;
    N_b=CacheSize;
    N_c=Producers(pp)-N_a-N_b;
    ProducersProbability_Uni(1:Producers(pp),pp)=[repmat(Prob_a./N_a,N_a,1);repmat(Prob_b./N_b,N_b,1);repmat(Prob_c./N_c,N_c,1)];
%     Freshness_Uni(1:Producers(pp),pp)=[repmat(F_a,N_a,1);repmat(F_b,N_b,1);repmat(F_c,N_a,1)];
end

Freshness_Uni=zeros(max(Producers),length(Producers));
for ii=1:length(Producers)
    Freshness_Uni(1:Producers(ii),ii)=sqrt(const*ProducersProbability_Uni(1,ii)./ProducersProbability_Uni(1:Producers(ii),ii));
end

Freshness_const=Freshness_Uni.*ProducersProbability_Uni;

producersRequest_Uni=zeros(count,length(Producers));

for pp=1:length(Producers)
    producersRequest_Uni(:,pp) =datasample(1:Producers(pp),count,'Weights',ProducersProbability_Uni(1:Producers(pp),pp));
end
clear r;



%% 3 Bucket Uniform Distribution Least Expected

for nn=1:length(Producers)
    N_min=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    N_max=zeros(length(ProbForSavingVectorR1),length(CacheSize));
%     global Probability_producers
    for cache=1:length(CacheSize)
        Probability_producers(1,1:Producers(nn))=ProducersProbability_Uni(1:Producers(nn),nn);
        N_a=CacheSize;
        N_b=CacheSize;
        N_c=Producers(nn)-N_a-N_b;
        Pop_producers=[N_a N_b N_c];
        Freshness_requirment=Freshness_Uni(:,nn);
        
        for jj=1:length(ProbForSavingVectorR1)
            ProbForSavingR1=ProbForSavingVectorR1(jj);
% memoryR1_LeastExpe and memoryR2_LeastExpe have following structure.
% First Column: Producer number
% Second Column: t_inst at which it is being fetched from producer

            memoryR1_LeastExpe=zeros(CacheSize(cache),2);

            Router1_hit_count=zeros(Producers(nn),1);

            count1=0;
            count2=0;

            message=sprintf('Running for Cache Size=%d and ProbForSavingR1=%f Probability_a=%f'...
                            ,CacheSize(cache),ProbForSavingR1,Prob_a);
            h=msgbox(message);
            clear message

            N_min_temp=0;
            N_max_temp=0;
            display('3 Bucket Uniform Distribution');
            for ii=1:length(time)
%                 display(t_inst);
                produ=producersRequest_Uni(ii,nn);
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
%     clear global Probability_producers;
    end
    
    R1_hit_count_Uni_LeastExpe(1:Producers(nn),nn)=Router1_hit_count;
    N_min_3Bucket_LeastExpe(:,nn)=N_min;
    N_max_3Bucket_LeastExpe(:,nn)=N_max;
    
end
toc

% delete(h1);
% clear('h1');
display('Done Uniform');


requests_Uni=zeros(max(Producers),length(Producers));
for pp=1:length(Producers)
    for ii=1:count
        requests_Uni(producersRequest_Uni(ii,pp),pp)=requests_Uni(producersRequest_Uni(ii,pp),pp)+1;
    end
end



% requests_Uni(1:3,1);
% R1_hit_count_Uni(1:3,1)
% temp1(:,1)=R1_hit_count(1:3,1);
% temp2(:,1)=requests(1:3,1);
% display('Hit rate Simulation for  2 Bucket Uniform Distribution')
% temp1./temp2
% temp2=repmat(requests_Uni,1,length(Producers));
hit_rate_Simul_Uni_LeastExpe=R1_hit_count_Uni_LeastExpe./requests_Uni;
clear temp1 temp2
hit_rate_total_Sim_Uni_LeastExpe=sum(R1_hit_count_Uni_LeastExpe)/count;


%% 3 Bucket Uniform Distribution LFU (Least Frequently Used)
Probability_producers=zeros(1,max(Producers));
Freshness_requirment=0;

for nn=1:length(Producers)
    N_min=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    N_max=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    for cache=1:length(CacheSize)
        Probability_producers(1,1:Producers(nn))=ProducersProbability_Uni(1:Producers(nn),nn);
        N_a=CacheSize;
        N_b=CacheSize;
        N_c=Producers(nn)-N_a-N_b;
        Pop_producers=[N_a N_b N_c];
        Freshness_requirment=Freshness_Uni(:,nn);

        for jj=1:length(ProbForSavingVectorR1)
            ProbForSavingR1=ProbForSavingVectorR1(jj);
% memoryR1_LeastExpe and memoryR2_LeastExpe have following structure.
% First Column: Producer number
% Second Column: t_inst at which it is being fetched from producer

            memoryR1_LFU=zeros(CacheSize(cache),2);

            Router1_hit_count=zeros(Producers(nn),1);

            count1=0;
            count2=0;

            message=sprintf('Running for Cache Size=%d and ProbForSavingR1=%f Probability_a=%f'...
                            ,CacheSize(cache),ProbForSavingR1,Prob_a);
            h=msgbox(message);
            clear message

            N_min_temp=0;
            N_max_temp=0;
            display('3 Bucket Uniform Distribution');
            for ii=1:length(time)
%                 display(t_inst);
                produ=producersRequest_Uni(ii,nn);
                t_inst=time(ii);
%                 Frshness=Freshness(t_inst,1);
                [N_min_temp,N_max_temp]=router1_LFU_3Bucket(produ,t_inst,...
                                                            ProbForSavingR1,N_min_temp,N_max_temp);

            end
            N_min(jj,cache)=N_min_temp;
            N_max(jj,cache)=N_max_temp;

            delete(h);
            clear('h');
        end        

    end
    R1_hit_count_Uni_LFU(1:Producers(nn),nn)=Router1_hit_count;
    N_min_3Bucket_LFU(:,nn)=N_min;
    N_max_3Bucket_LFU(:,nn)=N_max;
end
toc

% delete(h1);
% clear('h1');
display('Done Uniform');

% requests_Uni(1:3,1);
% R1_hit_count_Uni(1:3,1)
% temp1(:,1)=R1_hit_count(1:3,1);
% temp2(:,1)=requests(1:3,1);
% display('Hit rate Simulation for  2 Bucket Uniform Distribution')
% temp1./temp2
% temp2=repmat(requests_Uni,1,length(Producers));
hit_rate_Simul_Uni_LFU=R1_hit_count_Uni_LFU./requests_Uni;
clear temp1 temp2
hit_rate_total_Sim_Uni_LFU=sum(R1_hit_count_Uni_LFU)/count;


%% 3 Bucket Uniform Distribution RAND (Random)
Probability_producers=zeros(1,max(Producers));
Freshness_requirment=0;

for nn=1:length(Producers)
    N_min=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    N_max=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    for cache=1:length(CacheSize)
        Probability_producers(1,1:Producers(nn))=ProducersProbability_Uni(1:Producers(nn),nn);
        N_a=CacheSize(cache);
        N_b=CacheSize(cache);
        N_c=Producers(nn)-N_a-N_b;
        Pop_producers=[N_a N_b N_c];
        Freshness_requirment=Freshness_Uni(:,nn);
        
        for jj=1:length(ProbForSavingVectorR1)
            ProbForSavingR1=ProbForSavingVectorR1(jj);
% memoryR1_LeastExpe and memoryR2_LeastExpe have following structure.
% First Column: Producer number
% Second Column: t_inst at which it is being fetched from producer

            memoryR1_RAND=zeros(CacheSize(cache),2);

            Router1_hit_count=zeros(Producers(nn),1);

            count1=0;
            count2=0;

            message=sprintf('Running for Cache Size=%d and ProbForSavingR1=%f Probability_a=%f'...
                            ,CacheSize(cache),ProbForSavingR1,Prob_a);
            h=msgbox(message);
            clear message

            N_min_temp=0;
            N_max_temp=0;
            display('3 Bucket Uniform Distribution');
            for ii=1:length(time)
%                 display(t_inst);
                produ=producersRequest_Uni(ii,nn);
                t_inst=time(ii);
%                 Frshness=Freshness(t_inst,1);
                [N_min_temp,N_max_temp]=router1_RAND_plain_3Bucket(produ,t_inst,...
                                                            ProbForSavingR1,N_min_temp,N_max_temp);

            end
            N_min(jj,cache)=N_min_temp;
            N_max(jj,cache)=N_max_temp;

            delete(h);
            clear('h');
        end        

    end
    R1_hit_count_Uni_RAND(1:Producers(nn),nn)=Router1_hit_count;
    N_min_3Bucket_RAND(:,nn)=N_min;
    N_max_3Bucket_RAND(:,nn)=N_max;
end
toc

% delete(h1);
% clear('h1');
display('Done Uniform');

% requests_Uni(1:3,1);
% R1_hit_count_Uni(1:3,1)
% temp1(:,1)=R1_hit_count(1:3,1);
% temp2(:,1)=requests(1:3,1);
% display('Hit rate Simulation for  2 Bucket Uniform Distribution')
% temp1./temp2
% temp2=repmat(requests_Uni,1,length(Producers));
hit_rate_Simul_Uni_RAND=R1_hit_count_Uni_RAND./requests_Uni;
clear temp1 temp2
hit_rate_total_Sim_Uni_RAND=sum(R1_hit_count_Uni_RAND)/count;


%% 3 Bucket Uniform Distribution LRU (Least Recently Used)
Probability_producers=zeros(1,max(Producers));
Freshness_requirment=0;

for nn=1:length(Producers)
    N_min=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    N_max=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    for cache=1:length(CacheSize)
        Probability_producers(1,1:Producers(nn))=ProducersProbability_Uni(1:Producers(nn),nn);
        N_a=CacheSize(cache);
        N_b=CacheSize(cache);
        N_c=Producers(nn)-N_a-N_b;
        Pop_producers=[N_a N_b N_c];
        Freshness_requirment=Freshness_Uni(:,nn);
        
        for jj=1:length(ProbForSavingVectorR1)
            ProbForSavingR1=ProbForSavingVectorR1(jj);
% memoryR1_LRU have following structure.
% First Column: t_stamp
% Second Column: Producer number
% Third Column: t_inst at which it is being fetched from producer
            memoryR1_LRU=zeros(CacheSize(cache),3);

            Router1_hit_count=zeros(Producers(nn),1);

            message=sprintf('Running for Cache Size=%d and ProbForSavingR1=%f and Probability_a=%f'...
                            ,CacheSize(cache),ProbForSavingR1,Prob_a);
            h=msgbox(message);
            clear message

            N_min_temp=0;
            N_max_temp=0;
            display('3 Bucket Uniform Distribution');
            for ii=1:length(time)
%                 display(t_inst);
                produ=producersRequest_Uni(ii,nn);
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

    end
    R1_hit_count_Uni_LRU(1:Producers(nn),nn)=Router1_hit_count;
    N_min_3Bucket_LRU(:,nn)=N_min;
    N_max_3Bucket_LRU(:,nn)=N_max;
end
toc


display('Done Uniform');

% temp2=repmat(requests_Uni,1,length(Producers));
hit_rate_Simul_Uni_LRU=R1_hit_count_Uni_LRU./requests_Uni;
clear temp1 temp2
hit_rate_total_Sim_Uni_LRU=sum(R1_hit_count_Uni_LRU)/count;


%% ###################################### Zipf Distribution with parameter beta #######################################
clear nn
ProducersProbability_Zipf=zeros(length(Producers),max(Producers));
% Freshness_Zipf=zeros(max(Producers),length(Producers));
for pp=1:length(Producers)
    nn=1:Producers(pp);
    ProducersProbability_Zipf(pp,1:Producers(pp))=(nn.^-beta)/sum((nn.^-beta));
%     Freshness_Zipf(1:Producers(pp),pp)=Freshness_const(1:Producers(pp),pp)./(ProducersProbability_Zipf(pp,1:Producers(pp))');
end

Freshness_Zipf=zeros(max(Producers),length(Producers));
for ii=1:length(Producers)
    Freshness_Zipf(1:Producers(ii),ii)=sqrt(const*ProducersProbability_Zipf(ii,1)./ProducersProbability_Zipf(ii,1:Producers(ii)));
end

% Freshness_Zipf=Freshness_const./ProducersProbability_Zipf';
producersRequest_Zipf=zeros(count,length(Producers));

clear temp1;
for pp=1:length(Producers)
    temp1(:,1)=ProducersProbability_Zipf(pp,1:Producers(pp));
    producersRequest_Zipf(:,pp)=datasample(1:Producers(pp),count,'Weights',temp1);
    clear temp1
end
clear r;

%% Zipf distribution with parameter beta Least Expected
Probability_producers=zeros(1,max(Producers));
Freshness_requirment=0;

for nn=1:length(Producers)
    N_min=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    N_max=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    Probability_producers(1,1:Producers(nn))=ProducersProbability_Zipf(nn,1:Producers(nn));
    for cache=1:length(CacheSize)
        N_a=CacheSize(cache);
        N_b=CacheSize(cache);
        N_c=Producers(nn)-N_a-N_b;
        Pop_producers=[N_a N_b N_c];
        Freshness_requirment=Freshness_Zipf(:,nn);

        for jj=1:length(ProbForSavingVectorR1)
            ProbForSavingR1=ProbForSavingVectorR1(jj);
    % memoryR1_LeastExpe and memoryR2_LeastExpe have following structure.
    % First Column: latest time_instant when data was being used under
    % condition it was fresh.
    % Second Column: Producer number
    % Third Column: time_stamp at ehich data for corresponding producer was
    % being fetched and stored.
            memoryR1_LeastExpe=zeros(CacheSize(cache),2);

            Router1_hit_count=zeros(Producers(nn),1);

            count1=0;
            count2=0;

            message=sprintf('Running for Cache Size=%d and ProbForSavingR1=%f and beta=%f'...
                            ,CacheSize(cache),ProbForSavingR1,beta);
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
                                                            ProbForSavingR1,N_min_temp,N_max_temp);

            end
            N_min(jj,cache)=N_min_temp; % jj-> row number; kk-> column number
            N_max(jj,cache)=N_max_temp;

            delete(h);
            clear('h');
        end

    end
    R1_hit_count_Zipf_LeastExpe(1:Producers(nn),nn)=Router1_hit_count;
    N_min_Zipf_LeastExpe(:,nn)=N_min;
    N_max_Zipf_LeastExpe(:,nn)=N_max;
end
toc
display('Done Zipf');

requests_Zipf=zeros(max(Producers),length(Producers));
for pp=1:length(Producers)
    for ii=1:count
        requests_Zipf(producersRequest_Zipf(ii,pp),pp)=requests_Zipf(producersRequest_Zipf(ii,pp),pp)+1;
    end
end


% requests_Uni(1:3,1)
% R1_hit_count_Uni(1:3,1)
% temp1(:,1)=R1_hit_count(1:3,1);
% temp2(:,1)=requests(1:3,1);
% display('Hit rate Simulation for Zipf Distribution')
% temp1./temp2
% R1_hit_count_Zipf./requests_Zipf
% temp2=repmat(requests_Zipf,1,1,length(Producers));
hit_rate_Simul_Zipf_LeastExpe=R1_hit_count_Zipf_LeastExpe./requests_Zipf;
clear temp1 temp2
hit_rate_total_Sim_Zipf_LeastExpe=sum(R1_hit_count_Zipf_LeastExpe)/count;

%% Zipf distribution with parameter beta LFU (Least Frequently Used)
Probability_producers=zeros(1,max(Producers));
Freshness_requirment=0;

for nn=1:length(Producers)
    N_min=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    N_max=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    Probability_producers(1,1:Producers(nn))=ProducersProbability_Zipf(nn,1:Producers(nn));

    for cache=1:length(CacheSize)
        N_a=CacheSize(cache);
        N_b=CacheSize(cache);
        N_c=Producers(nn)-N_a-N_b;
        Pop_producers=[N_a N_b N_c];
        Freshness_requirment=Freshness_Zipf(:,nn);

        for jj=1:length(ProbForSavingVectorR1)
            ProbForSavingR1=ProbForSavingVectorR1(jj);
    % memoryR1_LeastExpe and memoryR2_LeastExpe have following structure.
    % First Column: latest time_instant when data was being used under
    % condition it was fresh.
    % Second Column: Producer number
    % Third Column: time_stamp at ehich data for corresponding producer was
    % being fetched and stored.
            memoryR1_LFU=zeros(CacheSize(cache),2);

            Router1_hit_count=zeros(Producers(nn),1);

            count1=0;
            count2=0;

            message=sprintf('Running for Cache Size=%d and ProbForSavingR1=%f and beta=%f'...
                            ,CacheSize(cache),ProbForSavingR1,beta);
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
                [N_min_temp,N_max_temp]=router1_LFU_3Bucket(produ,t_inst,...
                                                            ProbForSavingR1,N_min_temp,N_max_temp);

            end
            N_min(jj,cache)=N_min_temp; % jj-> row number; kk-> column number
            N_max(jj,cache)=N_max_temp;

            delete(h);
            clear('h');
        end

    end
    R1_hit_count_Zipf_LFU(1:Producers(nn),nn)=Router1_hit_count;    
    N_min_Zipf_LFU(:,nn)=N_min;
    N_max_Zipf_LFU(:,nn)=N_max;
end
toc
display('Done Zipf');


% requests_Uni(1:3,1)
% R1_hit_count_Uni(1:3,1)
% temp1(:,1)=R1_hit_count(1:3,1);
% temp2(:,1)=requests(1:3,1);
% display('Hit rate Simulation for Zipf Distribution')
% temp1./temp2
% R1_hit_count_Zipf./requests_Zipf
% temp2=repmat(requests_Zipf,1,1,length(Producers));
hit_rate_Simul_Zipf_LFU=R1_hit_count_Zipf_LFU./requests_Zipf;
clear temp1 temp2
hit_rate_total_Sim_Zipf_LFU=sum(R1_hit_count_Zipf_LFU)/count;

%% Zipf distribution with parameter beta RAND (RANDOM)
Probability_producers=zeros(1,max(Producers));
Freshness_requirment=0;

for nn=1:length(Producers)
    N_min=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    N_max=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    Probability_producers(1,1:Producers(nn))=ProducersProbability_Zipf(nn,1:Producers(nn));

    for cache=1:length(CacheSize)
        N_a=CacheSize(cache);
        N_b=CacheSize(cache);
        N_c=Producers(nn)-N_a-N_b;
        Pop_producers=[N_a N_b N_c];
        Freshness_requirment=Freshness_Zipf(:,nn);
        
        for jj=1:length(ProbForSavingVectorR1)
            ProbForSavingR1=ProbForSavingVectorR1(jj);
    % memoryR1_LeastExpe and memoryR2_LeastExpe have following structure.
    % First Column: latest time_instant when data was being used under
    % condition it was fresh.
    % Second Column: Producer number
    % Third Column: time_stamp at ehich data for corresponding producer was
    % being fetched and stored.
            memoryR1_RAND=zeros(CacheSize(cache),2);

            Router1_hit_count=zeros(Producers(nn),1);

            count1=0;
            count2=0;

            message=sprintf('Running for Cache Size=%d and ProbForSavingR1=%f and beta=%f'...
                            ,CacheSize(cache),ProbForSavingR1,beta);
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
                [N_min_temp,N_max_temp]=router1_RAND_plain_3Bucket(produ,t_inst,...
                                                            ProbForSavingR1,N_min_temp,N_max_temp);

            end
            N_min(jj,cache)=N_min_temp; % jj-> row number; kk-> column number
            N_max(jj,cache)=N_max_temp;

            delete(h);
            clear('h');
        end

    end
    R1_hit_count_Zipf_RAND(1:Producers(nn),nn)=Router1_hit_count;
    N_min_Zipf_RAND(:,nn)=N_min;
    N_max_Zipf_RAND(:,nn)=N_max;
end
toc
display('Done Zipf');


% requests_Uni(1:3,1)
% R1_hit_count_Uni(1:3,1)
% temp1(:,1)=R1_hit_count(1:3,1);
% temp2(:,1)=requests(1:3,1);
% display('Hit rate Simulation for Zipf Distribution')
% temp1./temp2
% R1_hit_count_Zipf./requests_Zipf
% temp2=repmat(requests_Zipf,1,1,length(Producers));
hit_rate_Simul_Zipf_RAND=R1_hit_count_Zipf_RAND./requests_Zipf;
clear temp1 temp2
hit_rate_total_Sim_Zipf_RAND=sum(R1_hit_count_Zipf_RAND)/count;

%% Zipf distribution with parameter beta LRU (Least Recently Used)
Probability_producers=zeros(1,max(Producers));
Freshness_requirment=0;

for nn=1:length(Producers)
    N_min=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    N_max=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    Probability_producers(1,1:Producers(nn))=ProducersProbability_Zipf(nn,1:Producers(nn));

    for cache=1:length(CacheSize)
        N_a=CacheSize(cache);
        N_b=CacheSize(cache);
        N_c=Producers(nn)-N_a-N_b;
        Pop_producers=[N_a N_b N_c];
        Freshness_requirment=Freshness_Zipf(:,nn);
        
        for jj=1:length(ProbForSavingVectorR1)
            ProbForSavingR1=ProbForSavingVectorR1(jj);
% memoryR1_LRU and memoryR2_LRU have following structure.
% First Column: latest time_instant when data was being used under
% condition it was fresh.
% Second Column: Producer number
% Third Column: time_stamp at which data for corresponding producer was
% being fetched and stored.
            memoryR1_LRU=zeros(CacheSize(cache),3);

            Router1_hit_count=zeros(Producers(nn),1);

            message=sprintf('Running for Cache Size=%d and ProbForSavingR1=%f and beta=%f'...
                            ,CacheSize(cache),ProbForSavingR1,beta);
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
                [N_min_temp,N_max_temp]=router1_LRU_plain_3Bucket(produ,t_inst,...
                                                            ProbForSavingR1,N_min_temp,N_max_temp);%,Sele_1,Sele_2);

            end
            N_min(jj,cache)=N_min_temp;
            N_max(jj,cache)=N_max_temp;

            delete(h);
            clear('h');
        end
        
    end
    R1_hit_count_Zipf_LRU(1:Producers(nn),nn)=Router1_hit_count;
    N_min_Zipf_LRU(:,nn)=N_min;
    N_max_Zipf_LRU(:,nn)=N_max;
end
toc
display('Done Zipf');

% temp2=repmat(requests_Zipf,1,1,length(Producers));
hit_rate_Simul_Zipf_LRU=R1_hit_count_Zipf_LRU./requests_Zipf;
clear temp1 temp2
hit_rate_total_Sim_Zipf_LRU=sum(R1_hit_count_Zipf_LRU)/count;

%% Therotical Upper Bound 3 Bucket Distribution
% N_a=CacheSize;
% N_b=CacheSize;
% upperBound1_Uni=zeros(1,length(Producers));
% for ii=1:length(Producers)
%     N_c=Producers(ii)-2*CacheSize;
% %     Freshness=[repmat(F_a,N_a,1);repmat(F_b,N_b,1);repmat(F_c,N_c,1)];
%     temp1(:,1)=ProducersProbability_Uni(1:Producers(ii),ii);
%     upperBound1_Uni(1,ii)=sum(((temp1.^2).*Freshness)./(ones(Producers(ii),1)+temp1.*Freshness));
%     clear temp1;
% end
upperBound1_Uni(1,:)=sum((Freshness_Uni.*(ProducersProbability_Uni.^2))./(1+Freshness_Uni.*ProducersProbability_Uni));

upperBound2_Uni=zeros(1,length(Producers));
for nn=1:length(Producers)
    upperBound2_Uni(1,nn)=sum(ProducersProbability_Uni(1:CacheSize,nn));
end


upperBoundMin_Uni=min(upperBound1_Uni,upperBound2_Uni);
%% Therotical Upper Bound Zipf Distribution
% N_a=CacheSize;
% N_b=CacheSize;
% upperBound1_Zipf=zeros(1,length(Producers));
% for ii=1:length(Producers)
%     N_c=Producers(ii)-2*CacheSize;
%     Freshness=[repmat(F_a,N_a,1);repmat(F_b,N_b,1);repmat(F_c,N_c,1)];
%     temp1(:,1)=ProducersProbability_Zipf(ii,1:Producers(ii));
%     upperBound1_Zipf(1,ii)=sum(((temp1.^2).*Freshness)./(ones(Producers(ii),1)+temp1.*Freshness));
%     clear temp1;
% end
upperBound1_Zipf(1,:)=sum((Freshness_Zipf.*(ProducersProbability_Zipf.^2)')./1+(Freshness_Zipf.*ProducersProbability_Zipf'));

upperBound2_Zipf=zeros(1,length(Producers));
for nn=1:length(Producers)
    upperBound2_Zipf(1,nn)=sum(ProducersProbability_Zipf(nn,1:CacheSize));
end


upperBoundMin_Zipf=min(upperBound1_Zipf,upperBound2_Zipf);

%% Result Plot
% myplot(xinput,yinputMatrix,xlabel1,ylabel1,title1,legend1,saveFigAs)
clear temp1;
temp1=cd;
xinput(:,1)=Producers;
yinputMatrix=horzcat(upperBoundMin_Uni',hit_rate_total_Sim_Uni_LeastExpe',hit_rate_total_Sim_Uni_LRU',hit_rate_total_Sim_Uni_RAND',hit_rate_total_Sim_Uni_LFU');
xlabel1=sprintf('Producers');
ylabel1=sprintf('Hit rate');
% title1=sprintf('Hit rate (p_{hit}) Vs Producers');
directory='D:\IoT\IoT\31Jan\LeastExpected\New_Results24092017';
legend1={sprintf('UB'),sprintf('LU'),sprintf('LRU'),sprintf('RAND'),sprintf('LFU')};
saveFigAs=sprintf('Hit_rate_Vs_Producers_policies_Uniform');
myplot(xinput,yinputMatrix,xlabel1,ylabel1,legend1,saveFigAs,directory);
cd(temp1);

yinputMatrix=horzcat(upperBoundMin_Zipf',hit_rate_total_Sim_Zipf_LeastExpe',hit_rate_total_Sim_Zipf_LRU',hit_rate_total_Sim_Zipf_RAND',hit_rate_total_Sim_Zipf_LFU');
xlabel1=sprintf('Producers');
ylabel1=sprintf('Hit rate');
% title1=sprintf('Hit rate (p_{hit}) Vs Producers');
directory='D:\IoT\IoT\31Jan\LeastExpected\New_Results24092017';
legend1={sprintf('UB'),sprintf('LU'),sprintf('LRU'),sprintf('RAND'),sprintf('LFU')};
saveFigAs=sprintf('Hit_rate_Vs_Producers_policies_Zipf');
myplot(xinput,yinputMatrix,xlabel1,ylabel1,legend1,saveFigAs,directory);

%% always change the dataname for saving. Keep it simple and discriptive.
% temp1=cd;
cd('D:\IoT\IoT\31Jan\LeastExpected\New_Results24092017\Data')
save('cmp_producers_policies_3Class_Sqrt_Freshness'); % F_L is freshenss requirement for Bucket 1 and Bucket 3. 
