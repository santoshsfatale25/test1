function [Value,index2]=FindLeastExpe(producers,t_stamps,t_inst)
global Probability_producers Pop_producers FreshnessMin FreshnessMax
% producers
temp1(:,1)=Probability_producers(producers);
temp2=zeros(length(producers),1);
index1=find(producers<Pop_producers+1);
temp2(index1,1)=t_stamps(index1)+FreshnessMin-t_inst;
index1=find(producers>Pop_producers);
temp2(index1,1)=t_stamps(index1)+FreshnessMax-t_inst;
% temp1.*temp2
[Value,index2]=min(temp1.*temp2);

