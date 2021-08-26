% function [gbest,gbestval,fitcount,fit_cut,get_flag]= DMS_PSO_func(fhd,me,err0,err1,Max_FES,group_num,Particle_Number,Dimension,Xmin,Xmax,VRmin,VRmax,varargin)
function [gbest,gbestval,fitcount,suc,suc_fes]= MSPSO_func(jingdu,func_num,fhd,Dimension,Particle_Number,me,Max_FES,Xmin,Xmax,varargin)
%                                                            (jingdu,func_num,fhd,Dimension,Particle_Number,Max_Gen,Max_FES,VRmin,VRmax,varargin)
global orthm best_f best_keep initial_flag  fbias  gbias
%usage:     [gbest,gbestval,fitcount,fit_cut,get_flag]= DMS_PSO_func(fhd,err0,err1,Max_FES,group_num,Particle_Number,Dimension,Xmin,Xmax,VRmin,VRmax,varargin)
%this program used to solve minimization problem, if you want wot solve a
%maximization problem please add "-" in your function.
%Outputs:
%           gbest:      the final result
%           gbestval:   the corresponding 
%           fitcount:   the corresponding used FEs
%           fit_cut:    the used FEs when achieve gbestval<err0
%           get_flag:   get_flag=1 means get the predefined criterion gbestval<err0
%Inputs:
%           fhd:        your problem function 
%           err0:       criterion, the FEs used for achieve this criterion will be recorded
%           err1:       stopping criterion  (please set err0/err1 to -Inf if you do not want the algorithm stop before getting Max_FEs)
%           Max_FEs:    The number of max function evaluations
%           group_num:  The group number of DMS-PSO
%           Particle_Number:    The number of particles in EACH group
%           Dimension:  the dimension of your problem
%           Xmin,Xmax:  the bounds of initilization
%           VRmin,VRmax:the bounds of search space
%           varargin:   all the parameters you need in your problem function

rand('state',sum(100*clock));

PN = [2,3,5,6,10,15,30];

shreshold=[12,7,4,3,2,1,1]; % 停滞子种群个数的阈值，与当前子种群个数对应
select_no = 1;
group_ps=PN(select_no);
group_num=Particle_Number/group_ps;
ps=group_num*group_ps;

recorded = 0;  % 达到精度时记录相关信息
suc = 0;
suc_fes = 0;

D=Dimension;
L=100;
L_FES=200;
L_num=ceil(0.25.*group_num);
cc=[1.49445 1.49445];   %acceleration constants
iwt=0.729*ones(1,me);
% cc=[2 2];
% iwt=0.9-(1:me)*(0.5/me);
R=10; % Regrouping period
Cycle=3;

VRmin = Xmin;   % 这2句为自己根据原文添加的
VRmax = Xmax;

%% 将变量空间进行等分处理及与信息提取相关的变量初始化 -------------------------------------------------------
Rn = 10;
Len(1:D) = (VRmax - VRmin)/Rn; % 每个区间的长度。Len(i)表示第i维变量等分后每个区间的长度。
Region = zeros(Rn+1,D);  % Region 的第j列表示第j维变量的整个区间。即第i行j列表示第j维变量的第i个区间。
Region(1,:) = VRmin;     % 分为Rn个区间，则应该有 Rn+1 个分割点（包括首尾），Region为 Rn+1行D列。 
for i = 2:Rn+1
   Region(i,:) = VRmin + (i-1).*Len; 
end

Bestsection =zeros(Rn,D);  % 记录pbest在的第j维变量在第i个区间上出现的次数
Goodsection = zeros(Rn,D); % 第i行j列的值表示优良粒子的第j维变量在第i个区间上出现的次数。
Bedsection = zeros(Rn,D);  % 第i行j列的值表示劣质粒子的第j维变量在第i个区间上出现的次数。

% %   每个子种群 gbest 分别进行探测时需要用到的变量
% DetectMaxsection =zeros(group_num,Rn,D);  % 记录每个子种群的gbest在第j维变量在第i个区间上是否探测过：0:未探测； 1：已探测
% DetectMinsection =zeros(group_num,Rn,D); 
% Staysection=zeros(group_num,1,D);  % 记录每个子种群的gbest在多维上停滞在某一维的次数
% gbestsection =zeros(group_num,Rn,D); % 记录每个子种群的gbest在的第j维变量在第i个区间上出现的次数
% gbesttected = zeros(group_num,D);    % 记录每个子种群的gbest在的第j维变量上是否探测过

%   整个种群 gbest 分别进行探测时需要用到的变量
DetectMaxsection =zeros(Rn,D);  % 记录种群的gbest在第j维变量在第i个区间上是否探测过：0:未探测； 1：已探测
DetectMinsection =zeros(Rn,D); 
Staysection=zeros(1,D);  % 记录种群的gbest在多维上停滞在某一维的次数
gbestsection =zeros(Rn,D); % 记录种群的gbest在的第j维变量在第i个区间上出现的次数
gbesttected = zeros(1,D);    % 记录种群的gbest在的第j维变量上是否探测过

Attractor = Goodsection - Bedsection;
Regbest = zeros(10,D);  % 记录最近10代种群最优值在每维上出现的次数

% 记录多维变量的绑定信息。0表示为绑定，即独立的；非0表示已绑定。
% 如Binding[1]=Binding[3]=1，Binding[4]=Binding[5]=2分别表示第1,3维绑定在一起、第4,5维绑定在一起
Binding = zeros(1,D); 
maxDensity = ps/2;     % 当Bestsection中有多维变量在同一区间上出现的密度大于maxDensity时，则对这多个变量进行绑定。

Memory = 3;  % 最优粒子连续Memory代更新
BestFly = zeros(Memory,Rn,D); % 记录Memory代gbest更新的飞行轨迹，用section的方式记录。

maxSave = 4;   % 每个子种群在迭代过程中保留最优值的最大个数
ms = 1;
SaveBestPos1 = ones(maxSave,D)*inf; % 子种群1保留的最优位置：初始化
SaveBestVal1 = ones(maxSave,1)*inf;
SaveBestVel1 = ones(maxSave,D)*inf;
count1 = 0; % 记录已保存的最优个体数目
SaveBestPos2 = ones(maxSave,D)*inf; % 子种群2保留的最优位置：初始化
SaveBestVal2 = ones(maxSave,1)*inf;
SaveBestVel2 = ones(maxSave,D)*inf;
count2 = 0; % 记录已保存的最优个体数目

abc = 0;
 
if length(VRmin)==1
    VRmin=repmat(VRmin,1,D);
    VRmax=repmat(VRmax,1,D);
end
%%
if length(VRmin)==1
    VRmin=repmat(VRmin,1,D);
    VRmax=repmat(VRmax,1,D);
end

VRmin=repmat(VRmin,ps,1);
VRmax=repmat(VRmax,ps,1);

if length(Xmin)==1
    Xmin=repmat(Xmin,1,D);
    Xmax=repmat(Xmax,1,D);
end
Xmin=repmat(Xmin,ps,1);
Xmax=repmat(Xmax,ps,1);

mv=0.20*(Xmax-Xmin);
alpha=1*ones(1,me);
alpha2=5-(1:me).*(4.9./me);


Vmin=-0.1*(Xmax-Xmin);
Vmax=-Vmin;

pos=Xmin+(Xmax-Xmin).*rand(ps,D);

for i=1:ps;
    e(i,1)=feval(fhd,pos(i,:)',varargin{:})-fbias(func_num);
end

fitcount=ps;
vel=Vmin+2.*Vmax.*rand(ps,D);%initialize the velocity of the particles
pbest=pos;
pbestval=e; %initialize the pbest and the pbest's fitness value
for i=1:group_num
    group_id(i,:)=[((i-1)*group_ps+1):i*group_ps];
    pos_group(group_id(i,:))=i;
    [gbestval(i),gbestid]=min(pbestval(group_id(i,:)));
    gbest(i,:)=pbest(group_id(i,gbestid),:);%initialize the gbest and the gbest's fitness value
   
    old_gbest(i,:) = gbest(i,:);
    old_gbestval(i) = gbestval(i);
end
stop_group = zeros(1,group_num); % 记录子种群最优解停滞代数，初始均为0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
group_ps=PN(select_no);
group_num=Particle_Number/group_ps;
neighbor=[];
for ii=1:group_num
    % % Ring topological：  LPSO
    neighbor(group_id(ii,1),:)=[group_id(ii,group_ps),group_id(ii,2)];
    for jj=2:(group_ps-1)
        neighbor(group_id(ii,jj),:)=[group_id(ii,jj-1),group_id(ii,jj+1)];
    end
    neighbor(group_id(ii,group_ps),:)=[group_id(ii,group_ps-1),group_id(ii,1)];

    old_gbest(ii,:) = gbest(ii,:);
    old_gbestval(ii) = gbestval(ii);
end
stop_group = zeros(1,group_num); % 记录子种群最优解停滞代数，初始均为0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%% 以下为画图用
% old = 1;
% new = fitcount;
% yyy(old:new) = min(gbestval);
% old = new;

get_flag=0;
i=0;
group_num_changed = 0;  

% 记录进化过程中整个种群的 gbest 停滞代数
globalgbeststop = 0;
 % 记录进化前的整个种群的 gbest 及其对应的适应值。
[C I] = min(pbestval);   
oldglobalgbest = pbest(I,:);
oldglobalgbestv = C; 
globalgbest = pbest(I,:);
globalgbestv = C;
globalgbest_changed = 0; 

while fitcount<Max_FES%0.95*Max_FES    
    i=i+1;
    
    % 记录进化前的整个种群的 gbest 及其对应的适应值。
    [C I] = min(pbestval);   
    oldglobalgbest = pbest(I,:);
    oldglobalgbestv = C;  
    globalgbest_changed = 0;  
    
    avgfit_group = [];  % 每个子种群当前的平均适应值
    avgdis_group = [];  % 每个子种群每维变量距离 0 点的距离
    
%     if group_num_changed ==1  % 若子种群个数发生改变
%         gbestsection =zeros(group_num,Rn,D);  % 记录每个子种群的gbest在的第j维变量在第i个区间上出现的次数
%     end
    
    for k=1:ps
        if 0%mod(i,2)==0
            [tmp,tmpid]=min(pbestval(neighbor(group_id(ceil(k/group_ps),:)))); 
            aa(k,:)=cc(1).*rand(1,D).*(pbest(k,:)-pos(k,:))+cc(2).*rand(1,D).*(pbest(neighbor(group_id(ceil(k/group_ps),tmpid)),:)-pos(k,:));
        else
            aa(k,:)=cc(1).*rand(1,D).*(pbest(k,:)-pos(k,:))+cc(2).*rand(1,D).*(gbest(pos_group(k),:)-pos(k,:)); 
        end
        vel(k,:)=iwt(i)*vel(k,:)+aa(k,:); 
        vel(k,:)=(vel(k,:)>(alpha(i)*mv(k,:))).*(alpha(i)*mv(k,:))+(vel(k,:)<=alpha(i)*mv(k,:)).*vel(k,:);  % 默认的方法，与 line 121 配合使用
        vel(k,:)=(vel(k,:)<(-alpha(i)*mv(k,:))).*(-alpha(i)*mv(k,:))+(vel(k,:)>=(-alpha(i)*mv(k,:))).*(vel(k,:));    

        pos(k,:)=pos(k,:)+vel(k,:); 
       
        if rand>0.5
            pos(k,:)=(pos(k,:)>VRmax(1,:)).*VRmax(1,:)+(pos(k,:)<=VRmax(1,:)).*pos(k,:); 
            pos(k,:)=(pos(k,:)<VRmin(1,:)).*VRmin(1,:)+(pos(k,:)>=VRmin(1,:)).*pos(k,:);
        else
            pos(k,:)=((pos(k,:)>=VRmin(1,:))&(pos(k,:)<=VRmax(1,:))).*pos(k,:)...
                        +(pos(k,:)<VRmin(1,:)).*(VRmin(1,:)+0.10.*(VRmax(1,:)-VRmin(1,:)).*rand(1,D))...
                        +(pos(k,:)>VRmax(1,:)).*(VRmax(1,:)-0.10.*(VRmax(1,:)-VRmin(1,:)).*rand(1,D));
        end
        
        e(k,1)=feval(fhd,pos(k,:)',varargin{:})-fbias(func_num);
        fitcount=fitcount+1;           
    end
     
    new_group = 0;
    for k=1:ps        
        [G_id C]=find(group_id==k);  % G_id 为第 k 个粒子所处的子种群编号
        if new_group ~= G_id         % 统计某个子种群的历史最优解是否停滞
            new_group = G_id;
            stop_group(pos_group(k))=stop_group(pos_group(k))+1;  % 每次默认自增1，若下面的 if 语句成立，则再置为 0 ！
        end
        
        tmp=(pbestval(k)<e(k));
        temp=repmat(tmp,1,D);
        pbest(k,:)=temp.*pbest(k,:)+(1-temp).*pos(k,:);
        pbestval(k)=tmp.*pbestval(k)+(1-tmp).*e(k);%update the pbest
        if pbestval(k)<gbestval(pos_group(k))
            gbest(pos_group(k),:)=pbest(k,:);
            gbestval(pos_group(k))=pbestval(k);
            stop_group(pos_group(k))=0;
        end
    end    

    [CC II] = min(pbestval);
    newgb = pbest(II,:);
    newgbv = CC;
    newglobalgbestv = min(gbestval);
    
    %%%%%%%%%%%%%%%
     if CC <= jingdu && recorded == 0
         recorded = 1;
         gbestval(1)=CC;
         gbest(1,:)=pbest(II,:);
         suc = 1;
         suc_fes = fitcount; 
%          break;
     end
     %%%%%%%%%%%%%%%
     
%     %%%% 以下为画图用
%     new = fitcount;
%     yyy(old:new) = newglobalgbestv;
%     old = new;     
    
 %% 禁忌探测机制
    detected = 0;
    if mod(i,group_ps)==0  %mod(i,Cycle)==0 %  
      %%  禁忌、有目的地探测
        detected = 1;
        % 根据当前种群粒子的pbest更新Bestsection
        for ii=1:ps        
            for j = 1:D
               for k = 1:size(Region,1)-1
                   if pbest(ii,j) >= Region(k,j) && pbest(ii,j) <= Region(k+1,j)
                       Bestsection(k,j) =  Bestsection(k,j) + 1;  
                       break;
                   end
               end
            end
        end
       %%  对整个大种群的 gbest进行探测---------------- 1------------
        [C I] = min(pbestval);
        globalgbest = pbest(I,:);
        globalgbestv = C;  
         % 根据当前种群粒子的gbest更新gbestsection
        for j = 1:D
           for k = 1:size(Region,1)-1
               if globalgbest(1,j) >= Region(k,j) && globalgbest(1,j) <= Region(k+1,j)
                   gbestsection(k,j) =  gbestsection(k,j) + 1;  
                   break;
               end
           end
        end
        
        gbest_imp_flag=zeros(1,D);  % 记录 gbest 进行探测后是否得到改善的标志        
       
        for j=1:D

            [gval garea] = max(gbestsection(:,j));  % 找到gbest第j维所在区间
            [pvalmax pareamax] = sort(Bestsection(:,j),'descend');  % 找到pbest第j维出现最多的区域。记录在Bestsection第j维中
            [pvalmin pareamin] = sort(Bestsection(:,j),'ascend');  % 找到pbest第j维出现最少的区域。记录在Bestsection第j维中
            tmpbest=globalgbest;
            % 目前的问题：当多个密集区域或稀疏区域的值相等时，则gbest探测总是选择靠前的区域，应该使用循环探测机制，避免这种情况。
            if garea==pareamax(1)  % 如果gbest第j维处于pbest在该维出现最密集的区域，则向最pbest在该维出现最少的区域探测
                for jj=Rn:-1:1
                    if DetectMinsection(pareamin(jj),j) == 0  % 找到未探测到的pbest出现最少的区域。 0：未探测； 1：已探测
                        break;
                    end                
                end

                tmpbest(1,j)=Region(pareamin(jj),j)+rand*Len(j);
                DetectMinsection(pareamin(jj),j)=1;  % 设置为：已探测
%             else  % gbest第j维既不在对密集的区域也不在最少的区域,则随机选择向Bestsection中最密集或最少区域探测
%                 if rand<=0.5  
%                     for jj=Rn:-1:1
%                         if DetectMaxsection(pareamax(jj),j) == 0  % 找到未探测到的pbest出现最多的区域。 0：未探测； 1：已探测
%                             break;
%                         end                
%                     end
% 
%                     tmpbest(1,j)=Region(pareamax(jj),j)+rand*Len(j);
%                     DetectMaxsection(pareamax(jj),j)=1;  % 设置为：已探测
%                 else
%                     for jj=Rn:-1:1
%                         if DetectMinsection(pareamin(jj),j) == 0  % 找到未探测到的pbest出现最少的区域。 0：未探测； 1：已探测
%                             break;
%                         end                
%                     end
% 
%                     tmpbest(1,j)=Region(pareamin(jj),j)+rand*Len(j);
%                     DetectMinsection(pareamin(jj),j)=1;  % 设置为：已探测
%                 end
            end
            
            Staysection(1,j)=Staysection(1,j)+1;
            tmpbest=(tmpbest>VRmax(1,:)).*VRmax(1,:)+(tmpbest<=VRmax(1,:)).*tmpbest; 
            tmpbest=(tmpbest<VRmin(1,:)).*VRmin(1,:)+(tmpbest>=VRmin(1,:)).*tmpbest;
            fitness=feval(fhd,tmpbest',varargin{:})-fbias(func_num);            
            fitcount=fitcount+1;
           
            if fitness<globalgbestv
               globalgbest=tmpbest;
               globalgbestv=fitness;
               Staysection(1,j)=0;
               gbest_imp_flag(j)=1;  % 若探测后得到改善，则置该维为1
               %%%%%%%%%%%%%%%
                 if fitness <= jingdu && recorded == 0
                     recorded = 1;
                     gbestval(1)=fitness;
                     gbest(1,:)=pbest(II,:);
                     suc = 1;
                     suc_fes = fitcount; 
%                      break;
                 end
                 %%%%%%%%%%%%%%%
            end
        end
        pbest(I,:)=globalgbest(1,:);
        pbestval(I,:) = globalgbestv;
        
        for jj=1:D
            if sum(DetectMinsection(:,jj))==Rn  % 如果所有的区域都探测过
                DetectMinsection(:,jj)=0;       % 则将其复位为全0
            end
            if sum(DetectMaxsection(:,jj))==Rn  % 如果所有的区域都探测过
                DetectMaxsection(:,jj)=0;       % 则将其复位为全0
            end
        end        
%        
       %%  对各子种群中的 gbest进行探测---------------- 2------------
%         for ii=1:group_num
%            gbest_imp_flag=zeros(ii,D);  % 记录 gbest 进行探测后是否得到改善的标志,初始均置为 0 。
%         end
%        
%         for ii=1:group_num
%             
%             tb =  gbest(ii,:);  % 记录第 ii 个子种群的 gbest 
%             for j=1:D
%                 for k = 1:size(Region,1)-1
%                    if gbest(pos_group(ii),j) >= Region(k,j) && gbest(pos_group(ii),j) <= Region(k+1,j)
%                        gbestsection(ii,k,j) =  gbestsection(ii,k,j) + 1;  
%                        break;
%                    end
%                 end                
%                 
%                 [gval garea] = max(gbestsection(ii,:,j));  % 找到第ii个子种群的gbest第j维所在区间
%                 [pvalmax pareamax] = sort(Bestsection(:,j),'descend');  % 找到pbest第j维出现最多的区域。记录在Bestsection第j维中
%                 [pvalmin pareamin] = sort(Bestsection(:,j),'ascend');  % 找到pbest第j维出现最少的区域。记录在Bestsection第j维中
%                 tmpbest=gbest(ii,:);   % 探测前的子种群gbest先暂存
%                 
%                 % 目前的问题：当多个密集区域或稀疏区域的值相等时，则gbest探测总是选择靠前的区域，应该使用循环探测机制，避免这种情况。
%                 % 如果子种群gbest第j维处于pbest在该维出现最密集的区域,且该子种群未在该维上执行探测操作，则向最pbest在该维出现最少的区域探测
%                 if garea==pareamax(1) && gbesttected(ii,j)==0 
%                     for jj=Rn:-1:1
%                         if DetectMinsection(ii,pareamin(jj),j) == 0  % 找到未探测到的pbest出现最少的区域。 0：未探测； 1：已探测
%                             break;
%                         end                
%                     end
% 
%                     tmpbest(1,j)=Region(pareamin(jj),j)+rand*Len(j);
%                     DetectMinsection(ii,pareamin(jj),j)=1;  % 设置为：已探测
%                 elseif garea==pareamin(1) % % 如果子种群gbest第j维处于pbest在该维出现最少的区域，则向最pbest在该维出现最密集的区域探测
% %                     for jj=Rn:-1:1
% %                         if DetectMaxsection(pareamax(jj),j) == 0  % 找到未探测到的pbest出现最多的区域。 0：未探测； 1：已探测
% %                             break;
% %                         end                
% %                     end
% %                     tmpbest(1,j)=Region(pareamax(jj),j)+rand*Len(j);
% %                     DetectMaxsection(pareamax(jj),j)=1;  % 设置为：已探测
%                 else  % 子种群gbest第j维既不在对密集的区域也不在最少的区域,则随机选择向Bestsection中最密集或最少区域探测
%                     if rand<=0.5  
%                         for jj=Rn:-1:1
%                             if DetectMaxsection(ii,pareamax(jj),j) == 0  % 找到未探测到的pbest出现最多的区域。 0：未探测； 1：已探测
%                                 break;
%                             end                
%                         end
% 
%                         tmpbest(1,j)=Region(pareamax(jj),j)+rand*Len(j);
%                         DetectMaxsection(ii,pareamax(jj),j)=1;  % 设置为：已探测
%                     else
%                         for jj=Rn:-1:1
%                             if DetectMinsection(ii,pareamin(jj),j) == 0  % 找到未探测到的pbest出现最少的区域。 0：未探测； 1：已探测
%                                 break;
%                             end                
%                         end
% 
%                         tmpbest(1,j)=Region(pareamin(jj),j)+rand*Len(j);
%                         DetectMinsection(ii,pareamin(jj),j)=1;  % 设置为：已探测
%                     end
%                 end
%                 Staysection(ii,1,j)=Staysection(ii,1,j)+1;
%                 fitness=feval(fhd,tmpbest,varargin{:});%-fbias;
%                 if fitness < 0
%                      fitness=fitness; 
%                 end
%                 fitcount=fitcount+1;
% 
%                 if fitness<gbestval(ii)
%                    gbest(ii,:)=tmpbest;
%                    gbestval(ii)=fitness;
%                    Staysection(ii,1,j)=0;
%                    gbest_imp_flag(ii,j)=1;  % 若探测后得到改善，则置该维为1
%                 end
%             end
%         end
%         
%     %% 
%         for ii=1:group_num
%             for jj=1:D
%                 if sum(DetectMinsection(ii,:,jj))==Rn  % 如果所有的区域都探测过
%                     DetectMinsection(ii,:,jj)=0;       % 则将其复位为全0
%                 end
%                 if sum(DetectMaxsection(ii,:,jj))==Rn  % 如果所有的区域都探测过
%                     DetectMaxsection(ii,:,jj)=0;       % 则将其复位为全0
%                 end
%             end
%         end
%         
%         % 随机选择一粒子与gbest交换某一维.
%         gbest_changed = 0;  % 记录局部搜索前后 gbestval 是否变化的标志
%         for ii=1:group_num
%             bb=randperm(ps);
%             while pbestval(bb(1))==gbestval(ii)
%                 bb=randperm(ps);
%             end
%             tmppos=pos(bb(1),:);
%             for jj=1:D
%                if  gbest_imp_flag(ii,jj)==0   % 若 gbest 之前的探测未能使其得到改善
%                    tmpbest=gbest(ii,:);
%                    tmpbest(jj)=tmppos(jj);
%                    fitness=feval(fhd,tmpbest,varargin{:});
%                    fitcount=fitcount+1;
%                    if fitness < 0
%                       fitness=fitness; 
%                    end
%                    if fitness<gbestval(ii)
%                       gbest(ii,:)=tmpbest(1,:);
%                       gbestval(ii)=fitness;
%                       gbest_changed=1;
%                       gbest_imp_flag(ii,jj)=1;
%                    end
%                end
%             end 
%         end
%         
%         for ii=1:group_num
%             bb=randperm(ps);
%             while pbestval(bb(1))==gbestval(ii)
%                 bb=randperm(ps);
%             end
%             tmppos=pos(bb(1),:);
%             for jj=1:D
%                if  gbest_imp_flag(ii,jj)==0   % 若 gbest 之前的探测未能使其得到改善
%                    tmpbest=gbest(ii,:);
%                    tmpbest(jj)=tmppos(jj);
%                    fitness=feval(fhd,tmpbest,varargin{:});
%                    fitcount=fitcount+1;
%                    
%                    if fitness < 0
%                       fitness=fitness; 
%                    end
%                    
%                    if fitness<gbestval(ii)
%                       gbest(ii,:)=tmpbest(1,:);
%                       gbestval(ii)=fitness;
%                       gbest_changed=1;
%                       gbest_imp_flag(ii,jj)=1;
%                    end
%                end
%             end 
%         end
        
    end
        
    if detected == 1
       detected = 0;
       % 随机选择一粒子与gbest交换某一维. 
        gbest_changed = 0;  % 记录局部搜索前后 gbestval 是否变化的标志
        bb=randperm(ps);
        if pbestval(bb(1))==globalgbestv
            bb=randperm(ps);
        end
        tmppos=pos(bb(1),:);
        for ii=1:D
           if  gbest_imp_flag(ii)==0   % 若 gbest 之前的探测未能使其得到改善
                tmpbest=globalgbest(1,:);
%               tmpbest=newgb(1,:);
               tmpbest(ii)=tmppos(ii);              
               fitness=feval(fhd,tmpbest',varargin{:})-fbias(func_num);                
               fitcount=fitcount+1;
               
               if fitness < 0
                  tmpbest=(tmpbest>VRmax(1,:)).*VRmax(1,:)+(tmpbest<=VRmax(1,:)).*tmpbest; 
                  tmpbest=(tmpbest<VRmin(1,:)).*VRmin(1,:)+(tmpbest>=VRmin(1,:)).*tmpbest;
                  fitness=feval(fhd,tmpbest',varargin{:})-fbias(func_num);                  
                  fitcount=fitcount+1; 
               end
               
               if fitness<globalgbestv
                  globalgbest=tmpbest;
                  globalgbestv=fitness;
                  gbest_changed=1;
                  gbest_imp_flag(ii)=1;
                   %%%%%%%%%%%%%%%
                 if fitness <= jingdu && recorded == 0
                     recorded = 1;
                     gbestval(1)=fitness;
                     gbest(1,:)=pbest(II,:);
                     suc = 1;
                     suc_fes = fitcount; 
%                      break;
                 end
                 %%%%%%%%%%%%%%%
               end
           end
        end 
        
        bb=randperm(ps);
        if pbestval(bb(1))==globalgbestv
            bb=randperm(ps);
        end
        tmppos=pos(bb(1),:);
        for ii=1:D
           if  gbest_imp_flag(ii)==0   % 若 gbest 之前的探测未能使其得到改善
               tmpbest=globalgbest(1,:);
%               tmpbest=newgb(1,:);
               tmpbest(ii)=tmppos(ii);
               fitness=feval(fhd,tmpbest',varargin{:})-fbias(func_num);               
               fitcount=fitcount+1;

               if fitness < 0
                  tmpbest=(tmpbest>VRmax(1,:)).*VRmax(1,:)+(tmpbest<=VRmax(1,:)).*tmpbest; 
                  tmpbest=(tmpbest<VRmin(1,:)).*VRmin(1,:)+(tmpbest>=VRmin(1,:)).*tmpbest;
                  fitness=feval(fhd,tmpbest',varargin{:})-fbias(func_num);
                  fitcount=fitcount+1; 
               end
               
               if fitness<globalgbestv
                  globalgbest=tmpbest;
                  globalgbestv=fitness;
                  gbest_changed=1;
                  gbest_imp_flag(ii)=1;
                   %%%%%%%%%%%%%%%
                 if fitness <= jingdu && recorded == 0
                     recorded = 1;
                     gbestval(1)=fitness;
                     gbest(1,:)=pbest(II,:);
                     suc = 1;
                     suc_fes = fitcount; 
%                      break;
                 end
                 %%%%%%%%%%%%%%%
               end
           end
        end
        
        pbest(I,:)=globalgbest(1,:);
        pbestval(I,:) = globalgbestv;
    end
    %%
    % mod(i,50）替换为 mod(i,(floor(0.1*Max_FES/ps)))==0
    % 即：局部搜索可融入到子种群变更规模大小的代码中
%     if mod(i,(floor(0.1*Max_FES/ps)))==0  
%         [C I] = min(pbestval);
%         globalgbest = pbest(I,:);
%         globalgbestv = C;        
%         options = optimset('LargeScale','off','MaxFunEvals',ceil(0.10*fitcount),'Display','off');
%         [x,fval,exitflag,output] = fminunc(fhd,globalgbest(1,:),options,varargin{:});
%         fitcount=fitcount+output.funcCount;        
%         if fval < 0
%            fval=fval; 
%         end               
%         if fval<gbestval
%             globalgbest(1,:)=x;
%             globalgbestv=fval;            
%             pbest(I,:)=globalgbest(1,:);
%             pbestval(I,:) = globalgbestv;
%         end
%     end
    
    newglobalgbestv = min(gbestval);
%     %%%% 以下为画图用
%     new = fitcount;
%     yyy(old:new) = newglobalgbestv;
%     old = new;
%     
    % 当本代整个种群的 gbest 未得到优化时
    if oldglobalgbestv <= newglobalgbestv
        globalgbeststop = globalgbeststop + 1;
    else
        globalgbeststop = 0;
    end
    % ----------------------- 进行重组！
    maxstop_best = floor(group_ps/2);  % 此处考虑到了子种群中最优解的信息在Ring结构下的传播速度。但似乎本算法中子种群均采用的GPSO模式
    if  globalgbeststop >= maxstop_best %||mod(i,R)==0
        rc=randperm(ps);
        group_id=[];gbest=[];gbestval=[];pos_group=zeros(1,ps); neighbor=[];
        for k=1:group_num
            group_id(k,:)=rc(((k-1)*group_ps+1):k*group_ps);
            pos_group(group_id(k,:))=k;
            [gbestval(k),gbestid]=min(pbestval(group_id(k,:)));
            gbest(k,:)=pbest(group_id(k,gbestid),:);
            
             % % Ring topological：  LPSO
            if 1%group_ps >=3
                neighbor(group_id(k,1),:)=[group_id(k,group_ps),group_id(k,2)];
                for jj=2:(group_ps-1)
                    neighbor(group_id(k,jj),:)=[group_id(k,jj-1),group_id(k,jj+1)];
                end
                neighbor(group_id(k,group_ps),:)=[group_id(k,group_ps-1),group_id(k,1)];
            end
        end
        
        globalgbeststop = 0;
    end    
    
    % ----------------------- 进行合并，即扩大子种群规模！
    stop_ratio = 1/2; % 子种群停滞的比例，即有2/3的子种群出现停滞
     if select_no<length(PN) &&(mod(i,(floor(1/(length(PN)+10)*Max_FES/ps)))==0 ) %( (mod(i,(floor(0.10*Max_FES/ps)))==0 )|| sum(stop_group>=50)/length(stop_group)>=stop_ratio )  % 变更规模大小的周期到来或未进化的子种群个数达到总子种群个数的75%以上时
%    if select_no<length(PN) && (sum(stop_group>=group_ps)/length(stop_group))>=stop_ratio  % 变更规模大小的周期到来或未进化的子种群个数达到总子种群个数的stop_ratio以上时
%    if select_no<length(PN) && fitcount >= Max_FES*(select_no/length(PN)) && fitcount <= Max_FES*((select_no+1)/length(PN))
%--------  可考虑将拟牛顿法放在这里。动机：当子种群需要变更规模大小时，先对大种群的历史最优进行局部搜索，然后变更子种群规模。
        [C I] = min(pbestval);
        globalgbest = pbest(I,:);
        globalgbestv = C;        
        options = optimset('LargeScale','off','MaxFunEvals',ceil(0.10*fitcount),'Display','off');
        [x,fval,exitflag,output] = fminunc(fhd,globalgbest(1,:)',options,varargin{:});
        fval=fval-fbias(func_num);
        fitcount=fitcount+output.funcCount;        
        if fval < 0
           fval=fval; 
        end
        if fval<gbestval
            globalgbest(1,:)=x;
            globalgbestv=fval;            
            pbest(I,:)=globalgbest(1,:);
            pbestval(I,:) = globalgbestv;
        end
        
         %%%%%%%%%%%%%%%
         [CC II] = min(pbestval);

         if CC <= jingdu && recorded == 0
             recorded = 1;
             gbestval(1)=CC;
             gbest(1,:)=pbest(II,:);
             suc = 1;
             suc_fes = fitcount; 
%              break;
         end
         %%%%%%%%%%%%%%%
         %%%%%%%%%%%%%%%
     
%-------------------------------        
        select_no = select_no+1;
        group_ps=PN(select_no);
        group_num=Particle_Number/group_ps;
        group_id=[]; pos_group=[]; gbestval=[]; gbest=[]; old_gbest=[]; old_gbestval=[]; neighbor=[];
        for ii=1:group_num
            group_id(ii,:)=[((ii-1)*group_ps+1):ii*group_ps];
            pos_group(group_id(ii,:))=ii;
            [gbestval(ii),gbestid]=min(pbestval(group_id(ii,:)));
            gbest(ii,:)=pbest(group_id(ii,gbestid),:);%initialize the gbest and the gbest's fitness value

             % % Ring topological：  LPSO
            neighbor(group_id(ii,1),:)=[group_id(ii,group_ps),group_id(ii,2)];
            for jj=2:(group_ps-1)
                neighbor(group_id(ii,jj),:)=[group_id(ii,jj-1),group_id(ii,jj+1)];
            end
            neighbor(group_id(ii,group_ps),:)=[group_id(ii,group_ps-1),group_id(ii,1)];
    
            old_gbest(ii,:) = gbest(ii,:);
            old_gbestval(ii) = gbestval(ii);
        end
        stop_group = zeros(1,group_num); % 记录子种群最优解停滞代数，初始均为0
        
    end
    
%     if round(i/100)==i/100
% %     plot(pos(:,D-1),pos(:,D),'b*');
% %     hold on
% %     plot(gbest(:,D-1),gbest(:,D),'r*');   
% %     hold off
% %     axis([VRmin(1,D-1),VRmax(1,D-1),VRmin(1,D),VRmax(1,D)])
% %     title(['PSO: ',num2str(i),' generations, Gbestval=',num2str(min(gbestval))]);  
% %     drawnow
% i,gbestval
%     end

% if min(gbestval)<err0&&get_flag==0
%     fit_cut=fitcount;get_flag=1;
% end

% if min(gbestval)<err1
%     break;
% end
    if (i>=me) && (fitcount<=Max_FES)
        i=i-1;
    end
    %%%%%%%%%%%%%%%
    [C I]=min(gbestval);

     if C <= jingdu && recorded == 0
         recorded = 1;
         gbestval(1)=C;
         gbest(1,:)=pbest(I,:);
         suc = 1;
         suc_fes = fitcount; 
%          break;
     end
     %%%%%%%%%%%%%%%
     %%%%%%%%%%%%%%%
end

[tmp,tmpid]=sort(gbestval);
gbest=gbest(tmpid(1),:);
gbestval=gbestval(tmpid(1));

% % 
% % if gbestval>err1
% options = optimset('LargeScale','off','MaxFunEvals',0.05*Max_FES,'Display','off');
% % options = optimset('LargeScale','off','MaxFunEvals',max_gen,'DiffMinChange',1e-20);
% [x,fval,exitflag,output] = fminunc(fhd,gbest,options,varargin{:});
% fitcount=fitcount+output.funcCount;
% if fval<gbestval
%     gbest=x;
%     gbestval=fval;
% end
% % end

% while fitcount<Max_FES        
%     i=i+1;
%    
%     %--------- 以下为CEC2005-Dymanic Multi-Swarm Particle Swarm Optimizer 的思想
%     options = optimset('LargeScale','off','MaxFunEvals',0.05*Max_FES,'Display','off');
%     [tmp,tmpid]=sort(gbestval);
%     for k=1:L_num
%         [x,fval,exitflag,output] = fminunc(fhd,gbest(1,:),options,varargin{:});
% %         fval=fval-fbias;
%         fitcount=fitcount+output.funcCount;
%         if fval<gbestval
% %             [gbestval(tmpid(k)),gbestid]=min(pbestval(group_id(tmpid(k),:)));
% %             pbest(group_id(tmpid(k),gbestid),:)=x;
% %             pbestval(group_id(tmpid(k),gbestid))=fval;
%             gbest(1,:)=x;
%             gbestval=fval;
%         end
%     end
%     
%     
%    %%%% 以下为画图用
%     new = fitcount;
%     yyy(old:new) = min(gbestval);
%     old = new;
%     
% % if gbestval<err0&get_flag==0
% %     fit_cut=fitcount;get_flag=1;
% % end
% 
%    %%%%%%%%%%%%%%%  测试在设定精度下所耗的评价次数
%     [C I]=min(gbestval);
%      if C <= jingdu
%          gbestval(1)=C;
%          gbest(1,:)=pbest(I,:);
%         break;%return; 
%      end
%      %%%%%%%%%%%%%%%
%      
% %     if fitcount>=Max_FES
% %         break;
% %     end
%     if (i>=me) && (fitcount<=Max_FES)
%         i=i-1;
%     end
% end

% 
% if gbestval<err0 && get_flag==0
%     fit_cut=fitcount;get_flag=1; 
% end
% 
% if get_flag==0
%     fit_cut=fitcount;
% end
% 
yyy(Max_FES)=gbestval;
% 
% x=1:round(Max_FES/10):Max_FES;
% x(11)=Max_FES;
% yyy(Max_FES) = gbestval;
% hold on;
% plot(x,yyy(x),'-ks','MarkerFaceColor','k','MarkerSize',5);

