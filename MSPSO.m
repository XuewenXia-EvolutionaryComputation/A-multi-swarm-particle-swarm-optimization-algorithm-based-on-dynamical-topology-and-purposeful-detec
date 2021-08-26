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

shreshold=[12,7,4,3,2,1,1]; % ͣ������Ⱥ��������ֵ���뵱ǰ����Ⱥ������Ӧ
select_no = 1;
group_ps=PN(select_no);
group_num=Particle_Number/group_ps;
ps=group_num*group_ps;

recorded = 0;  % �ﵽ����ʱ��¼�����Ϣ
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

VRmin = Xmin;   % ��2��Ϊ�Լ�����ԭ����ӵ�
VRmax = Xmax;

%% �������ռ���еȷִ�������Ϣ��ȡ��صı�����ʼ�� -------------------------------------------------------
Rn = 10;
Len(1:D) = (VRmax - VRmin)/Rn; % ÿ������ĳ��ȡ�Len(i)��ʾ��iά�����ȷֺ�ÿ������ĳ��ȡ�
Region = zeros(Rn+1,D);  % Region �ĵ�j�б�ʾ��jά�������������䡣����i��j�б�ʾ��jά�����ĵ�i�����䡣
Region(1,:) = VRmin;     % ��ΪRn�����䣬��Ӧ���� Rn+1 ���ָ�㣨������β����RegionΪ Rn+1��D�С� 
for i = 2:Rn+1
   Region(i,:) = VRmin + (i-1).*Len; 
end

Bestsection =zeros(Rn,D);  % ��¼pbest�ڵĵ�jά�����ڵ�i�������ϳ��ֵĴ���
Goodsection = zeros(Rn,D); % ��i��j�е�ֵ��ʾ�������ӵĵ�jά�����ڵ�i�������ϳ��ֵĴ�����
Bedsection = zeros(Rn,D);  % ��i��j�е�ֵ��ʾ�������ӵĵ�jά�����ڵ�i�������ϳ��ֵĴ�����

% %   ÿ������Ⱥ gbest �ֱ����̽��ʱ��Ҫ�õ��ı���
% DetectMaxsection =zeros(group_num,Rn,D);  % ��¼ÿ������Ⱥ��gbest�ڵ�jά�����ڵ�i���������Ƿ�̽�����0:δ̽�⣻ 1����̽��
% DetectMinsection =zeros(group_num,Rn,D); 
% Staysection=zeros(group_num,1,D);  % ��¼ÿ������Ⱥ��gbest�ڶ�ά��ͣ����ĳһά�Ĵ���
% gbestsection =zeros(group_num,Rn,D); % ��¼ÿ������Ⱥ��gbest�ڵĵ�jά�����ڵ�i�������ϳ��ֵĴ���
% gbesttected = zeros(group_num,D);    % ��¼ÿ������Ⱥ��gbest�ڵĵ�jά�������Ƿ�̽���

%   ������Ⱥ gbest �ֱ����̽��ʱ��Ҫ�õ��ı���
DetectMaxsection =zeros(Rn,D);  % ��¼��Ⱥ��gbest�ڵ�jά�����ڵ�i���������Ƿ�̽�����0:δ̽�⣻ 1����̽��
DetectMinsection =zeros(Rn,D); 
Staysection=zeros(1,D);  % ��¼��Ⱥ��gbest�ڶ�ά��ͣ����ĳһά�Ĵ���
gbestsection =zeros(Rn,D); % ��¼��Ⱥ��gbest�ڵĵ�jά�����ڵ�i�������ϳ��ֵĴ���
gbesttected = zeros(1,D);    % ��¼��Ⱥ��gbest�ڵĵ�jά�������Ƿ�̽���

Attractor = Goodsection - Bedsection;
Regbest = zeros(10,D);  % ��¼���10����Ⱥ����ֵ��ÿά�ϳ��ֵĴ���

% ��¼��ά�����İ���Ϣ��0��ʾΪ�󶨣��������ģ���0��ʾ�Ѱ󶨡�
% ��Binding[1]=Binding[3]=1��Binding[4]=Binding[5]=2�ֱ��ʾ��1,3ά����һ�𡢵�4,5ά����һ��
Binding = zeros(1,D); 
maxDensity = ps/2;     % ��Bestsection���ж�ά������ͬһ�����ϳ��ֵ��ܶȴ���maxDensityʱ����������������а󶨡�

Memory = 3;  % ������������Memory������
BestFly = zeros(Memory,Rn,D); % ��¼Memory��gbest���µķ��й켣����section�ķ�ʽ��¼��

maxSave = 4;   % ÿ������Ⱥ�ڵ��������б�������ֵ��������
ms = 1;
SaveBestPos1 = ones(maxSave,D)*inf; % ����Ⱥ1����������λ�ã���ʼ��
SaveBestVal1 = ones(maxSave,1)*inf;
SaveBestVel1 = ones(maxSave,D)*inf;
count1 = 0; % ��¼�ѱ�������Ÿ�����Ŀ
SaveBestPos2 = ones(maxSave,D)*inf; % ����Ⱥ2����������λ�ã���ʼ��
SaveBestVal2 = ones(maxSave,1)*inf;
SaveBestVel2 = ones(maxSave,D)*inf;
count2 = 0; % ��¼�ѱ�������Ÿ�����Ŀ

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
stop_group = zeros(1,group_num); % ��¼����Ⱥ���Ž�ͣ�ʹ�������ʼ��Ϊ0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
group_ps=PN(select_no);
group_num=Particle_Number/group_ps;
neighbor=[];
for ii=1:group_num
    % % Ring topological��  LPSO
    neighbor(group_id(ii,1),:)=[group_id(ii,group_ps),group_id(ii,2)];
    for jj=2:(group_ps-1)
        neighbor(group_id(ii,jj),:)=[group_id(ii,jj-1),group_id(ii,jj+1)];
    end
    neighbor(group_id(ii,group_ps),:)=[group_id(ii,group_ps-1),group_id(ii,1)];

    old_gbest(ii,:) = gbest(ii,:);
    old_gbestval(ii) = gbestval(ii);
end
stop_group = zeros(1,group_num); % ��¼����Ⱥ���Ž�ͣ�ʹ�������ʼ��Ϊ0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%% ����Ϊ��ͼ��
% old = 1;
% new = fitcount;
% yyy(old:new) = min(gbestval);
% old = new;

get_flag=0;
i=0;
group_num_changed = 0;  

% ��¼����������������Ⱥ�� gbest ͣ�ʹ���
globalgbeststop = 0;
 % ��¼����ǰ��������Ⱥ�� gbest �����Ӧ����Ӧֵ��
[C I] = min(pbestval);   
oldglobalgbest = pbest(I,:);
oldglobalgbestv = C; 
globalgbest = pbest(I,:);
globalgbestv = C;
globalgbest_changed = 0; 

while fitcount<Max_FES%0.95*Max_FES    
    i=i+1;
    
    % ��¼����ǰ��������Ⱥ�� gbest �����Ӧ����Ӧֵ��
    [C I] = min(pbestval);   
    oldglobalgbest = pbest(I,:);
    oldglobalgbestv = C;  
    globalgbest_changed = 0;  
    
    avgfit_group = [];  % ÿ������Ⱥ��ǰ��ƽ����Ӧֵ
    avgdis_group = [];  % ÿ������Ⱥÿά�������� 0 ��ľ���
    
%     if group_num_changed ==1  % ������Ⱥ���������ı�
%         gbestsection =zeros(group_num,Rn,D);  % ��¼ÿ������Ⱥ��gbest�ڵĵ�jά�����ڵ�i�������ϳ��ֵĴ���
%     end
    
    for k=1:ps
        if 0%mod(i,2)==0
            [tmp,tmpid]=min(pbestval(neighbor(group_id(ceil(k/group_ps),:)))); 
            aa(k,:)=cc(1).*rand(1,D).*(pbest(k,:)-pos(k,:))+cc(2).*rand(1,D).*(pbest(neighbor(group_id(ceil(k/group_ps),tmpid)),:)-pos(k,:));
        else
            aa(k,:)=cc(1).*rand(1,D).*(pbest(k,:)-pos(k,:))+cc(2).*rand(1,D).*(gbest(pos_group(k),:)-pos(k,:)); 
        end
        vel(k,:)=iwt(i)*vel(k,:)+aa(k,:); 
        vel(k,:)=(vel(k,:)>(alpha(i)*mv(k,:))).*(alpha(i)*mv(k,:))+(vel(k,:)<=alpha(i)*mv(k,:)).*vel(k,:);  % Ĭ�ϵķ������� line 121 ���ʹ��
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
        [G_id C]=find(group_id==k);  % G_id Ϊ�� k ����������������Ⱥ���
        if new_group ~= G_id         % ͳ��ĳ������Ⱥ����ʷ���Ž��Ƿ�ͣ��
            new_group = G_id;
            stop_group(pos_group(k))=stop_group(pos_group(k))+1;  % ÿ��Ĭ������1��������� if ��������������Ϊ 0 ��
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
     
%     %%%% ����Ϊ��ͼ��
%     new = fitcount;
%     yyy(old:new) = newglobalgbestv;
%     old = new;     
    
 %% ����̽�����
    detected = 0;
    if mod(i,group_ps)==0  %mod(i,Cycle)==0 %  
      %%  ���ɡ���Ŀ�ĵ�̽��
        detected = 1;
        % ���ݵ�ǰ��Ⱥ���ӵ�pbest����Bestsection
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
       %%  ����������Ⱥ�� gbest����̽��---------------- 1------------
        [C I] = min(pbestval);
        globalgbest = pbest(I,:);
        globalgbestv = C;  
         % ���ݵ�ǰ��Ⱥ���ӵ�gbest����gbestsection
        for j = 1:D
           for k = 1:size(Region,1)-1
               if globalgbest(1,j) >= Region(k,j) && globalgbest(1,j) <= Region(k+1,j)
                   gbestsection(k,j) =  gbestsection(k,j) + 1;  
                   break;
               end
           end
        end
        
        gbest_imp_flag=zeros(1,D);  % ��¼ gbest ����̽����Ƿ�õ����Ƶı�־        
       
        for j=1:D

            [gval garea] = max(gbestsection(:,j));  % �ҵ�gbest��jά��������
            [pvalmax pareamax] = sort(Bestsection(:,j),'descend');  % �ҵ�pbest��jά�����������򡣼�¼��Bestsection��jά��
            [pvalmin pareamin] = sort(Bestsection(:,j),'ascend');  % �ҵ�pbest��jά�������ٵ����򡣼�¼��Bestsection��jά��
            tmpbest=globalgbest;
            % Ŀǰ�����⣺������ܼ������ϡ�������ֵ���ʱ����gbest̽������ѡ��ǰ������Ӧ��ʹ��ѭ��̽����ƣ��������������
            if garea==pareamax(1)  % ���gbest��jά����pbest�ڸ�ά�������ܼ�������������pbest�ڸ�ά�������ٵ�����̽��
                for jj=Rn:-1:1
                    if DetectMinsection(pareamin(jj),j) == 0  % �ҵ�δ̽�⵽��pbest�������ٵ����� 0��δ̽�⣻ 1����̽��
                        break;
                    end                
                end

                tmpbest(1,j)=Region(pareamin(jj),j)+rand*Len(j);
                DetectMinsection(pareamin(jj),j)=1;  % ����Ϊ����̽��
%             else  % gbest��jά�Ȳ��ڶ��ܼ�������Ҳ�������ٵ�����,�����ѡ����Bestsection�����ܼ�����������̽��
%                 if rand<=0.5  
%                     for jj=Rn:-1:1
%                         if DetectMaxsection(pareamax(jj),j) == 0  % �ҵ�δ̽�⵽��pbest������������ 0��δ̽�⣻ 1����̽��
%                             break;
%                         end                
%                     end
% 
%                     tmpbest(1,j)=Region(pareamax(jj),j)+rand*Len(j);
%                     DetectMaxsection(pareamax(jj),j)=1;  % ����Ϊ����̽��
%                 else
%                     for jj=Rn:-1:1
%                         if DetectMinsection(pareamin(jj),j) == 0  % �ҵ�δ̽�⵽��pbest�������ٵ����� 0��δ̽�⣻ 1����̽��
%                             break;
%                         end                
%                     end
% 
%                     tmpbest(1,j)=Region(pareamin(jj),j)+rand*Len(j);
%                     DetectMinsection(pareamin(jj),j)=1;  % ����Ϊ����̽��
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
               gbest_imp_flag(j)=1;  % ��̽���õ����ƣ����ø�άΪ1
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
            if sum(DetectMinsection(:,jj))==Rn  % ������е�����̽���
                DetectMinsection(:,jj)=0;       % ���临λΪȫ0
            end
            if sum(DetectMaxsection(:,jj))==Rn  % ������е�����̽���
                DetectMaxsection(:,jj)=0;       % ���临λΪȫ0
            end
        end        
%        
       %%  �Ը�����Ⱥ�е� gbest����̽��---------------- 2------------
%         for ii=1:group_num
%            gbest_imp_flag=zeros(ii,D);  % ��¼ gbest ����̽����Ƿ�õ����Ƶı�־,��ʼ����Ϊ 0 ��
%         end
%        
%         for ii=1:group_num
%             
%             tb =  gbest(ii,:);  % ��¼�� ii ������Ⱥ�� gbest 
%             for j=1:D
%                 for k = 1:size(Region,1)-1
%                    if gbest(pos_group(ii),j) >= Region(k,j) && gbest(pos_group(ii),j) <= Region(k+1,j)
%                        gbestsection(ii,k,j) =  gbestsection(ii,k,j) + 1;  
%                        break;
%                    end
%                 end                
%                 
%                 [gval garea] = max(gbestsection(ii,:,j));  % �ҵ���ii������Ⱥ��gbest��jά��������
%                 [pvalmax pareamax] = sort(Bestsection(:,j),'descend');  % �ҵ�pbest��jά�����������򡣼�¼��Bestsection��jά��
%                 [pvalmin pareamin] = sort(Bestsection(:,j),'ascend');  % �ҵ�pbest��jά�������ٵ����򡣼�¼��Bestsection��jά��
%                 tmpbest=gbest(ii,:);   % ̽��ǰ������Ⱥgbest���ݴ�
%                 
%                 % Ŀǰ�����⣺������ܼ������ϡ�������ֵ���ʱ����gbest̽������ѡ��ǰ������Ӧ��ʹ��ѭ��̽����ƣ��������������
%                 % �������Ⱥgbest��jά����pbest�ڸ�ά�������ܼ�������,�Ҹ�����Ⱥδ�ڸ�ά��ִ��̽�������������pbest�ڸ�ά�������ٵ�����̽��
%                 if garea==pareamax(1) && gbesttected(ii,j)==0 
%                     for jj=Rn:-1:1
%                         if DetectMinsection(ii,pareamin(jj),j) == 0  % �ҵ�δ̽�⵽��pbest�������ٵ����� 0��δ̽�⣻ 1����̽��
%                             break;
%                         end                
%                     end
% 
%                     tmpbest(1,j)=Region(pareamin(jj),j)+rand*Len(j);
%                     DetectMinsection(ii,pareamin(jj),j)=1;  % ����Ϊ����̽��
%                 elseif garea==pareamin(1) % % �������Ⱥgbest��jά����pbest�ڸ�ά�������ٵ�����������pbest�ڸ�ά�������ܼ�������̽��
% %                     for jj=Rn:-1:1
% %                         if DetectMaxsection(pareamax(jj),j) == 0  % �ҵ�δ̽�⵽��pbest������������ 0��δ̽�⣻ 1����̽��
% %                             break;
% %                         end                
% %                     end
% %                     tmpbest(1,j)=Region(pareamax(jj),j)+rand*Len(j);
% %                     DetectMaxsection(pareamax(jj),j)=1;  % ����Ϊ����̽��
%                 else  % ����Ⱥgbest��jά�Ȳ��ڶ��ܼ�������Ҳ�������ٵ�����,�����ѡ����Bestsection�����ܼ�����������̽��
%                     if rand<=0.5  
%                         for jj=Rn:-1:1
%                             if DetectMaxsection(ii,pareamax(jj),j) == 0  % �ҵ�δ̽�⵽��pbest������������ 0��δ̽�⣻ 1����̽��
%                                 break;
%                             end                
%                         end
% 
%                         tmpbest(1,j)=Region(pareamax(jj),j)+rand*Len(j);
%                         DetectMaxsection(ii,pareamax(jj),j)=1;  % ����Ϊ����̽��
%                     else
%                         for jj=Rn:-1:1
%                             if DetectMinsection(ii,pareamin(jj),j) == 0  % �ҵ�δ̽�⵽��pbest�������ٵ����� 0��δ̽�⣻ 1����̽��
%                                 break;
%                             end                
%                         end
% 
%                         tmpbest(1,j)=Region(pareamin(jj),j)+rand*Len(j);
%                         DetectMinsection(ii,pareamin(jj),j)=1;  % ����Ϊ����̽��
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
%                    gbest_imp_flag(ii,j)=1;  % ��̽���õ����ƣ����ø�άΪ1
%                 end
%             end
%         end
%         
%     %% 
%         for ii=1:group_num
%             for jj=1:D
%                 if sum(DetectMinsection(ii,:,jj))==Rn  % ������е�����̽���
%                     DetectMinsection(ii,:,jj)=0;       % ���临λΪȫ0
%                 end
%                 if sum(DetectMaxsection(ii,:,jj))==Rn  % ������е�����̽���
%                     DetectMaxsection(ii,:,jj)=0;       % ���临λΪȫ0
%                 end
%             end
%         end
%         
%         % ���ѡ��һ������gbest����ĳһά.
%         gbest_changed = 0;  % ��¼�ֲ�����ǰ�� gbestval �Ƿ�仯�ı�־
%         for ii=1:group_num
%             bb=randperm(ps);
%             while pbestval(bb(1))==gbestval(ii)
%                 bb=randperm(ps);
%             end
%             tmppos=pos(bb(1),:);
%             for jj=1:D
%                if  gbest_imp_flag(ii,jj)==0   % �� gbest ֮ǰ��̽��δ��ʹ��õ�����
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
%                if  gbest_imp_flag(ii,jj)==0   % �� gbest ֮ǰ��̽��δ��ʹ��õ�����
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
       % ���ѡ��һ������gbest����ĳһά. 
        gbest_changed = 0;  % ��¼�ֲ�����ǰ�� gbestval �Ƿ�仯�ı�־
        bb=randperm(ps);
        if pbestval(bb(1))==globalgbestv
            bb=randperm(ps);
        end
        tmppos=pos(bb(1),:);
        for ii=1:D
           if  gbest_imp_flag(ii)==0   % �� gbest ֮ǰ��̽��δ��ʹ��õ�����
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
           if  gbest_imp_flag(ii)==0   % �� gbest ֮ǰ��̽��δ��ʹ��õ�����
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
    % mod(i,50���滻Ϊ mod(i,(floor(0.1*Max_FES/ps)))==0
    % �����ֲ����������뵽����Ⱥ�����ģ��С�Ĵ�����
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
%     %%%% ����Ϊ��ͼ��
%     new = fitcount;
%     yyy(old:new) = newglobalgbestv;
%     old = new;
%     
    % ������������Ⱥ�� gbest δ�õ��Ż�ʱ
    if oldglobalgbestv <= newglobalgbestv
        globalgbeststop = globalgbeststop + 1;
    else
        globalgbeststop = 0;
    end
    % ----------------------- �������飡
    maxstop_best = floor(group_ps/2);  % �˴����ǵ�������Ⱥ�����Ž����Ϣ��Ring�ṹ�µĴ����ٶȡ����ƺ����㷨������Ⱥ�����õ�GPSOģʽ
    if  globalgbeststop >= maxstop_best %||mod(i,R)==0
        rc=randperm(ps);
        group_id=[];gbest=[];gbestval=[];pos_group=zeros(1,ps); neighbor=[];
        for k=1:group_num
            group_id(k,:)=rc(((k-1)*group_ps+1):k*group_ps);
            pos_group(group_id(k,:))=k;
            [gbestval(k),gbestid]=min(pbestval(group_id(k,:)));
            gbest(k,:)=pbest(group_id(k,gbestid),:);
            
             % % Ring topological��  LPSO
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
    
    % ----------------------- ���кϲ�������������Ⱥ��ģ��
    stop_ratio = 1/2; % ����Ⱥͣ�͵ı���������2/3������Ⱥ����ͣ��
     if select_no<length(PN) &&(mod(i,(floor(1/(length(PN)+10)*Max_FES/ps)))==0 ) %( (mod(i,(floor(0.10*Max_FES/ps)))==0 )|| sum(stop_group>=50)/length(stop_group)>=stop_ratio )  % �����ģ��С�����ڵ�����δ����������Ⱥ�����ﵽ������Ⱥ������75%����ʱ
%    if select_no<length(PN) && (sum(stop_group>=group_ps)/length(stop_group))>=stop_ratio  % �����ģ��С�����ڵ�����δ����������Ⱥ�����ﵽ������Ⱥ������stop_ratio����ʱ
%    if select_no<length(PN) && fitcount >= Max_FES*(select_no/length(PN)) && fitcount <= Max_FES*((select_no+1)/length(PN))
%--------  �ɿ��ǽ���ţ�ٷ��������������������Ⱥ��Ҫ�����ģ��Сʱ���ȶԴ���Ⱥ����ʷ���Ž��оֲ�������Ȼ��������Ⱥ��ģ��
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

             % % Ring topological��  LPSO
            neighbor(group_id(ii,1),:)=[group_id(ii,group_ps),group_id(ii,2)];
            for jj=2:(group_ps-1)
                neighbor(group_id(ii,jj),:)=[group_id(ii,jj-1),group_id(ii,jj+1)];
            end
            neighbor(group_id(ii,group_ps),:)=[group_id(ii,group_ps-1),group_id(ii,1)];
    
            old_gbest(ii,:) = gbest(ii,:);
            old_gbestval(ii) = gbestval(ii);
        end
        stop_group = zeros(1,group_num); % ��¼����Ⱥ���Ž�ͣ�ʹ�������ʼ��Ϊ0
        
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
%     %--------- ����ΪCEC2005-Dymanic Multi-Swarm Particle Swarm Optimizer ��˼��
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
%    %%%% ����Ϊ��ͼ��
%     new = fitcount;
%     yyy(old:new) = min(gbestval);
%     old = new;
%     
% % if gbestval<err0&get_flag==0
% %     fit_cut=fitcount;get_flag=1;
% % end
% 
%    %%%%%%%%%%%%%%%  �������趨���������ĵ����۴���
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

