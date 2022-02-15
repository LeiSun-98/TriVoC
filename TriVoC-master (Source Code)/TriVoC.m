%% TriVoC (Source Code)
% Please refer to paper: TriVoC: Efficient Voting-based Consensus Maximization for Robust Point Cloud Registration with Extreme Outlier Ratios
% Copyright by Lei Sun & Lu Deng, only for academic use

%% Please adjust the key input parameters as follows and directly run the code

n_ele=1000; % correspondence number

noise=0.01; % noise (standard deviation)

outlier_ratio=0.98; % outlier ratio: (0~0.99) in percentage

show_figure=1; % whether to show the correspondences

%% Build Environment
[pts_3d,pts_3d_,R_gt,t_gt]=Build_Scenario_A(n_ele,noise,outlier_ratio,1,show_figure);


%% TriVoc
    
tic;

adj_mat=zeros(n_ele,n_ele);
for i=1:n_ele-1
    for j=i+1:n_ele
           if abs(norm(pts_3d((i),:)-pts_3d((j),:))-norm(pts_3d_((i),:)-pts_3d_((j),:)))<=5*noise
               adj_mat(i,j)=1;adj_mat(j,i)=1;
           end
    end
end

votes=zeros(1,n_ele);
for i=1:n_ele
    votes(i)=sum(adj_mat(i,:));
end
[~,voters]=sort(votes,'descend');

best_set=[];best_size=0;max_num=459;rand_samp=randperm(n_ele);max_num2=45;max_num3=45;

for seq=1:n_ele*0.1

ele_this=voters(seq);

cou=0;potential_set=zeros(1);
for i=1:n_ele
    if adj_mat(i,ele_this)==1 && i~=ele_this
        cou=cou+1;
        potential_set(cou)=i;
    end
end

potential_size=length(potential_set);

if potential_size>best_size-1

adj_mat2=adj_mat(potential_set,potential_set);
    
votes2=zeros(1,potential_size);
for i=1:potential_size
    votes2(i)=sum(adj_mat2(i,:));
end
[~,voters2]=sort(votes2,'descend');
voters2=potential_set(voters2);

for seq2=1:potential_size*0.1
    
ele_this2=voters2(seq2);
    
cou=0;potential_set2=zeros(1);
for i=1:potential_size
    if adj_mat(potential_set(i),ele_this2)==1 && potential_set(i)~=ele_this && potential_set(i)~=ele_this2
        cou=cou+1;
        potential_set2(cou)=potential_set(i);
    end
end
potential_size2=length(potential_set2);

if potential_size2>best_size-2 && potential_set2(1)>0
    
adj_mat3=adj_mat(potential_set2,potential_set2);

votes3=zeros(1,potential_size2);
for i=1:potential_size2
    votes3(i)=sum(adj_mat3(i,:));
end
[~,voters3]=sort(votes3,'descend');
voters3=potential_set2(voters3);

for seq3=1:potential_size2*0.2

ele_this3=voters3(seq3);

[R_raw,t_raw]=Horn_minimal_solver(pts_3d([ele_this,ele_this2,ele_this3],:),pts_3d_([ele_this,ele_this2,ele_this3],:));

    res=zeros(1,potential_size2);count=0;consensus_set=zeros(1,1);
    pot_set=[potential_set2,ele_this,ele_this2];
    for i=1:potential_size2+2
        res(i)=norm(R_raw*((pts_3d(pot_set(i),:)))'+t_raw-(pts_3d_(pot_set(i),:))');
        if res(i)<=6*noise 
            count=count+1;
            consensus_set(count)=pot_set(i);
        end
    end
    
    if count>best_size
        best_set=consensus_set;
        best_ele=[ele_this,ele_this2,ele_this3];
        best_size=length(best_set);
        max_num=log(1-0.99)/log(1-(best_size/n_ele));
        max_num2=log(1-0.99)/log(1-((best_size-1)/potential_size));
        max_num3=log(1-0.99)/log(1-((best_size-2)/potential_size2));
    end
    
if seq3>=max_num3
    break
end

end
    
end

if seq2>=max_num2
    break
end
    
end


end

if seq>=max_num
    break
end

end
    

% end
    

opt_set=best_set;

q_=zeros(3,1);
p_=zeros(3,1);

opt_set=unique(opt_set);

len_o=length(opt_set);

for i=1:len_o

q_(1)=q_(1)+pts_3d_(opt_set(i),1);
q_(2)=q_(2)+pts_3d_(opt_set(i),2);
q_(3)=q_(3)+pts_3d_(opt_set(i),3);

p_(1)=p_(1)+pts_3d(opt_set(i),1);
p_(2)=p_(2)+pts_3d(opt_set(i),2);
p_(3)=p_(3)+pts_3d(opt_set(i),3);

end

p_=p_/len_o;
q_=q_/len_o;

s_best=1;

H=zeros(3,3);
for i=1:len_o
    H=H+(pts_3d(opt_set(i),:)'-p_)*(pts_3d_(opt_set(i),:)'-q_)';
end

[U,~,V]=svd(H);

R_opt=V*U';

t_opt=q_-s_best*R_opt*p_;


best_set=zeros(1,1);

res=zeros(1,n_ele);s_set=zeros(1,1);coun=0;
    
    for i=1:n_ele
        res(i)=(s_best*R_opt*((pts_3d(i,:)))'+t_opt-(pts_3d_(i,:))')'*(s_best*R_opt*((pts_3d(i,:)))'+t_opt-(pts_3d_(i,:))');
        if sqrt(res(i))<=6*noise
            coun=coun+1;
            best_set(coun)=i;
        end
    end


q_=zeros(3,1);
p_=zeros(3,1);

best_size=length(best_set);

for i=1:best_size

q_(1)=q_(1)+pts_3d_(best_set(i),1);
q_(2)=q_(2)+pts_3d_(best_set(i),2);
q_(3)=q_(3)+pts_3d_(best_set(i),3);

p_(1)=p_(1)+pts_3d(best_set(i),1);
p_(2)=p_(2)+pts_3d(best_set(i),2);
p_(3)=p_(3)+pts_3d(best_set(i),3);

end

p_=p_/best_size;
q_=q_/best_size;

s_best=1;

H=zeros(3,3);
for i=1:best_size
    H=H+(pts_3d(best_set(i),:)'-p_)*(pts_3d_(best_set(i),:)'-q_)';
end

[U,~,V]=svd(H);

R_opt=V*U';

t_opt=q_-s_best*R_opt*p_;

time=toc();

iteration=seq;

R_error=rotError(R_gt,R_opt)*180/pi;

t_error=norm(t_opt - t_gt');

display(['Rotation error in degrees:']);
R_error

display(['Translation error in meters:']);
t_error

display(['Runtime in seconds:']);
time

