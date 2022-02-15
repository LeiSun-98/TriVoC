function [pts_3d,pts_3d_,R,t]=Build_Scenario_A(n_ele,noise,outlier_ratio,scale_gt,show_figure)
%% synthetic data
% %set length of the cube
% Length=1;
% 
% n_pts=n_ele;
% 
% %Create all points
% pts_3d=Length*(rand(n_pts,3)-0.5);

%% real data (bunny)
% load source point cloud
pcd = pcread('bunny.ply'); %data/bunny T-rex_high_0
pcd_down = pcdownsample(pcd, 'random', 0.029); %0.029 for bunny 0.0065 for Ar 0.0001 for Lucy 0.006 for T-rex 

pts_3d = 6.3*pcd_down.Location; % source points %6.3 for bunny 0.0066 for Ar 0.00063 for Lucy 0.005 for T-rex 

pts_3d=double(pts_3d(1:n_ele,:));

pts_3d=pts_3d-0.5*[max(pts_3d(:,1))+min(pts_3d(:,1)), max(pts_3d(:,2))+min(pts_3d(:,2)), max(pts_3d(:,3))+min(pts_3d(:,3))];

n_pts=n_ele;

%% Transformation
%Generate random rotation
axis1= rand(1,3)-0.5;

axis1=axis1/norm(axis1);

angle=2*pi*rand(1);

aa = angle * axis1;

if norm(aa) < 2e-16
         R=eye(3);
else k = aa / norm(aa);
K1=[0, -k(3), k(2); k(3), 0, -k(1); -k(2), k(1), 0];

%construct a rotation matrix from an axis angle representation
R = eye(3) + sin(norm(aa)) * K1 + (1 -cos(norm(aa))) * K1 * K1;
end

%Generate random translation

t_raw=rand(1,3)-0.5;

t=0.0*[1,1,1]+(3*rand(1))*t_raw/norm(t_raw);


%transform by R & t
d_med= scale_gt * pts_3d * R';

    for i=1:n_pts
    pts_(i,:) = d_med(i,:) + t ;
    end

%add noise
pts_3d_=pts_+noise*randn(n_ele,3);

%add outliers

for i=1:round(n_ele*outlier_ratio)
    
    for iii=1:1e+18
    rand_vec=2*scale_gt*rand(1,3)-scale_gt;
        if norm(rand_vec)<= 1*scale_gt
            break
        end
    end
    
pts_3d_(i,:)=[0.0,0.0,0.0]+t+rand_vec;%5*mean(pts_3d_(i,:))*rand(1,3);

end


%% show figure

if show_figure==1
    
figure(1);

pc1=pointCloud(pts_3d(1:n_ele,:));
pc2_out=pointCloud(pts_3d_(1:round(outlier_ratio*n_ele),:));
pc2_in=pointCloud(pts_3d_(1+round(outlier_ratio*n_ele):n_ele,:));

pcshow([pc1.Location(:,1),pc1.Location(:,2),pc1.Location(:,3)],[0 0 1],'MarkerSize',90);

hold on;

pcshow([pc2_out.Location(:,1),pc2_out.Location(:,2),pc2_out.Location(:,3)],[1 0 0],'MarkerSize',200);

hold on;

pcshow([pc2_in.Location(:,1),pc2_in.Location(:,2),pc2_in.Location(:,3)],[0 1 0],'MarkerSize',200);

hold on;

for i=round(n_ele*outlier_ratio)+1:n_ele
    
    plot3([pts_3d(i,1),pts_3d_(i,1)],[pts_3d(i,2),pts_3d_(i,2)],[pts_3d(i,3),pts_3d_(i,3)],'g','LineWidth',2.5);
    
end

for i=1:round(n_ele*outlier_ratio)
    
    pp=plot3([pts_3d(i,1),pts_3d_(i,1)],[pts_3d(i,2),pts_3d_(i,2)],[pts_3d(i,3),pts_3d_(i,3)],'r','LineWidth',1);

    pp.Color(4) = 0.15;
    
end

grid off;
axis off;


end

end