function [R,t]=Horn_minimal(pts_3d,pts_3d_)

    v12=pts_3d(2,:)-pts_3d(1,:);
    X_axis=v12'/norm(v12);
    v13=pts_3d(3,:)-pts_3d(1,:);
    v23=cross(v12,v13);
    Y_axis=v23'/norm(v23);
    Z_axis=cross(X_axis,Y_axis);
    
    v12=pts_3d_(2,:)-pts_3d_(1,:);
    X_axis_=v12'/norm(v12);
    v13=pts_3d_(3,:)-pts_3d_(1,:);
    v23=cross(v12,v13);
    Y_axis_=v23'/norm(v23);
    Z_axis_=cross(X_axis_,Y_axis_);

    R=[X_axis_,Y_axis_,Z_axis_]*[X_axis,Y_axis,Z_axis]';
    
q_=zeros(3,1);
p_=zeros(3,1);
    
for i=1:3

q_(1)=q_(1)+pts_3d_(i,1);
q_(2)=q_(2)+pts_3d_(i,2);
q_(3)=q_(3)+pts_3d_(i,3);

p_(1)=p_(1)+pts_3d(i,1);
p_(2)=p_(2)+pts_3d(i,2);
p_(3)=p_(3)+pts_3d(i,3);

end

p_=p_/3;
q_=q_/3;

t=q_-R*p_;

end