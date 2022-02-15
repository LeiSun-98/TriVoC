function error = getAngularError(R_pre, R_gt)

error = abs(acos((trace(R_pre'*R_gt)-1)/2));

end