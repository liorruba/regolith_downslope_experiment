function [x_comp, y_comp, z_comp] = roty_deg(alpha, x_comp, y_comp, z_comp)
    % Rotation around the y axis.
    
    MAT = [cosd(alpha) 0 sind(alpha); 0 1 0; -sind(alpha) 0  cosd(alpha)] * [x_comp'; y_comp'; z_comp'];
    x_comp = MAT(1,:)';
    y_comp = MAT(2,:)';
    z_comp = MAT(3,:)';
end