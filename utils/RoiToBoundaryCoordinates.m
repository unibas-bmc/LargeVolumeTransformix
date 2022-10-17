function [coords] = RoiToBoundaryCoordinates(roi)
%roi_to_boundary_coordinates Takes a roi defined as [x1,x2,y1,y2,z1,z2] and
%gives back the coordinates (x,y,z) of the 8 corners of the roi

coords = [roi(1),roi(3),roi(5);
    roi(2),roi(3),roi(5);
    roi(1),roi(4),roi(5);
    roi(2),roi(4),roi(5);
    roi(1),roi(3),roi(6);
    roi(2),roi(3),roi(6);
    roi(1),roi(4),roi(6);
    roi(2),roi(4),roi(6)];

end

