function a = vecangle360(v1,v2,n)

%v1=xyz position at t1
%v2=xyz position at t2
%n=plane of reference for clockwise (positive) and counterclockwise (neg)
%relative to v1

x = cross(v1,v2);
c = sign(dot(x,n)) * norm(x);
a = atan2d(c,dot(v1,v2));
end