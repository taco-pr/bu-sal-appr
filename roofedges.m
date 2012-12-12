function [edge orient] = roofedges(nx,ny,nz,thresholds)

if nargin < 4
    thresholds = [0.5 0.95];
end
thresholds = sort(thresholds);

n1 = zeros(size(nx,1),size(nx,2),3);
n2 = zeros(size(nx,1),size(nx,2),3);
%a1 = zeros(size(nx,1),size(nx,2));
%a2 = zeros(size(nx,1),size(nx,2));


f = [1 0 0;1 0 0;1 0 0]/3;
n1(:,:,1) = imfilter(nx,f,'replicate','same','conv');
n1(:,:,2) = imfilter(ny,f,'replicate','same','conv');
n1(:,:,3) = imfilter(nz,f,'replicate','same','conv');
f = [0 0 1;0 0 1;0 0 1]/3;
n2(:,:,1) = imfilter(nx,f,'replicate','same','conv');
n2(:,:,2) = imfilter(ny,f,'replicate','same','conv');
n2(:,:,3) = imfilter(nz,f,'replicate','same','conv');
a1 = dot(n1,n2,3);

f = [1 0 0;1 0 0;1 0 0]'/3;
n1(:,:,1) = imfilter(nx,f,'replicate','same','conv');
n1(:,:,2) = imfilter(ny,f,'replicate','same','conv');
n1(:,:,3) = imfilter(nz,f,'replicate','same','conv');
f = [0 0 1;0 0 1;0 0 1]'/3;
n2(:,:,1) = imfilter(nx,f,'replicate','same','conv');
n2(:,:,2) = imfilter(ny,f,'replicate','same','conv');
n2(:,:,3) = imfilter(nz,f,'replicate','same','conv');
a2 = dot(n1,n2,3);


a1=1-(a1+1)/2;
a2=1-(a2+1)/2;

%edge = sqrt(a1.^2+a2.^2)/sqrt(2);
edge = max(a1,a2);
orient = atan(a2./a1)*180/pi;
orient(orient<0) = orient(orient<0)+180;

h=[];%zeros(size(edge,1)*size(edge,2));
for i=1:size(edge,1);
    h=[h (edge(i,1:size(edge,2)))];
end
h=hist(h,1000);
h=h/sum(h);
%figure;plot(h);hold on;

sumh = 0;
thresmin=0;thresmax=size(h,2);
for i=1:size(h,2);
    sumh=sumh+h(i);
    if sumh > thresholds(2)
        thresmax=i;
        break;
    elseif sumh > thresholds(1) && thresmin==0;
        thresmin=i;
    end
    
end
dif = (max(max(edge))-min(min(edge)));
thresmin=thresmin/size(h,2)*dif;
thresmax=thresmax/size(h,2)*dif;

edge = (edge-thresmin)/(thresmax-thresmin);
edge(edge>1)=1;
edge(edge<0)=0;




end
