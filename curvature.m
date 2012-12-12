function [H, K] = curvature(range_image)

r=range_image;

% for i = 1:1:20
%     r=imfilter(r,fspecial('gaussian',[5 5], 3));
% end

r=medfilt2(r,[3 3]);


hx = 1/50* [-2 -1 0 1 2;
        -2 -1 0 1 2;
        -2 -1 0 1 2;
        -2 -1 0 1 2;
        -2 -1 0 1 2];

hy = 1/50* [2  2  2  2  2;
      1  1  1  1  1;
      0  0  0  0  0;
     -1 -1 -1 -1 -1;
     -2 -2 -2 -2 -2];
 
hxx = 1/35* [2 -1 -2 -1 2;
        2 -1 -2 -1 2;
        2 -1 -2 -1 2;
        2 -1 -2 -1 2;
        2 -1 -2 -1 2];
    
hxy = 1/100* [-4 -2  0  2  4;
       -2 -1  0  1  2;
        0  0  0  0  0;
        2  1  0 -1 -2;
        4  2  0 -2 -4];

hyy = 1/35* [ 2  2  2  2  2;
       -1 -1 -1 -1 -1;
       -2 -2 -2 -2 -2;
       -1 -1 -1 -1 -1;
        2  2  2  2  2];

fx = imfilter(r,hx,'replicate','same','conv');
fy = imfilter(r,hy,'replicate','same','conv');
fxx = imfilter(r,hxx,'replicate','same','conv');
fyy = imfilter(r,hyy,'replicate','same','conv');
fxy = imfilter(r,hxy,'replicate','same','conv');

H = (fxx.*(1+fy.^2)+fyy.*(1+fx.^2)-2*fx.*fy.*fxy)./(2*(1+fx.^2+fy.^2).^(3/2));
K = (fxx.*fyy-fxy.^2)./(1+fx.^2+fy.^2).^2;

end

