function chan = getchan(imgs,channel,param)

chan = {};

% iterate over center scales
for scaleI = 1 : length(imgs)  

  img = imgs{scaleI};

  % Intensity
  if ( channel == 'I' )    
    chan{end+1} = squeeze( mean(img(:,:,1:3),3) );
    
  % Color
  elseif ( channel == 'C' )    
    if ( size(img,3) >= 3 )
      lum = squeeze( mean(img(:,:,1:3),3) ) + 0.01;    
      chan{end+1} = squeeze(abs(img(:,:,3,:) - min( img(:,:,1,:) , img(:,:,2,:)) )) ./ lum; % b - y (r,g)
      chan{end+1} = squeeze(abs(img(:,:,1,:) - img(:,:,2,:)))  ./ lum; % r - g 
    end
    
  % Color (dlk-colorspace)  
  elseif ( channel == 'D' )
    if ( size(img,3) >= 3 )
      map = img(:,:,1:3);      
      for fr = 1 : size(img,4)
        map(:,:,:,fr) = rgb2dkl( map(:,:,:,fr) );
      end      
      for i = 1 : 3
        chan{end+1} = squeeze( map(:,:,i,:) );
      end
    end
    
  % Orientation  
  elseif ( channel == 'O' )
    img = squeeze( mean(img(:,:,1:3),3) );
    
    for angi = 1 : param.nGaborAngles
      map = img;
      for fr = 1 : size(img,3)
        f0 = myconv2(img(:,:,fr),param.gabor{angi}.g0);
        f90 = myconv2(img(:,:,fr),param.gabor{angi}.g90);
        map(:,:,fr) = attenuateBordersGBVS( abs(f0) + abs(f90) , 13 );
      end
      chan{end+1} = map;
    end
    
  % Range for Intensity 
  elseif ( channel == 'R' )
      
      chan{end+1} = img(:,:,4);
    
  % Range for Orientation
  elseif ( channel == 'S' ) 
    img = img(:,:,4);
    
    for angi = 1 : param.nGaborAngles
      map = img;
      for fr = 1 : size(img,3)
        f0 = myconv2(img(:,:,fr),param.gabor{angi}.g0);
        f90 = myconv2(img(:,:,fr),param.gabor{angi}.g90);
        
        map(:,:,fr) = abs(f0)+abs(f90);
        map(:,:,fr) = attenuateBordersGBVS( abs(f0) + abs(f90) , 13 );
      end
      chan{end+1} = map;
    end

  % Jump-Edges
  elseif ( channel == 'T' )
      
      chan{end+1} = jumpedges(img(:,:,4),param.JumpEdgesThres);
    
  % Roof-Edges
  elseif ( channel == 'U' )
      
      h= [1 0 0 0 -1; 1 1 0 -1 -1; 1 1 0 -1 -1; 1 1 0 -1 -1; 1 0 0 0 -1];
      [nx ny nz] = normals(img(:,:,5),img(:,:,6),img(:,:,7),h);
      chan{end+1} = roofedges(nx,ny,nz,param.RoofEdgesThres);

  % normals
  elseif ( channel == 'V' )
      
      h= [1 0 0 0 -1; 1 1 0 -1 -1; 1 1 0 -1 -1; 1 1 0 -1 -1; 1 0 0 0 -1];
      [nx ny nz] = normals(img(:,:,5),img(:,:,6),img(:,:,7),h);
      chan{end+1} = nx;
      chan{end+1} = ny;
      chan{end+1} = nz;

  % x y z as itself
  elseif ( channel == 'W' )
      
      chan{end+1} = img(:,:,5);
      chan{end+1} = img(:,:,6);
      chan{end+1} = img(:,:,7);

  % harris-corner on jumpedges
  elseif ( channel == 'X' )
      j = jumpedges(img(:,:,4),param.JumpEdgesThres);
      h = harris(j,param.HarrisCornerDetectionSigma);
      chan{end+1} = h;

      

  end  % end if chan == M

end % end scale loop
