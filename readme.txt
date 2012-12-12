how to use:

[map,chanmaps,maps,chans] = my_simpsal( img , param );

input: img... matrix image_height x image_width x [r g b x y z]
				r g b ... red, green, blue channels
				x y z ... coordinates of each pixel

	param ... parameter settings (see default_fast_param)

output: map ... 	final saliency map
	chanmaps ...	final saliency map for each channel
	maps ...	saliency map for each sub-channel, for each channel
	chans ...	feature maps for each sub-channel, for each channel





-------------------------------------------------------
based on 
Jonathan Harel, A Saliency Implementation in MATLAB: http://www.klab.caltech.edu/~harel/share/gbvs.php

with extenstions for 2.5D images (like kinect data)


