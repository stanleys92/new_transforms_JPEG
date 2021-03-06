This software was created by Independent JPEG Group.
- modified 2016-2017 by Stanislav Svoboda
- supervisor: Ing. David Barina, Ph.D.

This software enables to use even other transform than DCT, which is used in original JPEG standard. You can choose from 12 different transforms like local discrete cosine transform, discrete Chebyshev transform, discrete wavelet transforms, ... Experiments on dataset of images shows that the local discrete cosine transform performs better than the discrete cosine transform. Besides, it removes blocking artifacts. Blocking artifacts can be reduced even by wavelet transforms but they dont perform better than local discrete cosine transform. The results of PSNR and SSIM for few transforms can be seen on images below. With this extension comes a restriction. These new transforms were only optimized for block size 8x8, for other parameters of block size new transforms will not work. 

modifications:
  - cjpeg.c
      - changed command line arguments
      - added loading whole image to 2D array
  - djpeg.c
      - changed command line arguments
      - added loading whole image to 2D array
  - jcdctmgr.c
      - added init whole image 2D array
  - jddctmgr.c
      - added init whole image 2D array
  - jfdctint.c
      - added new transforms
  - jidctint.c
      - added new transforms

additions:
  - jfdctint.h
      - added functions for new transforms (FORWARD)
  - jidctint.h
      - added functions for new transforms (INVERSE)
  - trans.h
      - added functions for allocation of whole image


- example usage:
	- compression: ./cjpeg -outfile testimg.jpeg -quality 70 -nosmooth -scale 1/1 -trans ldct testimg.bmp
	- decompression: ./djpeg -bmp -outfile newtestimg.bmp -nosmooth -trans ildct testimg.jpeg
- notes and restrictions:
	- when compressing the image, you must run with these parameters: "-outfile output_file_name.jpeg", "-nosmooth" and "-scale 1/1"
	- these transforms can be used for compression:
		- "-trans dct" - discrete cosine transform
		- "-trans dst" - discrete sine transform
		- "-trans dht" - discrete Hartley transform
		- "-trans wht" - Walsh-Hadamard transform
		- "-trans wav" - discrete wavelet transform - separable CDF 5/3
		- "-trans wav97" - discrete wavelet transform - separable CDF 9/7
		- "-trans wavrb" - discrete wavelet transform - nonseparable CDF 5/3
		- "-trans wavrb97" - discrete wavelet transform - nonseparable CDF 9/7
		- "-trans wavall" - discrete wavelet transform - modified CDF 5/3
		- "-trans wavall97" - discrete wavelet transform - modified CDF 9/7
		- "-trans dcht" - discrete Chebyshev transform
		- "-trans ldct" - local discrete cosine transform
		- "-trans smrt" - mapped real transform
		- "-trans adct" - approximate discrete cosine transform

	- when decompressing the image, you must run with these parameters: "-bmp" "-outfile output_file_name.bmp" and "-nosmooth"
	- these transforms can be used for decompression:
		- "-trans dct" - discrete cosine transform
		- "-trans dst" - discrete sine transform
		- "-trans dht" - discrete Hartley transform
		- "-trans wht" - Walsh-Hadamard transform
		- "-trans wav" - discrete wavelet transform - separable CDF 5/3
		- "-trans wav97" - discrete wavelet transform - separable CDF 9/7
		- "-trans wavrb" - discrete wavelet transform - nonseparable CDF 5/3
		- "-trans wavrb97" - discrete wavelet transform - nonseparable CDF 9/7
		- "-trans wavall" - discrete wavelet transform - modified CDF 5/3
		- "-trans wavall97" - discrete wavelet transform - modified CDF 9/7
		- "-trans dcht" - discrete Chebyshev transform
		- "-trans ldct" - local discrete cosine transform
		- "-trans smrt" - mapped real transform
		- "-trans adct" - approximate discrete cosine transform



![Comparison of DCT, WHT, LDCT, SMRT and CDF 9/7, using PSNR metric](https://cloud.githubusercontent.com/assets/9696012/26224697/dfcb5c28-3c23-11e7-9c84-fbbbe01e1de9.png)
![Comparison of DCT, WHT, LDCT, SMRT, DChT and CDF 9/7, using SSIM metric](https://cloud.githubusercontent.com/assets/9696012/26224701/e38bbdf8-3c23-11e7-97a3-ad259a409dee.png)
