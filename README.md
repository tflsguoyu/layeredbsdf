# Position-Free Monte Carlo Simulation for Arbitrary Layered BSDFs

[Yu Guo](https://tflsguoyu.github.io/), [Miloš Hašan](http://miloshasan.net/), [Shuang Zhao](https://shuangz.com/). 
In ACM Transactions on Graphics (SIGGRAPH Asia 2018). 

<img src="https://github.com/tflsguoyu/layeredbsdf_suppl/blob/master/github/images/teaser.jpg" height="300px"><img src="https://github.com/tflsguoyu/layeredbsdf_suppl/blob/master/github/images/green.jpg" height="300px"><img src="https://github.com/tflsguoyu/layeredbsdf_suppl/blob/master/github/images/yellow.jpg" height="300px"><img src="https://github.com/tflsguoyu/layeredbsdf_suppl/blob/master/github/images/blue.jpg" height="300px"><img src="https://github.com/tflsguoyu/layeredbsdf_suppl/blob/master/github/images/magenta.jpg" height="300px">

[[Paper](https://github.com/tflsguoyu/layeredbsdf_paper/blob/master/layeredbsdf.pdf)]
[[Code](https://github.com/tflsguoyu/layeredbsdf)]
[[Supplemental Materials](https://tflsguoyu.github.io/layeredbsdf_suppl/)]
[[Poster](https://github.com/tflsguoyu/layeredbsdf_poster/blob/master/layeredbsdf_poster.pdf)]
[[Fastforward on Siggraph Asia 2018](https://youtu.be/v5u6LYCN_PU) ([Slides](https://www.dropbox.com/s/zirw16peipdtq70/layeredbsdf_ff.pptx?dl=0))]

[[Presentation Slides on Siggraph Asia 2018](https://www.dropbox.com/s/i8h4h9jph1np3dt/layeredbsdf_main.pptx?dl=0)]
[[Two Minute Papers](https://youtu.be/Bv3yat484aQ)]

## Overview
This is a branch of the Mitsuba (0.6.0) renderer (official repo: https://github.com/mitsuba-renderer/mitsuba)

## Install

  ### Linux (Preferred, Tested on Ubuntu 16.04, 18.04)
   - `$ sudo apt update`
   - `$ sudo apt upgrade`
   - `$ sudo apt install git scons libboost-all-dev libpng-dev libjpeg-dev libopenexr-dev libxerces-c-dev libeigen3-dev libfftw3-dev libglewmx-dev freeglut3-dev`
   - `$ git clone https://github.com/tflsguoyu/layeredbsdf.git`
   - `$ cd layeredbsdf/`
   - `$ mv config_linux.py config.py`
   - `$ scons -j x` (x = # of cpu cores)
   - `$ source setpath.sh`
   Now you can render scenes
   - `$ mitsuba xxx.xml`
   
  ### Windows (Tested on Windows 10 x64)
   - install visual studio 2017
   - clone this git to local folder
   - go to folder ..\layeredbsdf\
   - rename config_windows.py to config.py
   - download [dependencies](https://github.com/tflsguoyu/layeredbsdf_suppl/blob/master/github/dependencies.zip)
   - 
   
## Examples in paper (click image to download scene files)

<a href="https://www.dropbox.com/s/jbethredc7r9742/teaser.zip?dl=0">
  <img src="https://github.com/tflsguoyu/layeredbsdf_suppl/blob/master/github/images/teaser.jpg" title="teaser" height="128px">
</a>
  
<a href="https://github.com/tflsguoyu/layeredbsdf_suppl/blob/master/github/scenes/figure2.zip">
  <img src="https://github.com/tflsguoyu/layeredbsdf_suppl/blob/master/github/images/figure2.jpg" title="figure2" height="128px">
</a>
  
<a href="https://github.com/tflsguoyu/layeredbsdf_suppl/blob/master/github/scenes/figure3.zip">
  <img src="https://github.com/tflsguoyu/layeredbsdf_suppl/blob/master/github/images/figure3.jpg" title="figure3" height="128px">
</a>  

<a href="https://github.com/tflsguoyu/layeredbsdf_suppl/blob/master/github/scenes/figure8.zip">
  <img src="https://github.com/tflsguoyu/layeredbsdf_suppl/blob/master/github/images/figure8.jpg" title="figure8" height="128px">
</a>  

<a href="https://www.dropbox.com/s/ds4ic1fw1cypqn7/figure11.zip?dl=0">
  <img src="https://github.com/tflsguoyu/layeredbsdf_suppl/blob/master/github/images/figure11.jpg" title="figure11" height="128px">
</a>  

</br>

<a href="https://github.com/tflsguoyu/layeredbsdf_suppl/blob/master/github/scenes/figure12t.zip">
  <img src="https://github.com/tflsguoyu/layeredbsdf_suppl/blob/master/github/images/figure12t.jpg" title="figure12t" height="128px">
</a>  

<a href="https://github.com/tflsguoyu/layeredbsdf_suppl/blob/master/github/scenes/figure12b.zip">
  <img src="https://github.com/tflsguoyu/layeredbsdf_suppl/blob/master/github/images/figure12b.jpg" title="figure12b" height="128px">
</a>  

<a href="https://www.dropbox.com/s/ssof650irx7blug/figure13.zip?dl=0">
  <img src="https://github.com/tflsguoyu/layeredbsdf_suppl/blob/master/github/images/figure13.jpg" title="figure13" height="128px">
</a>  

<a href="https://github.com/tflsguoyu/layeredbsdf_suppl/blob/master/github/scenes/figure14.zip">
  <img src="https://github.com/tflsguoyu/layeredbsdf_suppl/blob/master/github/images/figure14.jpg" title="figure14" height="128px">
</a>  

<a href="https://github.com/tflsguoyu/layeredbsdf_suppl/blob/master/github/scenes/figure15.zip">
  <img src="https://github.com/tflsguoyu/layeredbsdf_suppl/blob/master/github/images/figure15.jpg" title="figure15" height="128px">
</a>  

## Scene file (.xml) explanation
 - `<scene version="0.6.0">` (Here using 0.6.0, but not 0.5.0)
 - `<integrator type="path_layered">` (`path_layered` preferred, but can still use `path`, `volpath` or `volpath_simple` instead)
 - `<bsdf type="multilayered"> ... </bsdf>` (BSDF type is `multilayered`, both our `uni-dir` and `bi-dir` methods are implemented here)
 - 
## Notes
 - Default precision in `config.py` is `single`. If you find too many warnings or even it is crashed when rendering scenes, you should try `double` precision instead. (Already provided in `config.py`)
 - `conductor` and `dielectric` are not supported now, use `roughconductor (a=0.001)` and `roughdielectric (a=0.001)` instead.
 - Welcome to report bugs and leave comments (Yu Guo: tflsguoyu@gmail.com)
