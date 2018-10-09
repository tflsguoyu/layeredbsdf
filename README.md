# Position-Free Monte Carlo Simulation for Arbitrary Layered BSDFs

[Yu Guo](https://www.ics.uci.edu/~yug10/), [Miloš Hašan](http://miloshasan.net/), [Shuang Zhao](https://shuangz.com/). 
In ACM Transactions on Graphics (SIGGRAPH Asia 2018). 
[[Project page]](https://shuangz.com/projects/layered-sa18/)

<img src="https://www.ics.uci.edu/~yug10/projects/SiggraphAsia2018/git-readme/images/teaser.jpg" width="500px">

## Overview
This is a branch of the Mitsuba (0.6.0) renderer (official repo: https://github.com/mitsuba-renderer/mitsuba)

## Install

  ### Linux (Preferred, Tested on Ubuntu 16.04)
   - `$ git clone https://github.com/tflsguoyu/layeredbsdf.git`
   - `$ cd layeredbsdf/`
   - `$ mv config_linux.py config.py`
   - `$ sudo apt install build-essential scons mercurial libpng12-dev libjpeg-dev libilmbase-dev libxerces-c-dev libboost-all-dev libopenexr-dev libglewmx-dev libxxf86vm-dev libpcrecpp0v5 libeigen3-dev libfftw3-dev`
   - `$ scons -j x` (x = # of cpu cores)
   - `source setpath.sh`
   Now you can render scenes
   - `$ mitsuba xxx.xml`
   
  ### Windows (Tested on Windows 10 x64)
   - install visual studio 2017
   - clone this git to local folder
   - go to folder ..\layeredbsdf\
   - rename config_windows.py to config.py
   - download [dependencies](https://www.ics.uci.edu/~yug10/projects/SiggraphAsia2018/git-readme/dependencies.zip)
   - 
   - 
   
## Example scenes (click image to download scene files)
(Will have more..)

<a href="https://www.ics.uci.edu/~yug10/projects/SiggraphAsia2018/git-readme/scenes/teaser.zip">
  <img src="https://www.ics.uci.edu/~yug10/projects/SiggraphAsia2018/git-readme/images/teaser.jpg" title="teaser" height="128px">
</a>
  
<a href="https://www.ics.uci.edu/~yug10/projects/SiggraphAsia2018/git-readme/scenes/figure2.zip">
  <img src="https://www.ics.uci.edu/~yug10/projects/SiggraphAsia2018/git-readme/images/figure2.jpg" title="figure2" height="128px">
</a>
  
<a href="https://www.ics.uci.edu/~yug10/projects/SiggraphAsia2018/git-readme/scenes/figure3.zip">
  <img src="https://www.ics.uci.edu/~yug10/projects/SiggraphAsia2018/git-readme/images/figure3.jpg" title="figure2" height="128px">
</a>  

## Scene file (.xml) explanation
 - `<scene version="0.6.0">` (Here using 0.6.0, but not 0.5.0)
 - `<integrator type="path_layered">` (This integrator is based on `path`)
 - `<bsdf type="multilayered"> ... </bsdf>` (BSDF type is `multilayered`, both our `uni-dir` and `bi-dir` methods are implemented here)
 - 
## Notes
 - Default precision in `config.py` is `single`. If you find too many warnings or even crashed when rendering scenes, you should try `double` precision instead. (Already provided in `config.py`)
 - Welcome to report bugs and leave comments (Yu Guo: tflsguoyu@gmail.com)
