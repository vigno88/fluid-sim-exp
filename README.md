# fluid-sim-exp
Bunch of fluid sim expirements to learn fluid simulation

## 2D

The first 2D fluid sim requires the skia library. Just put the skia git project in the same folder
and place the `libskia.so` in `/usr/lib/`. 

The skia project is build with:
```
bin/gn gen out/Static --args='is_official_build=true cc="clang" cxx="clang++" is_component_build=true skia_use_system_expat=false skia_use_system_libjpeg_turbo=false skia_use_system_libpng=false skia_use_system_libwebp=false skia_use_system_zlib=false skia_use_system_icu=false skia_use_system_harfbuzz=false skia_use_gl=true'

ninja -C out/Static
```

Then the initial project that uses skia is compiled using:
```
g++ -g -std=c++1z glfw_ship.cpp -lskia -ldl -lpthread -ljpeg -lfreetype -lz -lpng -lglfw -lfontconfig -lwebp -lwebpmux -lwebpdemux -lGL -Iskia -Lskia/out/Static/
```
