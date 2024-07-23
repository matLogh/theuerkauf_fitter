# Theuerkauf fitter
Theuerkauf is based on the [HDTV code](https://github.com/janmayer/hdtv), that uses Theuerkauf peak model for it's fits. This model primarily made for to fit the gamma spectra, it is a Gaussian peak with options of adding left and right tail, as well as background step underneath the peak, detailed description is available in the HDTV repository. Purpose of this repo is to have a reliable, pure C++ fitter for this peak model, that can be integrated in a projects - e.g. for optimization purposes.

This code (re)introduces two classes `TheuerkaufPeak` and `TheuerkaufFitter`. Former contain the peak model with switched to option to switch on/off parameters, that can be in 4 states:
* `FREE`, free parameter to be fitted,
* `FIXED`, fixed to a given value,
* `SAME`, free parameter, but if peak is part of the `TheuerkaufFitter`, then given type parameter type (e.g. sigma) will have common parameter in the fit, e.i. all will share the same value,
* `NONE`, disable the parameter, do not use it in the fit.

`TheuerkaufFitter` is used to fit one or more Theuerkauf peaks with background, which is by default polynomial (polynomial order is adjustable) or can be supplied by user. After fitting, `Theuerkauf::Analyze(TH1&)` is useful to see the fitted spectrum with total fit function, separated total background function (including peak steps) and individual peaks. It also shows on residuals with 2sigma interval (2*sqrt(bin_content)) and 95% confidence interval calculated by ROOT fitter. 

For *optimization purposes*, the ```TheuerkaufPeak::GetFWxM(const double width_multiple)``` function is useful, as it calculates width including also the tails, not just Gaussian sigma. 
## Dependencies & install

CERN ROOT framework, tested on v6.24 compiled against C++17.

### Example code

Code with an example can be build as
```bash
git clone https://github.com/matLogh/theuerkauf_fitter.git
cd theuerkauf_fitter
mkdir build
cd build
cmake ../example
make
```

### Include to your project

To add it to your project, you can add this project as a submodule 
```bash
git submodule add https://github.com/matLogh/theuerkauf_fitter.git <relative/path/where/to/place/submodule>
git submodule init
git submodule update
```
In your main ```CMakeLists.txt``` you need to have have ```${ROOT_LIBRARIES}``` and you should add
```cmake
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/<relative/path/where/to/place/submodule>)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/<relative/path/where/to/place/submodule>)
```
<!--  -->