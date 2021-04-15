// This is:
// A __cdecl x64 Python module can be imported as `import explotest`
// means you need rename from .dll to .pyd using pybind11
//

// The DLL will only build if DEBUG preprocessor is not set.
// pybind11 only works in release mode
// and all code in the project should also be in release



////////////////////////////////////////////////
///////////////// Python API //////////////////
///////////////////////////////////////////////


#pragma once

#define BUILDING_DLL
#include <limits> // std::numeric_limits
#include "pybind11/embed.h"
#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"
#include "pybind11/cast.h"
#include "pybind11/stl.h"
#include "pybind11/functional.h"
#include "pybind11/stl_bind.h"
#include "GeographicLib/Geodesic.hpp"


namespace py = pybind11;
namespace gd = GeographicLib;

gd::Geodesic geodesic = gd::Geodesic::WGS84(); // Geodesic currently in use, so dont need to export this class to python

void WGS84(){ // set current Geodesic class to WGS84
    geodesic = gd::Geodesic::WGS84();
}

void Geodesic(double a, double f) { // set current Geodesic class to specified a, f params
    geodesic = gd::Geodesic(a, f);
}

std::pair<double, double> Direct(double lat1, double lon1, double az1, double s12) {
    double lat2, lon2;
    geodesic.Direct(lat1, lon1, az1, s12, lat2, lon2);    
    return std::make_pair(lat2, lon2);    
}

double Inverse(double lat1, double lon1, double lat2, double lon2) {
    double s12;
    geodesic.Inverse(lat1, lon1, lat2, lon2, s12);
    return s12;
}

std::pair<double, double> InverseAngle(double lat1, double lon1, double lat2, double lon2) {
    double az1, az2;
    geodesic.Inverse(lat1, lon1, lat2, lon2, az1, az2);
    return std::make_pair(az1, az2);
}


// for Python use  __cdecl
PYBIND11_MODULE(geolib, m){

    // optional module docstring
    m.doc() = "geographiclib - python api pybind11";

    m.def("WGS84", &WGS84, "Set current elipsoid to WGS84");

    m.def("Geodesic", &Geodesic, "Set current elipsoid to specified a, f params",
        py::arg("equator radius"), py::arg("f flattening"));
    
    m.def("Direct", &Direct, "Solve the direct geodesic problem where the length of the geodesic is specified in terms of distance.",
        py::arg("lat1"), py::arg("lon1"), py::arg("az1"), py::arg("s12"));

    m.def("Inverse", &Inverse, "Solve the inverse geodesic problem return the length between two points.",
        py::arg("lat1"), py::arg("lon1"), py::arg("lat2"), py::arg("lon2"));

    m.def("InverseAngle", &InverseAngle, "Solve the inverse geodesic problem return the angles az1/az2 in the two points.",
        py::arg("lat1"), py::arg("lon1"), py::arg("lat2"), py::arg("lon2"));
}
