#include <portinfo>

#include "SliderBlockApp.hh"

#include "QuasistaticPrescribedSlip.hh"
// #include "QuasistaticSpontaneousRupture.hh"
#include "DynamicPrescribedSlip.hh"
// #include "DynamicSpontaneousRupture.hh"

#include "petsc.h"

#include <getopt.h> // USES getopt_long()
#include <sstream> // USES std::ostringstream, std::istringstream
#include <cassert> // USES assert()
#include <iostream> // USES std::cout

// --------------------------------------------------------------------------------------------------
SliderBlockApp::SliderBlockApp() :
    _optEquations("quasistatic"),
    _optRupture("prescribed_slip"),
    _outputFilename("sliderblock.h5"),
    _formulation(NULL)
{}


// --------------------------------------------------------------------------------------------------
SliderBlockApp::~SliderBlockApp(void) {
    delete _formulation;_formulation = NULL;
}


// --------------------------------------------------------------------------------------------------
int
SliderBlockApp::run(int argc,
                    char* argv[]) {
    PetscErrorCode err = PetscInitialize(&argc, &argv, NULL, NULL);CHKERRQ(err);

    _getOptions();
    if (_showHelp) {
        _printHelp();
        return 0;
    } // if

    _initialize();
    _solve();

    delete _formulation;_formulation = NULL;
    err = PetscFinalize();CHKERRQ(err);

    return 0;
} // run


// --------------------------------------------------------------------------------------------------
void
SliderBlockApp::_getOptions(void) {
    PetscErrorCode err = 0;

    PetscBool showHelp = PETSC_FALSE;
    err = PetscOptionsGetBool(NULL, NULL, "-show-help", &showHelp, NULL);
    _showHelp = PETSC_TRUE == showHelp;

    char equations[32];
    PetscBool setEquations = PETSC_FALSE;
    err = PetscOptionsGetString(NULL, NULL, "-equations", equations, sizeof(equations), &setEquations);
    if (setEquations) { _optEquations = equations; }

    char rupture[32];
    PetscBool setRupture = PETSC_FALSE;
    err = PetscOptionsGetString(NULL, NULL, "-rupture", rupture, sizeof(rupture), &setRupture);
    if (setRupture) { _optRupture = rupture; }

    char filename[32];
    PetscBool setFilename = PETSC_FALSE;
    err = PetscOptionsGetString(NULL, NULL, "-output-filename", filename, sizeof(filename), &setFilename);
    if (setFilename) { _outputFilename = filename; }
} // _getOptions


// --------------------------------------------------------------------------------------------------
void
SliderBlockApp::_printHelp(void) {
    std::cout
        << "Usage: sliderblocks "
        << "[-show-help] "
        << "[-equations quasistatic|dynamic] "
        << "[-rupture prescribed_slip|spontaneous_rupture] "
        << "[-output-filename FILENAME] "
        << "[PETSc options]\n\n"
        << "Defaults:\n"
        << "    -equations quasistatic\n"
        << "    -rupture prescribed_slip\n"
        << "    -output-filename sliderblock.h5\n"
        << std::endl;
}


// --------------------------------------------------------------------------------------------------
void
SliderBlockApp::_initialize(void) {
    delete _formulation;_formulation = NULL;
    if (( _optEquations == "quasistatic") && ( _optRupture == "prescribed_slip") ) {
        _formulation = new QuasistaticPrescribedSlip();
    } else if (( _optEquations == "dynamic") && ( _optRupture == "prescribed_slip") ) {
        _formulation = new DynamicPrescribedSlip();
#if 0
    } else if (( _optEquations == "quasistatic") && ( _optRupture == "spontaneous_rupture") ) {
        _formulation = new QuasistaticSpontaneousRupture();
    } else if (( _optEquations == "dynamic") && ( _optRupture == "spontaneous_rupture") ) {
        _formulation = new DynamicSpontaneousRupture();
#endif
    } else {
        std::ostringstream msg;
        msg << "Cannot parse simulation options.\n"
            << "  --equations=" << _optEquations << ". Choices: quasistatic, dynamic\n"
            << "  --rupture=" << _optRupture << ". Choices: prescribed_slip, spontaneous_rupture\n";
        throw std::runtime_error(msg.str().c_str());
    } // if/else

    assert(_formulation);
    _formulation->initialize(_outputFilename.c_str());
}


// --------------------------------------------------------------------------------------------------
void
SliderBlockApp::_solve(void) {
    assert(_formulation);
    _formulation->solve();
}


// End of file
