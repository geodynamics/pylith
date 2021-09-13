#if !defined(sliderblockapp_hh)
#define sliderblockapp_hh

#include <portinfo>

#include <string> // HASA std::string

class Formulation;

class SliderBlockApp {
public:

    SliderBlockApp(void);
    ~SliderBlockApp(void);

    int run(int argc,
            char* argv[]);

private:

    void _getOptions(void);

    void _printHelp(void);

    void _initialize(void);

    void _solve(void);

    std::string _optEquations;
    std::string _optRupture;
    std::string _optFriction;
    std::string _outputFilename;
    Formulation* _formulation;
    bool _showHelp;

private:

    // Not implemented
    SliderBlockApp(const SliderBlockApp&);
    const SliderBlockApp& operator=(const SliderBlockApp&);

};

#endif

// End of file
