#define hypot _hypot  //this is only due to some changes in MSCV2010

//#include <math.h>

//#ifndef hypot
//#define hypot _hypot
//#endif
#include <io.h>

#include "Toolkit.h"

#include "mainwindow.h"
#include "customSplashScreen.h"
#include <QApplication>
#include <QPixmap>
#include <QSplashScreen>
#include <QCoreApplication>
#include <QTimer>
#include <QThread>
//#include <shellscalingapi.h>

class I : public QThread
{
public:
static void sleep(unsigned long secs) {
QThread::sleep(secs);
}
};


//BOOL dpi_result = SetProcessDPIAware();
//BOOL dpi_result = SetProcessDpiAwareness();

// MAIN method
// the main URANOS code has been transferred to one window application w.
// the window can be closed and URANOS continues to run with the command line output.
// the window w is used to generate the settings, which are saved in the uranos.cfg and in the geometryConfig.dat
// upon clicking 'simulate' URNAOS loads the settings from these files and executes the main simulation which loops over n initial neutrons
// during the simulation it shows the live output of many different parameters including the spatial distribution
int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    //QApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
    //QFont f;
    //int defaultFontSize = f.pointSize();
    //QScreen *screen = QGuiApplication::primaryScreen();
    //QRect  screenGeometry = screen->geometry();
    //int height = screenGeometry.height();
    //int width = screenGeometry.width();
    //cout<<"height"<<height;
    //cout<<"width"<<width;
    bool disableGUI = false;
    bool silent = false;
    string configFilePath = "";  

    QPixmap splashImage("splashScreen.png");

    QSplashScreen splash(splashImage,Qt::WindowStaysOnTopHint);


    /* To intercept mousclick to hide splash screen. Since the
    splash screen is typically displayed before the event loop
    has started running, it is necessary to periodically call. */
    //app.processEvents();

    //qApp->processEvents();
    a.processEvents();


    //splash->showStatusMessage(QObject::tr("Initializing…"));
    //splash->showStatusMessage(QObject::tr("Loading something…"));

    rootlogon();

    MainWindow w;


    // these are the command line options for starting URANOS
    if (argc>1)
    {
        // activates thermal neutron transport
        if ((std::string(argv[1])=="thermal"))
        {
            w.activateThermal();
        }

        // hides the GUI and starts from the config file in the URANOS folder
        if ((std::string(argv[1])=="noGUI"))
        {
            disableGUI = true;
        }

        // starts the batch run for n monoenergetic runs with neutron downward
        if ((std::string(argv[1])=="detectorBatchRun"))
        {
            w.activateThermalSkyEvaporation();
            w.activateDetectorBatchRun();
        }
    }

    if ((argc>2)&&(argc<4))
    {
       
        // activates the source control panel which allows other sources than the cosmic spectrum
        if ((std::string(argv[1])=="nuclear")&&(std::string(argv[2])=="warfare"))
        //if (true)
        {
            w.activateSkyEvaporation();
        }

         // hides the GUI and starts from the config file in the folder specified by the second parameter
        if ((std::string(argv[1])=="noGUI")&&(std::string(argv[2])!=""))
        {
            disableGUI = true;
            configFilePath = std::string(argv[2]);
            std::replace(configFilePath.begin(), configFilePath.end(), '\\', '/'); cout<<" +";


            if (!access( configFilePath.c_str(), 0 ) == 0 )            
            {
                cout<<"Config File not found"<<endl;
                return 0;
            }
        }
    }

    if  ((argc>3)&&(argc<5))
    {
        // activates the source control panel which allows other sources than the cosmic spectrum including thermal neutron transport
        if ((std::string(argv[1])=="thermal")&&(std::string(argv[2])=="nuclear")&&(std::string(argv[3])=="warfare"))
        //if (true)
        {
            w.activateThermalSkyEvaporation();
        }

        //  hides the GUI and starts from the config file in the folder specified by the second parameter and allows a third parameter
        if ((std::string(argv[1])=="noGUI")&&(std::string(argv[2])!=""))
        {
            disableGUI = true;
            configFilePath = std::string(argv[2]);
            std::replace( configFilePath.begin(), configFilePath.end(), '\\', '/'); cout<<" +";

            ifstream input_stream(configFilePath,ios::in);

            //if (!access( textHere.c_str(), 0 ) == 0 )

            if (!access(configFilePath.c_str(), 0 ) == 0)
            //if(input_stream.bad())
            {
                cout<<"Config File not found"<<endl;
                return 0;
            }
        }

        // disables command line output
        if (std::string(argv[3])=="silent")
        {
            w.beSilent();
            silent = true;
        }

        // activates the source control panel which allows other sources than the cosmic spectrum including thermal neutron transport
        if (std::string(argv[3])=="source")
        {
            w.activateThermalSkyEvaporation();
        }

        // activates thermal neutron transport
        if (std::string(argv[3])=="thermal")
        {
            w.activateThermal();
        }

        // activates the 2D batchrun over two different parameter sets like soil moisture and air humidity
        if ((std::string(argv[1])=="batchrun")&&(std::string(argv[2])!="")&&(std::string(argv[3])!=""))
        {
            w.activateParameterBatchRun(atoi(argv[2]),atoi(argv[3]));
            cout<<"Batchrun from Parameter "<<atoi(argv[2])<<" to "<<atoi(argv[3])<<endl;

            /*
            if ((std::string(argv[4])=="noGUI"))
            {
                disableGUI = true;
            }
            else
            {

            }
            */
        }
    }

    // same as above for arguments <5 used otherwise
    if ((argc>4)&&(argc<7))
    {
        if ((std::string(argv[1])=="batchrun")&&(std::string(argv[2])!="")&&(std::string(argv[3])!=""))
        {
            w.activateParameterBatchRun(atoi(argv[2]),atoi(argv[3]));
            cout<<"Batchrun from Parameter "<<atoi(argv[2])<<" to "<<atoi(argv[3])<<endl;

            if ((std::string(argv[4])=="noGUI"))
            {
                disableGUI = true;
            }
            else
            {

            }

            if (argc>5)
            {
                if ((std::string(argv[4])=="noGUI")&&(std::string(argv[5])!=""))
                {
                     disableGUI = true;

                     configFilePath = std::string(argv[5]);
                     std::replace( configFilePath.begin(), configFilePath.end(), '\\', '/'); cout<<" +";

                     ifstream input_stream(configFilePath,ios::in);

                     if (!access(configFilePath.c_str(), 0 ) == 0)
                     {
                         cout<<"Config File not found"<<endl;
                         return 0;
                     }
                }
            }
        }
    }


    //w.activateSkyEvaporation();
    //w.activateThermalSkyEvaporation();

    //disableGUI = true;

    // command line output logos
    if (false)
    {
        cout<<endl;
        cout<<"                                                                                                                                        .,*//*.      "<<endl;
        cout<<"                                                                                                                                    .,/(((((*.       "<<endl;
        cout<<"                                                              .                                           ....*((((*,...          .*(//*/((/.        "<<endl;
        cout<<" .*/**//,       .,.   .,,*/((((((((((((((((/,.              ,(#(,         /(((((((*.      .*(((((((,    ,/#(, ,(%%(, ,(#/,       ,//**((((/,         "<<endl;
        cout<<" .*#%%%(,      .,,*****,  .,/#%%%#/,...,*(%%%#/.           ,(%%%(,        .,/#%%%%%(,      ..*(%#/,   ...*##/**(##(*,/##*...    ,***//*/(#(,         "<<endl;
        cout<<" .*#%%%(,       ,,******,   *#%%%#,      ,#%%%#/.         ./#%%%#/.         ./#%%%%%#/.       *#/.  ,//. ./#(,,/##/,,(#/. ./,  .,**///#(/*,.         "<<endl;
        cout<<" .*#%%%(,       ,/(#%%#/*   *#%%%#,      ,(%%%%(.        ./#(#%%%#/.        ./#//#%%%%#*.     *(/. ,/#%(,,/##//(##(//#(*.,(%/. .*/(/((*.             "<<endl;
        cout<<" .*#%%%(,         ./#(,     *#%%%#,     ,/#%%%(*        .*(/,*#%%%#*        ./(*.,(#%%%#(,    *(/. ..*(####*. .*##*  ./####(.. .*#%%%##(//*,..       "<<endl;
        cout<<" .*#%%%(,          *#(.     *#%%%#(///((###((*.         *(/. ./#%%%(,       ./(*  .*(%%%%#/.  *(/. *, ,(%#*   .*##*   .*#%(,,*, ./#%%%%%%%%%%#(/,    "<<endl;
        cout<<" .*#%%%(,          *(/.     *#%%%#/*,*/#%#(,.          ,((,   .(%%%%(,      ./(*    ./#%%%%#*.*(/. (/.,(%(,   .*##*    ,(%(,/(*   ,/(#%%%%%%%%%%#/,  "<<endl;
        cout<<" ./#%%%(,          *(/.     *#%%%#,    ,(%%#(,        ,(#(/***/(#%%%#/.     ./(*      ,(#%%%#((#/. *..*(%#*   .*##*    *#%(,.*,...   .,*/((##%%%%#/. "<<endl;
        cout<<" ./#%%%(*......   .*(/,    .*#%%%#*.   ,(%%%%(,      ,(#(/******/#%%%#/.    ./(*.      .*#%%%%%#/. ,,*(###(*. ./##*. .*(###(,,../#/.       ..,/#%%(, "<<endl;
        cout<<"../#%%%(*..........*(/,...../#%%%#*...../#%%%#/,,/,.*##/.........*(%%%%(,...,(#/..........*(#%%#(,.,/##*.,/##//(##(//##/..*#/,../#%#/,........,(%#/.."<<endl;
        cout<<"../#%%%(*..........*(/,..,,*(#%%%#(*,,..*(%%%%#(((/(#%%#/*,.....,*(%%%%#(/**(#%#/*,.........*#%#(,...*,..*(#(**(##(**(#(*..,..../#%###/*,...,*(#(*..."<<endl;
        cout<<"../#%%%(*..........*(/,..*///////////*...,*((#((/*////////*....,//////////////////*..........,/(/......,*(%#*,*#%%#*,*#%(*,.....*//,,*/(((((((/,....."<<endl;
        cout<<".,/#%%%#*.........,/#/,.................................................................................,*//,,/#%%#/,,//*,..........................."<<endl;
        cout<<"  ,(%%%#/,       .*((,.  .......,,,,,,,*******////////((((((((((((((((((((((((((((((((((((((((((((((((((((((((######((((((((((((((((((((/////****,,.."<<endl;
        cout<<"   ,(#%%%#/*,..,*(#(,       .......,,,,,,*******/////////(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((//////***,,,.. "<<endl;
        cout<<"     ,*((######((*,             .........,,,,,,,,,****************************************************************************************,,,,,...   "<<endl;
        cout<<"           ..                          .........................,,,,,,,,,,,,,,,,,,,,,,,........................................................      "<<endl;
    }
	
	
    if (!silent)
    {
        cout<<"                                                                       ./%%/** "<<endl;
        cout<<"**/,         /((((((((/.        ((     /((((.    ((((   /( (#( (/    ./%%/*    "<<endl;
        cout<<"*%%,   *****  *%%,   ,%%/      /%%/      /%%%/    */  /..#,/#/,#../  *%%*      "<<endl;
        cout<<"*%%,    ./(   *%%,   /%%*     */*%%*     /*,#%#,  */ .*##* *#* *##*. %%%%/*    "<<endl;
        cout<<"*%%,    .*/   *%%(((%%*      *#  #%#.    /* .%%%/.*/  .(#  *#*  #(.   /%%%%%%/ "<<endl;
        cout<<"/%%,    .*/   *%%*  ,%%*    *#((((%%(    /*   *%%%%/  .##. *#* .##.        /#%/"<<endl;
        cout<<"/%%,    .*/   /%%*   #%%   *#(    *%%(   /*     (%%/ ,#*,#.*#*.#.*#, #%/    %%/"<<endl;
        cout<<"/%%,.....*/..//%%*..../#%%/#(....../%%/../*......,%/   ,(#,###,#(,   /%%%%%%/  "<<endl;
        cout<<",%%/    */, ...**.....**%%%(.......((((((((......((((..((.(###(.((..((((%%%/.. "<<endl;
        cout<<" ,#%#%%#/,     ........//////////////////////////////////////////////////***..."<<endl;
        cout<<"                        ...................................................    "<<endl;
    }

    // Main Window control
    if (disableGUI)
    {
        w.disabledGUIRun(configFilePath);
        //if (w.isFinished())  return 0;
        I::sleep(0.5);

        w.close();

        return a.closeAllWindows();
    }
    else
    {
        splash.show();
        if(splashImage.isNull());
        else  I::sleep(1.5);

        w.show();

        splash.finish(&w);
    }


    //w.close();
    //w.cosmicNSimulator(&w);


    return a.exec();
}
