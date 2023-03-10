/***************************************************************************
**                                                                        **
**  URANOS - Ultra RApid Neutron-Only Simulation                          **
**  designed for Environmental Research                                   **
**  Copyright (C) 2015-2023 Markus Koehli,                                **
**  Physikalisches Institut, Heidelberg University, Germany               **
**                                                                        **
****************************************************************************/

#ifdef _WIN32
    #include <io.h>
    #define access _access
#elif __linux__
    #include <inttypes.h>
    #include <unistd.h>
    #define __int64 int64_t
    #define _close close
    #define _read read
    #define _lseek64 lseek64
    #define _O_RDONLY O_RDONLY
    #define _open open
    #define _lseeki64 lseek64
    #define _lseek lseek
    #define stricmp strcasecmp
#endif

#ifdef _WIN32
#pragma warning (disable: 4018 4100 4101 4189 4305)
#elif __linux_
#endif



#include "Toolkit.h"

#include "mainwindow.h"
//#include "customSplashScreen.h"
#include <QApplication>
#include <QPixmap>
#include <QSplashScreen>
#include <QCoreApplication>
#include <QTimer>
#include <QThread>
//#include <shellscalingapi.h>

string versionStringMain = "v1.11 (09.03.2023)";

class I : public QThread
{
    public:
    static void sleep(unsigned long secs)
    {
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

    QPixmap splashImage(":/resources/splashScreen.png");

    QSplashScreen splash(splashImage,Qt::WindowStaysOnTopHint);

    #ifdef _WIN32
    #elif __linux__
    QFont custom_font = QFont();
    //custom_font.setWeight(20);
    //custom_font.setPixelSize(9);
    custom_font.setPointSize(8);
    a.setFont(custom_font, "QLabel");
    a.setFont(custom_font, "QTabWidget");
    a.setFont(custom_font, "QPushButton");
    a.setFont(custom_font, "QLineEdit");
    a.setFont(custom_font, "QGroupBox");
    a.setFont(custom_font, "QCheckBox");
    a.setFont(custom_font, "QRadioButton");
    a.setFont(custom_font, "QTableView");
    a.setFont(custom_font, "QTableWidgetItem");
    a.setFont(custom_font, "QStandardItemModel");
    a.setFont(custom_font, "QStandardItem");
    a.setFont(custom_font, "QString");
    a.setFont(custom_font, "QSpinBox");
    #endif

    /* To intercept mousclick to hide splash screen. Since the
    splash screen is typically displayed before the event loop
    has started running, it is necessary to periodically call. */

    a.processEvents();


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

        // starts URANOS with the source control panel
        if ((std::string(argv[1])=="tnw"))
        {
            w.activateThermalSkyEvaporation();
        }

        // starts the batch run for n monoenergetic runs with neutrons directly downward
        if ((std::string(argv[1])=="detectorBatchRun"))
        {
            w.activateThermalSkyEvaporation();
            w.activateDetectorBatchRun();
        }

        // starts the batch run for n monoenergetic runs with neutrons having a random angular distribution
        if ((std::string(argv[1])=="detectorBatchRun2"))
        {
            w.activateThermalSkyEvaporation();
            w.activateDetectorBatchRun2();
        }

        // starts the batch run for n monoenergetic runs with batch run on angles
        if ((std::string(argv[1])=="detectorAngleBatchRun"))
        {
            w.activateThermalSkyEvaporation();
            w.activateDetectorAngleBatchRun();
        }

        if ((std::string(argv[1])=="--version") || ((std::string(argv[1])=="version")))
        {
            cout<<versionStringMain<<endl;
            return 0;
        }

        if ((std::string(argv[1])=="--buildversion") || ((std::string(argv[1])=="buildversion")))
        {
            cout<<QSysInfo::buildAbi().toStdString()<<endl;
            return 0;
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
        if (((std::string(argv[1])=="noGUI") || (std::string(argv[1])=="config")) && (std::string(argv[2])!=""))
        {
            if (std::string(argv[1])=="noGUI") disableGUI = true;
            configFilePath = std::string(argv[2]);
            std::replace(configFilePath.begin(), configFilePath.end(), '\\', '/'); cout<<" +";

            if (!(access(configFilePath.c_str(), 0) == 0))
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
        if (((std::string(argv[1])=="noGUI") || (std::string(argv[1])=="config")) && (std::string(argv[2])!=""))
        {
            if (std::string(argv[1])=="noGUI") disableGUI = true;
            configFilePath = std::string(argv[2]);
            std::replace( configFilePath.begin(), configFilePath.end(), '\\', '/'); cout<<" +";

            ifstream input_stream(configFilePath,ios::in);

            if (!(access(configFilePath.c_str(), 0) == 0))
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

                     if (!(access(configFilePath.c_str(), 0) == 0))
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

    if (!silent)
    {
        cout<<"Welcome to URANOS "<<versionStringMain<<endl;
    }

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
        a.closeAllWindows();
        return 0;
    }
    else
    {
        if (configFilePath.length() > 1)
        {
            cout<<"Using inline config file"<<endl;
            w.setConfigFilePath(configFilePath);
        }
        w.getGraphConfig();
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
