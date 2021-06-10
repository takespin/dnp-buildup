#include "mainwindow.h"
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QGridLayout>
#include <QFileDialog>
#include <QSettings>
#include <QDebug>
#include <QCoreApplication>
#include "math.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
{
    // By default, the path to the parameter file and the data dile
    // are set to the application directory path.
    // Later on, we will overwrite the path by the last used directory (in loadSettings)
    setParameterFilePath(QCoreApplication::applicationDirPath());
    setDataFilePath(QCoreApplication::applicationDirPath());

    setWindowTitle(tr("Spin Diffusion (11 Dec 2020)"));

    calcInProgress=false;
    resize(QSize(400,600));

    sdCube = new TSDCube;

    createWidgets();
    createPanel();
    createConnections();

}

MainWindow::~MainWindow()
{
    // We save settings for the next run, before quitting the application
    saveSettings();
    delete sdCube;
}

void MainWindow::saveSettings()
{

 ;
}

void MainWindow::loadSettings()
{
;
}

void MainWindow::createWidgets()
{

    loadParametersPushButton = new QPushButton(tr("Load"));
    loadParametersPushButton->setFixedSize(100,50);
    saveParametersPushButton = new QPushButton(tr("Save"));
    saveParametersPushButton->setFixedSize(100,50);
    outputFileNamePushButton = new QPushButton(tr("..."));
    startStopPushButton = new QPushButton(tr("Start"));
    startStopPushButton->setFixedSize(80,80);

    currentStatusLabel = new QLabel;

    dCoeffLineEdit = new QLineEdit();
    dCoeffLineEdit->setAlignment(Qt::AlignRight);
    lengthLineEdit = new QLineEdit();
    lengthLineEdit->setAlignment(Qt::AlignRight);
    spinDensityLineEdit = new QLineEdit();
    spinDensityLineEdit->setAlignment(Qt::AlignRight);
    nOfDivisionsLineEdit = new QLineEdit();
    nOfDivisionsLineEdit->setAlignment(Qt::AlignRight);
    dtLineEdit = new QLineEdit();
    dtLineEdit->setAlignment(Qt::AlignRight);
    nObsLineEdit = new QLineEdit();
    nObsLineEdit->setAlignment(Qt::AlignRight);
    nDNPLineEdit = new QLineEdit();
    nDNPLineEdit->setAlignment(Qt::AlignRight);
    exchangeProbabilityLineEdit = new QLineEdit();
    exchangeProbabilityLineEdit->setAlignment(Qt::AlignRight);
    finalTimeLineEdit = new QLineEdit();
    finalTimeLineEdit->setAlignment(Qt::AlignRight);
    sourcePolarizationLineEdit = new QLineEdit();
    sourcePolarizationLineEdit->setAlignment(Qt::AlignRight);
    outputFileNameLineEdit = new QLineEdit();
    outputFileNameLineEdit->setAlignment(Qt::AlignLeft);
    plainTextEdit = new QPlainTextEdit;
    R1LineEdit = new QLineEdit();
    R1LineEdit->setAlignment(Qt::AlignRight);
    polarizationScalingCheckBox = new QCheckBox(tr("Scaling of spin-diffusion rate with polarization."));
	

}


void MainWindow::createPanel()
{
  QWidget *widget1 = new QWidget(this);
  setCentralWidget(widget1);

  // We make a vertical layout, which is a child of widget1
  QVBoxLayout *vLayout1 = new QVBoxLayout(widget1);

    QHBoxLayout *hLayout1 = new QHBoxLayout;
    QHBoxLayout *hLayout1b = new QHBoxLayout;
    QGridLayout *gLayout1 = new QGridLayout;
    QHBoxLayout *hLayout2 = new QHBoxLayout;

    hLayout1->addWidget(new QLabel(tr("Parameters")));
    hLayout1->addWidget(loadParametersPushButton);
    hLayout1->addWidget(saveParametersPushButton);
    hLayout1->addStretch();

    gLayout1->addWidget(new QLabel(tr("Diffusion Coefficient")),0,0,1,1);
    gLayout1->addWidget(dCoeffLineEdit,0,1,1,1);
    gLayout1->addWidget(new QLabel(tr("m^2/s")),0,2,1,1);

    gLayout1->addWidget(new QLabel(tr("# of spins per unit volume")),1,0,1,1);
    gLayout1->addWidget(spinDensityLineEdit,1,1,1,1);
    gLayout1->addWidget(new QLabel(tr("m^-3")),1,2,1,1);

    gLayout1->addWidget(new QLabel(tr("Length")),2,0,1,1);
    gLayout1->addWidget(lengthLineEdit,2,1,1,1);
    gLayout1->addWidget(new QLabel(tr("m")),2,2,1,1);

    gLayout1->addWidget(new QLabel(tr("# of divisions")),3,0,1,1);
    gLayout1->addWidget(nOfDivisionsLineEdit,3,1,1,1);
    //gLayout1->addWidget(new QLabel(tr("")),2,2,1,1);

    gLayout1->addWidget(new QLabel(tr("Delta t")),4,0,1,1);
    gLayout1->addWidget(dtLineEdit,4,1,1,1);
    gLayout1->addWidget(new QLabel(tr("s")),4,2,1,1);

    gLayout1->addWidget(new QLabel(tr("# of steps between e-n polarization transfer")),5,0,1,1);
    gLayout1->addWidget(nDNPLineEdit,5,1,1,1);

    gLayout1->addWidget(new QLabel(tr("exchange probability")),6,0,1,1);
    gLayout1->addWidget(exchangeProbabilityLineEdit,6,1,1,1);

    gLayout1->addWidget(new QLabel(tr("# of steps between records")),7,0,1,1);
    gLayout1->addWidget(nObsLineEdit,7,1,1,1);
//    gLayout1->addWidget(new QLabel(tr("s")),4,2,1,1);

    gLayout1->addWidget(new QLabel(tr("Final time")),8,0,1,1);
    gLayout1->addWidget(finalTimeLineEdit,8,1,1,1);
    gLayout1->addWidget(new QLabel(tr("s")),8,2,1,1);

    gLayout1->addWidget(new QLabel(tr("Source polarization")),9,0,1,1);
    gLayout1->addWidget(sourcePolarizationLineEdit,9,1,1,1);
//    gLayout1->addWidget(new QLabel(tr("")),6,2,1,1);


    gLayout1->addWidget(new QLabel(tr("Spin-lattice relaxation rate (1/T1)")),10,0,1,1);
    gLayout1->addWidget(R1LineEdit,10,1,1,1);
    gLayout1->addWidget(new QLabel(tr("1/s")),10,2,1,1);

    //    gLayout1->addWidget(new QLabel(tr("Scaling of spin-diffusion rate")),10,0,1,1);
    gLayout1->addWidget(polarizationScalingCheckBox,11,0,1,2);

    hLayout1b->addWidget(new QLabel(tr("Output file name")));
    hLayout1b->addWidget(outputFileNameLineEdit);
    hLayout1b->addWidget(outputFileNamePushButton);



    hLayout2->addWidget(startStopPushButton);
    hLayout2->addWidget(currentStatusLabel);

    vLayout1->addLayout(hLayout1);
    vLayout1->addLayout(gLayout1);
    vLayout1->addLayout(hLayout1b);
    vLayout1->addLayout(hLayout2);
    vLayout1->addWidget(plainTextEdit);

}

void MainWindow::createConnections()
{

    connect(saveParametersPushButton,SIGNAL(clicked()),this,SLOT(onSaveParametersPushButtonClicked()));
    connect(loadParametersPushButton,SIGNAL(clicked()),this,SLOT(onLoadParametersPushButtonClicked()));
    connect(outputFileNamePushButton,SIGNAL(clicked()),this,SLOT(onOutputFileNamePushButtonClicked()));
    connect(startStopPushButton,SIGNAL(clicked()),this,SLOT(onStartStopButtonClicked()));
    connect(sdCube,SIGNAL(sendMessage(QString)),this->currentStatusLabel,SLOT(setText(QString)));
    connect(sdCube,SIGNAL(sendData(QString)),this,SLOT(onReceiveData(QString)));
    connect(sdCube,SIGNAL(calcComplete()),this,SLOT(onCalcComplete()));
}

void MainWindow::onReceiveData(QString qs)
{
    plainTextEdit->appendPlainText(qs);
    plainTextEdit->moveCursor(QTextCursor::End);
}

void MainWindow::onLoadParametersPushButtonClicked()
{
    QString fileName=QFileDialog::getOpenFileName(this,tr("Open Parameter File"),
                                                  parameterFilePath());
    if(fileName.isNull()) {return;}

    QFileInfo fi(fileName);
    setParameterFilePath(fi.absolutePath());

    QSettings settings(fileName, QSettings::IniFormat);
    settings.beginGroup("SDCube");
      dCoeffLineEdit->setText(settings.value("Diffusion coefficient","").toString());
      lengthLineEdit->setText(settings.value("Length","").toString());
      spinDensityLineEdit->setText(settings.value("Number of spins per unit volume","").toString());
      nOfDivisionsLineEdit->setText(settings.value("Number of divisions","").toString());
      dtLineEdit->setText(settings.value("Delta t","").toString());
      nObsLineEdit->setText(settings.value("NObs","").toString());
      nDNPLineEdit->setText(settings.value("NDNP","").toString());
      exchangeProbabilityLineEdit->setText(settings.value("Exchange probability","").toString());
      finalTimeLineEdit->setText(settings.value("Final time","").toString());
      sourcePolarizationLineEdit->setText(settings.value("Source polarization","").toString());
      polarizationScalingCheckBox->setChecked(settings.value("Polarization scaling","").toBool());
      R1LineEdit->setText(settings.value("Spin-lattice relaxation rate","").toString());
      outputFileNameLineEdit->setText(settings.value("Output filename","").toString());
    settings.endGroup();
}


void MainWindow::onSaveParametersPushButtonClicked()
{
    if(!setupParams())
    {

      return;
    }

    QString fileName=QFileDialog::getSaveFileName(this,tr("Save Parameter File"),
                                                  parameterFilePath());
    if(fileName.isNull()) {return;}

    QFileInfo fi(fileName);
    setParameterFilePath(fi.absolutePath());

    QSettings settings(fileName, QSettings::IniFormat);
    settings.beginGroup("SDCube");
      settings.setValue("Diffusion coefficient", dCoeffLineEdit->text());
      settings.setValue("Length", lengthLineEdit->text());
      settings.setValue("Number of spins per unit volume", spinDensityLineEdit->text());
      settings.setValue("Number of divisions", nOfDivisionsLineEdit->text());
      settings.setValue("Delta t", dtLineEdit->text());
      settings.setValue("NDNP", nDNPLineEdit->text());
      settings.setValue("Exchange probability", exchangeProbabilityLineEdit->text());
      settings.setValue("NObs", nObsLineEdit->text());
      settings.setValue("Final time", finalTimeLineEdit->text());
      settings.setValue("Source polarization", sourcePolarizationLineEdit->text());
      settings.setValue("Polarization scaling", polarizationScalingCheckBox->isChecked());
      settings.setValue("Spin-lattice relaxation rate", R1LineEdit->text());
      settings.setValue("Output filename", outputFileNameLineEdit->text());
    settings.endGroup();
    settings.sync();


}

void MainWindow::onOutputFileNamePushButtonClicked()
{
    QString fileName = QFileDialog::getSaveFileName(this, tr("Output file"));
    if(fileName.isNull()) {return;}

    QFileInfo fi(fileName);
    setDataFilePath(fi.absolutePath());

    outputFileNameLineEdit->setText(fileName);

}

bool MainWindow::setupParams()
{
    QString qs;
    double d;
    int i;
    bool ok;

    //----- Diffusion coeff. -----
    qs=dCoeffLineEdit->text();
    d=qs.toDouble(&ok);
    if(!ok)
    {
        currentStatusLabel->setText("Invalid expression: " + qs);
        return false;
    }
    if(d<0)
    {
        currentStatusLabel->setText("Diffusion coefficient cannot be negative");
        return false;
    }

    sdCube->setDCoeff_xx(d);
    sdCube->setDCoeff_yy(d);
    sdCube->setDCoeff_zz(d);

    //----- Length -----
    qs=lengthLineEdit->text();
    d=qs.toDouble(&ok);
    if(!ok)
    {
        currentStatusLabel->setText("Invalid expression: " + qs);
        return false;
    }
    if(d<0)
    {
        currentStatusLabel->setText("Length cannot be negative: " + qs);
        return false;
    }
    sdCube->setLength(d);

    //----- # of spins per unit volume (spinDensity) -----
    qs=spinDensityLineEdit->text();
    d=qs.toDouble(&ok);
    if(!ok)
    {
        currentStatusLabel->setText("Invalid expression: " + qs);
        return false;
    }
    if(d<0)
    {
        currentStatusLabel->setText("Number of spins per unit volume cannot be negative: " + qs);
        return false;
    }
    sdCube->setSpinDensity(d);

    //----- # of divisions -----
    qs=nOfDivisionsLineEdit->text();
    i=qs.toInt(&ok);
    if(!ok)
    {
        currentStatusLabel->setText("Invalid expression: " + qs);
        return false;
    }
    if(i<0)
    {
        currentStatusLabel->setText("Number of divisions cannot be negative: " + qs);
        return false;
    }
    if(i<3)
    {
        currentStatusLabel->setText("Number of divisions is too small:" + qs);
    }
    sdCube->setNOfDivisions(i);

    //----- delta t -----
    qs=dtLineEdit->text();
    d=qs.toDouble(&ok);
    if(!ok)
    {
        currentStatusLabel->setText("Invalid expression: " + qs);
        return false;
    }
    if(d<0)
    {
        currentStatusLabel->setText("Dt cannot be negative: " + qs);
        return false;
    }
    sdCube->setDt(d);
    //----- # of steps between e-n polarization transfer -----
    qs=nDNPLineEdit->text();
    i=qs.toInt(&ok);
    if(!ok)
    {
        currentStatusLabel->setText("Invalid expression: " + qs);
        return false;
    }
    if(i<1)
    {
        currentStatusLabel->setText("Number of steps (between e-n polarization transfer) cannot be smaller than 1: " + qs);
        return false;
    }
    sdCube->setNDNP(i);

    //----- Exchange probability -----
    qs=exchangeProbabilityLineEdit->text();
    d=qs.toDouble(&ok);
    if(!ok)
    {
        currentStatusLabel->setText("Invalid expression: " + qs);
        return false;
    }
    if(d<0 || d>1)
    {
        currentStatusLabel->setText("Exchange probability must be between 0 and 1: " + qs);
        return false;
    }
    sdCube->setExchangeProbability(d);

    //----- # of steps between records -----
    qs=nObsLineEdit->text();
    i=qs.toInt(&ok);
    if(!ok)
    {
        currentStatusLabel->setText("Invalid expression: " + qs);
        return false;
    }
    if(i<1)
    {
        currentStatusLabel->setText("Number of steps (between records) cannot be smaller than 1: " + qs);
        return false;
    }
    sdCube->setNObs(i);

    //----- Final t -----
    qs=finalTimeLineEdit->text();
    d=qs.toDouble(&ok);
    if(!ok)
    {
        currentStatusLabel->setText("Invalid expression: " + qs);
        return false;
    }
    if(d<0)
    {
        currentStatusLabel->setText("Dt cannot be negative: " + qs);
        return false;
    }
    sdCube->setFinalTime(d);


    //----- spin-lattice relaxation rate -----
    qs=R1LineEdit->text();
    d=qs.toDouble(&ok);
    if(!ok)
    {
        currentStatusLabel->setText("Invalid expression: " + qs);
        return false;
    }
    if(d<0)
    {
        currentStatusLabel->setText("Spin-lattice relaxation rate cannot be negative: " + qs);
        return false;
    }
    sdCube->setR1(d);

    //----- source polarization -----
    qs=sourcePolarizationLineEdit->text();
    d=qs.toDouble(&ok);
    if(!ok)
    {
        currentStatusLabel->setText("Invalid expression: " + qs);
        return false;
    }
    if(d<0)
    {
        currentStatusLabel->setText("Polarization cannot be negative: " + qs);
        return false;
    }
    if(d>1)
    {
        currentStatusLabel->setText("Polarization cannot exceed 1: " + qs);
        return false;
    }
    sdCube->setSourcePolarization(d);

    //----- scaling of diffusion rate with polarization -----
    sdCube->setPolarizationScaling(polarizationScalingCheckBox->isChecked());

    //----- Output filename
    sdCube->setOutputFileName(outputFileNameLineEdit->text());

    //----------------------
    if(!sdCube->setupParams())
    {
        currentStatusLabel->setText(sdCube->errorMessage());
        return false;
    }

    currentStatusLabel->setText("No error found.");
    return true;
}

void MainWindow::onStartStopButtonClicked()
{

    if(!calcInProgress)
    {
        if(!setupParams()) return;
        plainTextEdit->clear();
        startStopPushButton->setText("Stop");
        calcInProgress=true;
        sdCube->doCalc();
    }
    else
    {
        sdCube->stop();
        startStopPushButton->setText("Start");
        calcInProgress=false;
    }
}

void MainWindow::onCalcComplete()
{
    startStopPushButton->setText("Start");
    calcInProgress=false;
}
