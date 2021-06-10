#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QLineEdit>
#include <QLabel>
#include <QString>
#include <QPushButton>
#include <QCheckBox>
#include <QPlainTextEdit>

#include "sdcube.h"


class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

    TSDCube *sdCube;

    QString parameterFilePath() {return FParameterFilePath;}
    QString dataFilePath() {return FDataFilePath;}

public slots:
    void loadSettings();
    void saveSettings();


private slots:

    void onLoadParametersPushButtonClicked();
    void onSaveParametersPushButtonClicked();
    void onOutputFileNamePushButtonClicked();
    void onStartStopButtonClicked();
    void onCalcComplete();
    void setParameterFilePath(QString qs) {FParameterFilePath=qs;}
    void setDataFilePath(QString qs) {FDataFilePath=qs;}
    void onReceiveData(QString qs);

private:
    void createWidgets();
    void createPanel();
    void createConnections();
    void loadParameters();
    void saveParameters();
    bool setupParams();
    QString FParameterFilePath;
    QString FDataFilePath;
    bool calcInProgress;

    QPushButton *loadParametersPushButton;
    QPushButton *saveParametersPushButton;
    QPushButton *outputFileNamePushButton;
    QPushButton *startStopPushButton;

    QLabel *currentStatusLabel;

    QLineEdit *dCoeffLineEdit;
    QLineEdit *lengthLineEdit;
    QLineEdit *spinDensityLineEdit;
    QLineEdit *nOfDivisionsLineEdit;
    QLineEdit *dtLineEdit;
    QLineEdit *nObsLineEdit;
    QLineEdit *nDNPLineEdit;
    QLineEdit *exchangeProbabilityLineEdit;
    QLineEdit *finalTimeLineEdit;
    QLineEdit *sourcePolarizationLineEdit;
    QLineEdit *outputFileNameLineEdit;
    QCheckBox *polarizationScalingCheckBox;
    QPlainTextEdit *plainTextEdit;
    QLineEdit *R1LineEdit;


};
#endif // MAINWINDOW_H
