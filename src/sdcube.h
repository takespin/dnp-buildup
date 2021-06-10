#ifndef TSDCUBE_H
#define TSDCUBE_H

#include <QObject>
#include <QThread>
#include <QMutex>
#include <QWaitCondition>
#include <QVector>

class TSDCube : public QThread
{
    Q_OBJECT
public:
//    enum simMode {Continuous, Transient};
    TSDCube();
    ~TSDCube()
    {
      mutex.lock();
        stopped=true;
        condition.wakeOne();
      mutex.unlock();
      wait();


    }

    QVector< QVector< QVector<double>>> P, Q;
    bool setupParams();

    void doCalc();
    void stop() {QMutexLocker locker(&mutex); stopped=true; condition.wakeAll();}

  //  int divisions() {return FDivisions;}
  //  void setDivisions(int k);

    double lambda_xx() {return FLambda_xx;}
    double lambda_yy() {return FLambda_yy;}
    double lambda_zz() {return FLambda_zz;}

    double eta() {return FEta;}

    void setPolarizationScaling(bool b) {FPolarizationScaling=b;}
    bool polarizationScaling() {return FPolarizationScaling;}

    double length() {return FLength;}
    void setLength(double f) {FLength=f;}
    int NOfDivisions() {return FNOfDivisions;}
    void setNOfDivisions(int i);

    void setDt(double f) {FDt=f;}
    double dt() {return FDt;}
    void setNObs(int i) {FNObs=i;}
    int nObs() {return FNObs;}
    void setNDNP(int i) {FNDNP=i;}
    int nDNP() {return FNDNP;}
    void setExchangeProbability(double d) {FExchangeProbability=d;}
    double exchangeProbability() {return FExchangeProbability;}
    void setFinalTime(double f) {FFinalTime=f;}
    double finalTime() {return FFinalTime;}
    void setDx(double f) {FDx=f;}
    double dx() {return FDx;}

    double spinDensity() {return FSpinDensity;}
    void setSpinDensity(double f) {FSpinDensity=f;}

    double DCoeff_xx() {return FDCoeff_xx;}
    void setDCoeff_xx(double f) {FDCoeff_xx=f;}
    double DCoeff_yy() {return FDCoeff_yy;}
    void setDCoeff_yy(double f) {FDCoeff_yy=f;}
    double DCoeff_zz() {return FDCoeff_zz;}
    void setDCoeff_zz(double f) {FDCoeff_zz=f;}
    double sourcePolarization() {return FSourcePolarization;}
    void setSourcePolarization(double f) {FSourcePolarization=f;}

    double R1() {return FR1;}
    void setR1(double d) {FR1=d;}

    void setOutputFileName(QString qs) {FOutputFileName=qs;}
    QString outputFileName() {return FOutputFileName;}

    double averageP();
    QString errorMessage() {return FErrorMessage;}
    void setErrorMessage(QString qs) {FErrorMessage=qs;}
    void evolve();
    void evolve2();
    void Q2P();
    void P2Q();


    bool evolve_P();
    void evolve_NP();

    void spinLatticeRelaxation();

    void prepareCalc();



signals:
    void sendMessage(QString qs);
    void sendData(QString qs);
    void calcComplete();
    void sendCurrentTime(double f);

protected:
    void run();

private:
    QString FErrorMessage;
    bool FError;
    double FFinalTime;
   // int FDivisions;
    QMutex mutex;
    QWaitCondition condition;
    bool volatile stopped;
    double FLength;
    int FNOfDivisions;
    double FLambda_xx, FLambda_yy, FLambda_zz;
    double FEta;
    double FDx, FDt;
    int FNDNP;
    double FExchangeProbability;
    int FNObs;
    double FDCoeff_xx,FDCoeff_yy,FDCoeff_zz;
    QString FOutputFileName;
    double FSourcePolarization;
    double FSpinDensity;
    double FR1;

    bool FPolarizationScaling;

    enum cState {K1,L1,L2,L3,L4,E,M1,M2,M3,ERR};
    enum cState currentState;
    enum cState nextState;

    double nVoxel; // # of spin per voxel.


};

#endif // TSDCUBE_H
