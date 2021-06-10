#include <QDebug>
#include <QFile>
#include "sdcube.h"
#include "math.h"
#include "cstdlib"
#include <QDebug>

TSDCube::TSDCube()
{
    FPolarizationScaling=false;
    setR1(0);
}


void TSDCube::doCalc()
{
    QMutexLocker locker(&mutex);
    stopped=false;
    if(!isRunning())
    {
      prepareCalc();
      start(HighPriority);
    }
    else
    {
      condition.wakeOne();
    }
}

void TSDCube::prepareCalc()
{

    for(int i=0; i<NOfDivisions(); i++) {
    for(int j=0; j<NOfDivisions(); j++) {
    for(int k=0; k<NOfDivisions(); k++) {
      P[i][j][k]=0;
      Q[i][j][k]=0;
    }
    }
    }

    nVoxel=spinDensity()*(length()/NOfDivisions())
                        *(length()/NOfDivisions())
                        *(length()/NOfDivisions());

    qDebug() << "nVoxel: " << nVoxel;

    return;
}

void TSDCube::setNOfDivisions(int k)
{
    if(k>2) FNOfDivisions=k; else FNOfDivisions=3;
//    if(k>1000) FDivisions=1000;

    P.resize(FNOfDivisions);
    Q.resize(FNOfDivisions);
    for(int i=0; i<FNOfDivisions; i++)
    {
      P[i].resize(FNOfDivisions);
      Q[i].resize(FNOfDivisions);
      for(int j=0; j<FNOfDivisions; j++)
      {
          P[i][j].resize(FNOfDivisions);
          Q[i][j].resize(FNOfDivisions);
      }
    }

}



bool TSDCube::setupParams()
{
    if (NOfDivisions()<1)
    {
       FError=true;
       FErrorMessage= tr("Number of Division (") + QString::number(NOfDivisions()) + tr(") is invalid.");
       emit sendMessage(errorMessage());
       return false;
    }
    FDx=length()/NOfDivisions();
    if(dx()<=0)
    {
        FError=true;
        FErrorMessage= tr("dx (") + QString::number(dx()) + tr(") is invalid.");
        emit sendMessage(errorMessage());
        return false;
    }
    FLambda_xx = DCoeff_xx()*dt()/(dx()*dx());
    FLambda_yy = DCoeff_yy()*dt()/(dx()*dx());
    FLambda_zz = DCoeff_zz()*dt()/(dx()*dx());
    FEta = 1.0 - (2.0*lambda_xx()) - (2.0*lambda_yy()) - (2.0*lambda_zz());
    return true;
}

double TSDCube::averageP()
{
  double aP=0;
  int n=NOfDivisions();
  int nCubed=n*n*n;

  for(int i=0; i<n; i++) { for(int j=0; j<n; j++) { for(int k=0; k<n; k++) {
     aP += P.at(i).at(j).at(k);
  }}}

  aP = aP/nCubed;

  return aP;
}

// Spin diffusion without polarization scaling
void TSDCube::evolve_NP()
{
  int n=NOfDivisions();
  double diff;

  int ip,jp,kp;

  P2Q(); // copy P -> Q

  for(int k=0; k<n; k++) {
  for(int i=0; i<n; i++) {
  for(int j=0; j<n; j++) {
    ip=i+1; jp=j+1; kp=k+1;

    // We set periodic boundaries
    if(ip>n-1) {ip=0;}
    if(jp>n-1) {jp=0;}
    if(kp>n-1) {kp=0;}

    // flow along x
    diff=lambda_xx()*(P[ip][j][k]-P[i][j][k]);
    Q[ip][j][k] -= diff;
    Q[i][j][k]  += diff;

    // flow along +y
    diff=lambda_yy()*(P[i][jp][k]-P[i][j][k]);
    Q[i][jp][k] -= diff;
    Q[i][j][k]  += diff;

    // flow along +z
    diff=lambda_zz()*(P[i][j][kp]-P[i][j][k]);
    Q[i][j][kp] -= diff;
    Q[i][j][k]  += diff;

  } // j
  } // i
  } // k

  Q2P(); // copy Q -> P

}



// Polarization-scaled spin diffusion
bool TSDCube::evolve_P()
{
  int n=NOfDivisions();
  double fp;
//  double fm;
  double diff;

  int ip,jp,kp;

  P2Q(); // copy P -> Q

  for(int k=0; k<n; k++) {
  for(int i=0; i<n; i++) {
  for(int j=0; j<n; j++) {
    ip=i+1; jp=j+1; kp=k+1;

    // We set periodic boundaries
    if(ip>n-1) {ip=0;}
    if(jp>n-1) {jp=0;}
    if(kp>n-1) {kp=0;}

    // flow along x
    fp=sqrt( 1.0 - pow(0.5*(P[ip][j][k]+P[i][j][k]), 2) );
    if(fp==0)
    {
        FError=true;
        setErrorMessage("Error: division by 0 is not allowed.");
        return false;
    }
    diff=lambda_xx()*(P[ip][j][k]-P[i][j][k])/fp;
    Q[ip][j][k] -= diff;
    Q[i][j][k]  += diff;

    // flow along +y
    fp=sqrt( 1.0 - pow(0.5*(P[i][jp][k]+P[i][j][k]), 2) );
    if(fp==0)
    {
        FError=true;
        setErrorMessage("Error: division by 0 is not allowed.");
        return false;
    }
    diff=lambda_yy()*(P[i][jp][k]-P[i][j][k])/fp;
    Q[i][jp][k] -= diff;
    Q[i][j][k]  += diff;

    // flow along +z
    fp=sqrt( 1.0 - pow(0.5*(P[i][j][kp]+P[i][j][k]), 2) );
    if(fp==0)
    {
        FError=true;
        setErrorMessage("Error: division by 0 is not allowed.");
        return false;
    }
    diff=lambda_zz()*(P[i][j][kp]-P[i][j][k])/fp;
    Q[i][j][kp] -= diff;
    Q[i][j][k]  += diff;

  } // j
  } // i
  } // k

  Q2P(); // copy Q -> P

  return true;
}


void TSDCube::spinLatticeRelaxation()
{
    // spin-lattice relaxation (P)
    for(int i=0; i<NOfDivisions(); i++) {
    for(int j=0; j<NOfDivisions(); j++) {
    for(int k=0; k<NOfDivisions(); k++) {
      P[i][j][k] -= R1()*dt()*P[i][j][k];
      // Here we assume the thermal polarization is negligibly small.
    }
    }
    }

}


void TSDCube::evolve()
{
    int n=NOfDivisions();
    //+=+=+=+ CASE A: (k=0) +=+=+=+
    // CASE A1: (i=0, j=0, k=0)
      Q[0][0][0] = eta()*P[0][0][0]
                 + lambda_xx() * ( P[n-1][0][0] + P[1][0][0] )
                 + lambda_yy() * ( P[0][n-1][0] + P[0][1][0] )
                 + lambda_zz() * ( P[0][0][n-1] + P[0][0][1] );
      // CASE A2: (i=n-1, j=0, k=0)
        Q[n-1][0][0] = eta()*P[n-1][0][0]
                     + lambda_xx() * ( P[n-2][0][0]   + P[0][0][0] )
                     + lambda_yy() * ( P[n-1][n-1][0] + P[n-1][1][0] )
                     + lambda_zz() * ( P[n-1][0][n-1] + P[n-1][0][1] );
      // CASE A3: (i=n-1, j=n-1, k=0)
        Q[n-1][n-1][0] = eta()*P[n-1][n-1][0]
                       + lambda_xx() * ( P[n-2][n-1][0]   + P[0][n-1][0] )
                       + lambda_yy() * ( P[n-1][n-2][0]   + P[n-1][0][0] )
                       + lambda_zz() * ( P[n-1][n-1][n-1] + P[n-1][n-1][1] );
      // CASE A4: (i=0, j=n-1, k=0)
        Q[0][n-1][0] = eta()*P[0][n-1][0]
                     + lambda_xx() * ( P[n-1][n-1][0] + P[1][n-1][0] )
                     + lambda_yy() * ( P[0][n-2][0]   + P[0][0][0] )
                     + lambda_zz() * ( P[0][n-1][n-1] + P[0][n-1][1] );
      // CASE A5: (0<i<n-1, j=0, k=0)
      for(int i=1; i<n-1; i++) {
        Q[i][0][0] = eta()*P[i][0][0]
                   + lambda_xx() * ( P[i-1][0][0] + P[i+1][0][0] )
                   + lambda_yy() * ( P[i][n-1][0] + P[i][1][0] )
                   + lambda_zz() * ( P[i][0][n-1] + P[i][0][1] );
      }// i
      // CASE A6: (i=n-1, 0<j<n-1, k=0)
      for(int j=1; j<n-1; j++) {
        Q[n-1][j][0] = eta()*P[n-1][j][0]
                     + lambda_xx() * ( P[n-2][j][0]   + P[0][j][0] )
                     + lambda_yy() * ( P[n-1][j-1][0] + P[n-1][j+1][0] )
                     + lambda_zz() * ( P[n-1][j][n-1] + P[n-1][j][1] );
      } // j
      // CASE A7: (0<i<n-1, j=n-1, k=0)
      for(int i=1; i<n-1; i++) {
        Q[i][n-1][0] = eta()*P[i][n-1][0]
                     + lambda_xx() * ( P[i-1][n-1][0] + P[i+1][n-1][0] )
                     + lambda_yy() * ( P[i][n-2][0]   + P[i][0][0] )
                     + lambda_zz() * ( P[i][n-1][n-1] + P[i][n-1][1] );
      } // i
      // CASE A8: (i=0, 0<j<n-1, k=0)
      for(int j=1; j<n-1; j++) {
        Q[0][j][0] = eta()*P[0][j][0]
                   + lambda_xx() * ( P[n-1][j][0] + P[1][j][0] )
                   + lambda_yy() * ( P[0][j-1][0] + P[0][j+1][0] )
                   + lambda_zz() * ( P[0][j][n-1] + P[0][j][1] );
      } // j
      // CASE A9: (0<i<n-1, 0<j<n-1, k=0)
      for(int i=1; i<n-1; i++) {
      for(int j=1; j<n-1; j++) {
        Q[i][j][0] = eta()*P[i][j][0]
                   + lambda_xx() * ( P[i-1][j][0] + P[i+1][j][0] )
                   + lambda_yy() * ( P[i][j-1][0] + P[i][j+1][0] )
                   + lambda_zz() * ( P[i][j][n-1] + P[i][j][1] );
      }// j
      }// i


      //
      //+=+=+=+ CASE B: (0<k<n-1) +=+=+=+
      // CASE B1: (i=0, j=0, 0<k<n-1)
      for(int k=1; k<n-1; k++) {
        Q[0][0][k] = eta()*P[0][0][k]
                   + lambda_xx() * ( P[n-1][0][k] + P[1][0][k] )
                   + lambda_yy() * ( P[0][n-1][k] + P[0][1][k] )
                   + lambda_zz() * ( P[0][0][k-1] + P[0][0][k+1] );
      } // k
      // CASE B2: (i=n-1, j=0, 0<k<n-1)
      for(int k=1; k<n-1; k++) {
        Q[n-1][0][k] = eta()*P[n-1][0][k]
                     + lambda_xx() * ( P[n-2][0][k]   + P[0][0][k] )
                     + lambda_yy() * ( P[n-1][n-1][k] + P[n-1][1][k] )
                     + lambda_zz() * ( P[n-1][0][k-1] + P[n-1][0][k+1] );
      } // k
      // CASE B3: (i=n-1, j=n-1, 0<k<n-1)
      for(int k=1; k<n-1; k++) {
        Q[n-1][n-1][k] = eta()*P[n-1][n-1][k]
                       + lambda_xx() * ( P[n-2][n-1][k]   + P[0][n-1][k] )
                       + lambda_yy() * ( P[n-1][n-2][k]   + P[n-1][0][k] )
                       + lambda_zz() * ( P[n-1][n-1][k-1] + P[n-1][n-1][k+1] );
      } // k
      // CASE B4: (i=0, j=n-1, 0<k<n-1)
      for(int k=1; k<n-1; k++) {
        Q[0][n-1][k] = eta()*P[0][n-1][k]
                     + lambda_xx() * ( P[n-1][n-1][k] + P[1][n-1][k] )
                     + lambda_yy() * ( P[0][n-2][k]   + P[0][0][k] )
                     + lambda_zz() * ( P[0][n-1][k-1] + P[0][n-1][k+1] );
      } // k
      // CASE B5: (0<i<n-1, j=0, 0<k<n-1)
      for(int k=1; k<n-1; k++) {
      for(int i=1; i<n-1; i++) {
        Q[i][0][k] = eta()*P[i][0][k]
                   + lambda_xx() * ( P[i-1][0][k] + P[i+1][0][k] )
                   + lambda_yy() * ( P[i][n-1][k] + P[i][1][k] )
                   + lambda_zz() * ( P[i][0][k-1] + P[i][0][k+1] );
      } // i
      } // k
      // CASE B6: (i=n-1, 0<j<n-1, 0<k<n-1)
      for(int k=1; k<n-1; k++) {
      for(int j=1; j<n-1; j++) {
        Q[n-1][j][k] = eta()*P[n-1][j][k]
                     + lambda_xx() * ( P[n-2][j][k]   + P[0][j][k] )
                     + lambda_yy() * ( P[n-1][j-1][k] + P[n-1][j+1][k] )
                     + lambda_zz() * ( P[n-1][j][k-1] + P[n-1][j][k+1] );
      } // j
      } // k
      // CASE B7: (0<i<n-1, j=n-1, 0<k<n-1)
      for(int k=1; k<n-1; k++) {
      for(int i=1; i<n-1; i++) {
        Q[i][n-1][k] = eta()*P[i][n-1][k]
                     + lambda_xx() * ( P[i-1][n-1][k] + P[i+1][n-1][k] )
                     + lambda_yy() * ( P[i][n-2][k]   + P[i][0][k] )
                     + lambda_zz() * ( P[i][n-1][k-1] + P[i][n-1][k+1] );
      } // i
      } // k
      // CASE B8: (i=0, 0<j<n-1, 0<k<n-1)
      for(int k=1; k<n-1; k++) {
      for(int j=1; j<n-1; j++) {
        Q[0][j][k] = eta()*P[0][j][k]
                   + lambda_xx() * ( P[n-1][j][k] + P[1][j][k] )
                   + lambda_yy() * ( P[0][j-1][k] + P[0][j+1][k] )
                   + lambda_zz() * ( P[0][j][k-1] + P[0][j][k+1] );
      } // j
      } // k
      // CASE B9: (0<i<n-1, 0<j<n-1, 0<k<n-1)
      for(int k=1; k<n-1; k++) {
      for(int i=1; i<n-1; i++) {
      for(int j=1; j<n-1; j++) {
        Q[i][j][k] = eta()*P[i][j][k]
                   + lambda_xx() * ( P[i-1][j][k] + P[i+1][j][k] )
                   + lambda_yy() * ( P[i][j-1][k] + P[i][j+1][k] )
                   + lambda_zz() * ( P[i][j][k-1] + P[i][j][k+1] );
      } // j
      } // i
      } // k

      //+=+=+=+ CASE C: (k=n-1) +=+=+=+
      // CASE C1: (i=0, j=0, k=n-1)
        Q[0][0][n-1] = eta()*P[0][0][n-1]
                     + lambda_xx() * ( P[n-1][0][n-1] + P[1][0][n-1] )
                     + lambda_yy() * ( P[0][n-1][n-1] + P[0][1][n-1] )
                     + lambda_zz() * ( P[0][0][n-2]   + P[0][0][0] );
      // CASE C2: (i=n-1, j=0, k=n-1)
        Q[n-1][0][n-1] = eta()*P[n-1][0][n-1]
                       + lambda_xx() * ( P[n-2][0][n-1]   + P[0][0][n-1] )
                       + lambda_yy() * ( P[n-1][n-1][n-1] + P[n-1][1][n-1] )
                       + lambda_zz() * ( P[n-1][0][n-2]   + P[n-1][0][0] );
      // CASE C3: (i=n-1, j=n-1, k=n-1)
        Q[n-1][n-1][n-1] = eta()*P[n-1][n-1][n-1]
                         + lambda_xx() * ( P[n-2][n-1][n-1] + P[0][n-1][n-1] )
                         + lambda_yy() * ( P[n-1][n-2][n-1] + P[n-1][0][n-1] )
                         + lambda_zz() * ( P[n-1][n-1][n-2] + P[n-1][n-1][0] );
      // CASE C4: (i=0, j=n-1, k=n-1)
        Q[0][n-1][n-1] = eta()*P[0][n-1][n-1]
                       + lambda_xx() * ( P[n-1][n-1][n-1] + P[1][n-1][n-1] )
                       + lambda_yy() * ( P[0][n-2][n-1]   + P[0][0][n-1] )
                       + lambda_zz() * ( P[0][n-1][n-2]   + P[0][n-1][0] );
      // CASE C5: (0<i<n-1, j=0, k=n-1)
      for(int i=1; i<n-1; i++) {
        Q[i][0][n-1] = eta()*P[i][0][n-1]
                     + lambda_xx() * ( P[i-1][0][n-1] + P[i+1][0][n-1] )
                     + lambda_yy() * ( P[i][n-1][n-1] + P[i][1][n-1] )
                     + lambda_zz() * ( P[i][0][n-2]   + P[i][0][0] );
      } // i
      // CASE C6: (i=n-1, 0<j<n-1, k=n-1)
      for(int j=1; j<n-1; j++) {
        Q[n-1][j][n-1] = eta()*P[n-1][j][n-1]
                       + lambda_xx() * ( P[n-2][j][n-1]   + P[0][j][n-1] )
                       + lambda_yy() * ( P[n-1][j-1][n-1] + P[n-1][j+1][n-1] )
                       + lambda_zz() * ( P[n-1][j][n-2]   + P[n-1][j][0] );
      } // j
      // CASE C7: (0<i<n-1, j=n-1, k=n-1)
      for(int i=1; i<n-1; i++) {
        Q[i][n-1][n-1] = eta()*P[i][n-1][n-1]
                       + lambda_xx() * ( P[i-1][n-1][n-1] + P[i+1][n-1][n-1] )
                       + lambda_yy() * ( P[i][n-2][n-1]   + P[i][0][n-1] )
                       + lambda_zz() * ( P[i][n-1][n-2]   + P[i][n-1][0] );
      } // i
      // CASE C8: (i=0, 0<j<n-1, k=n-1)
      for(int j=1; j<n-1; j++) {
        Q[0][j][n-1] = eta()*P[0][j][n-1]
                     + lambda_xx() * ( P[n-1][j][n-1] + P[1][j][n-1] )
                     + lambda_yy() * ( P[0][j-1][n-1] + P[0][j+1][n-1] )
                     + lambda_zz() * ( P[0][j][n-2]   + P[0][j][0] );
      } // j
      // CASE C9: (0<i<n-1, 0<j<n-1, k=n-1)
      for(int i=1; i<n-1; i++) {
      for(int j=1; j<n-1; j++) {
        Q[i][j][n-1] = eta()*P[i][j][n-1]
                     + lambda_xx() * ( P[i-1][j][n-1] + P[i+1][j][n-1] )
                     + lambda_yy() * ( P[i][j-1][n-1] + P[i][j+1][n-1] )
                     + lambda_zz() * ( P[i][j][n-2]   + P[i][j][0] );
      } // j
      } // i

    return;
}

void TSDCube::evolve2()
{
    int n=NOfDivisions();
    //+=+=+=+ CASE A: (k=0) +=+=+=+
    // CASE A1: (i=0, j=0, k=0)
      P[0][0][0] = eta()*Q[0][0][0]
                 + lambda_xx() * ( Q[n-1][0][0] + Q[1][0][0] )
                 + lambda_yy() * ( Q[0][n-1][0] + Q[0][1][0] )
                 + lambda_zz() * ( Q[0][0][n-1] + Q[0][0][1] );
      // CASE A2: (i=n-1, j=0, k=0)
        P[n-1][0][0] = eta()*Q[n-1][0][0]
                     + lambda_xx() * ( Q[n-2][0][0]   + Q[0][0][0] )
                     + lambda_yy() * ( Q[n-1][n-1][0] + Q[n-1][1][0] )
                     + lambda_zz() * ( Q[n-1][0][n-1] + Q[n-1][0][1] );
      // CASE A3: (i=n-1, j=n-1, k=0)
        P[n-1][n-1][0] = eta()*Q[n-1][n-1][0]
                       + lambda_xx() * ( Q[n-2][n-1][0]   + Q[0][n-1][0] )
                       + lambda_yy() * ( Q[n-1][n-2][0]   + Q[n-1][0][0] )
                       + lambda_zz() * ( Q[n-1][n-1][n-1] + Q[n-1][n-1][1] );
      // CASE A4: (i=0, j=n-1, k=0)
        P[0][n-1][0] = eta()*Q[0][n-1][0]
                     + lambda_xx() * ( Q[n-1][n-1][0] + Q[1][n-1][0] )
                     + lambda_yy() * ( Q[0][n-2][0]   + Q[0][0][0] )
                     + lambda_zz() * ( Q[0][n-1][n-1] + Q[0][n-1][1] );
      // CASE A5: (0<i<n-1, j=0, k=0)
      for(int i=1; i<n-1; i++) {
        P[i][0][0] = eta()*Q[i][0][0]
                   + lambda_xx() * ( Q[i-1][0][0] + Q[i+1][0][0] )
                   + lambda_yy() * ( Q[i][n-1][0] + Q[i][1][0] )
                   + lambda_zz() * ( Q[i][0][n-1] + Q[i][0][1] );
      }// i
      // CASE A6: (i=n-1, 0<j<n-1, k=0)
      for(int j=1; j<n-1; j++) {
        P[n-1][j][0] = eta()*Q[n-1][j][0]
                     + lambda_xx() * ( Q[n-2][j][0]   + Q[0][j][0] )
                     + lambda_yy() * ( Q[n-1][j-1][0] + Q[n-1][j+1][0] )
                     + lambda_zz() * ( Q[n-1][j][n-1] + Q[n-1][j][1] );
      } // j
      // CASE A7: (0<i<n-1, j=n-1, k=0)
      for(int i=1; i<n-1; i++) {
        P[i][n-1][0] = eta()*Q[i][n-1][0]
                     + lambda_xx() * ( Q[i-1][n-1][0] + Q[i+1][n-1][0] )
                     + lambda_yy() * ( Q[i][n-2][0]   + Q[i][0][0] )
                     + lambda_zz() * ( Q[i][n-1][n-1] + Q[i][n-1][1] );
      } // i
      // CASE A8: (i=0, 0<j<n-1, k=0)
      for(int j=1; j<n-1; j++) {
        P[0][j][0] = eta()*Q[0][j][0]
                   + lambda_xx() * ( Q[n-1][j][0] + Q[1][j][0] )
                   + lambda_yy() * ( Q[0][j-1][0] + Q[0][j+1][0] )
                   + lambda_zz() * ( Q[0][j][n-1] + Q[0][j][1] );
      } // j
      // CASE A9: (0<i<n-1, 0<j<n-1, k=0)
      for(int i=1; i<n-1; i++) {
      for(int j=1; j<n-1; j++) {
        P[i][j][0] = eta()*Q[i][j][0]
                   + lambda_xx() * ( Q[i-1][j][0] + Q[i+1][j][0] )
                   + lambda_yy() * ( Q[i][j-1][0] + Q[i][j+1][0] )
                   + lambda_zz() * ( Q[i][j][n-1] + Q[i][j][1] );
      }// j
      }// i


      //
      //+=+=+=+ CASE B: (0<k<n-1) +=+=+=+
      // CASE B1: (i=0, j=0, 0<k<n-1)
      for(int k=1; k<n-1; k++) {
        P[0][0][k] = eta()*Q[0][0][k]
                   + lambda_xx() * ( Q[n-1][0][k] + Q[1][0][k] )
                   + lambda_yy() * ( Q[0][n-1][k] + Q[0][1][k] )
                   + lambda_zz() * ( Q[0][0][k-1] + Q[0][0][k+1] );
      } // k
      // CASE B2: (i=n-1, j=0, 0<k<n-1)
      for(int k=1; k<n-1; k++) {
        P[n-1][0][k] = eta()*Q[n-1][0][k]
                     + lambda_xx() * ( Q[n-2][0][k]   + Q[0][0][k] )
                     + lambda_yy() * ( Q[n-1][n-1][k] + Q[n-1][1][k] )
                     + lambda_zz() * ( Q[n-1][0][k-1] + Q[n-1][0][k+1] );
      } // k
      // CASE B3: (i=n-1, j=n-1, 0<k<n-1)
      for(int k=1; k<n-1; k++) {
        P[n-1][n-1][k] = eta()*Q[n-1][n-1][k]
                       + lambda_xx() * ( Q[n-2][n-1][k]   + Q[0][n-1][k] )
                       + lambda_yy() * ( Q[n-1][n-2][k]   + Q[n-1][0][k] )
                       + lambda_zz() * ( Q[n-1][n-1][k-1] + Q[n-1][n-1][k+1] );
      } // k
      // CASE B4: (i=0, j=n-1, 0<k<n-1)
      for(int k=1; k<n-1; k++) {
        P[0][n-1][k] = eta()*Q[0][n-1][k]
                     + lambda_xx() * ( Q[n-1][n-1][k] + Q[1][n-1][k] )
                     + lambda_yy() * ( Q[0][n-2][k]   + Q[0][0][k] )
                     + lambda_zz() * ( Q[0][n-1][k-1] + Q[0][n-1][k+1] );
      } // k
      // CASE B5: (0<i<n-1, j=0, 0<k<n-1)
      for(int k=1; k<n-1; k++) {
      for(int i=1; i<n-1; i++) {
        P[i][0][k] = eta()*Q[i][0][k]
                   + lambda_xx() * ( Q[i-1][0][k] + Q[i+1][0][k] )
                   + lambda_yy() * ( Q[i][n-1][k] + Q[i][1][k] )
                   + lambda_zz() * ( Q[i][0][k-1] + Q[i][0][k+1] );
      } // i
      } // k
      // CASE B6: (i=n-1, 0<j<n-1, 0<k<n-1)
      for(int k=1; k<n-1; k++) {
      for(int j=1; j<n-1; j++) {
        P[n-1][j][k] = eta()*Q[n-1][j][k]
                     + lambda_xx() * ( Q[n-2][j][k]   + Q[0][j][k] )
                     + lambda_yy() * ( Q[n-1][j-1][k] + Q[n-1][j+1][k] )
                     + lambda_zz() * ( Q[n-1][j][k-1] + Q[n-1][j][k+1] );
      } // j
      } // k
      // CASE B7: (0<i<n-1, j=n-1, 0<k<n-1)
      for(int k=1; k<n-1; k++) {
      for(int i=1; i<n-1; i++) {
        P[i][n-1][k] = eta()*Q[i][n-1][k]
                     + lambda_xx() * ( Q[i-1][n-1][k] + Q[i+1][n-1][k] )
                     + lambda_yy() * ( Q[i][n-2][k]   + Q[i][0][k] )
                     + lambda_zz() * ( Q[i][n-1][k-1] + Q[i][n-1][k+1] );
      } // i
      } // k
      // CASE B8: (i=0, 0<j<n-1, 0<k<n-1)
      for(int k=1; k<n-1; k++) {
      for(int j=1; j<n-1; j++) {
        P[0][j][k] = eta()*Q[0][j][k]
                   + lambda_xx() * ( Q[n-1][j][k] + Q[1][j][k] )
                   + lambda_yy() * ( Q[0][j-1][k] + Q[0][j+1][k] )
                   + lambda_zz() * ( Q[0][j][k-1] + Q[0][j][k+1] );
      } // j
      } // k
      // CASE B9: (0<i<n-1, 0<j<n-1, 0<k<n-1)
      for(int k=1; k<n-1; k++) {
      for(int i=1; i<n-1; i++) {
      for(int j=1; j<n-1; j++) {
        P[i][j][k] = eta()*Q[i][j][k]
                   + lambda_xx() * ( Q[i-1][j][k] + Q[i+1][j][k] )
                   + lambda_yy() * ( Q[i][j-1][k] + Q[i][j+1][k] )
                   + lambda_zz() * ( Q[i][j][k-1] + Q[i][j][k+1] );
      } // j
      } // i
      } // k

      //+=+=+=+ CASE C: (k=n-1) +=+=+=+
      // CASE C1: (i=0, j=0, k=n-1)
        P[0][0][n-1] = eta()*Q[0][0][n-1]
                     + lambda_xx() * ( Q[n-1][0][n-1] + Q[1][0][n-1] )
                     + lambda_yy() * ( Q[0][n-1][n-1] + Q[0][1][n-1] )
                     + lambda_zz() * ( Q[0][0][n-2]   + Q[0][0][0] );
      // CASE C2: (i=n-1, j=0, k=n-1)
        P[n-1][0][n-1] = eta()*Q[n-1][0][n-1]
                       + lambda_xx() * ( Q[n-2][0][n-1]   + Q[0][0][n-1] )
                       + lambda_yy() * ( Q[n-1][n-1][n-1] + Q[n-1][1][n-1] )
                       + lambda_zz() * ( Q[n-1][0][n-2]   + Q[n-1][0][0] );
      // CASE C3: (i=n-1, j=n-1, k=n-1)
        P[n-1][n-1][n-1] = eta()*Q[n-1][n-1][n-1]
                         + lambda_xx() * ( Q[n-2][n-1][n-1] + Q[0][n-1][n-1] )
                         + lambda_yy() * ( Q[n-1][n-2][n-1] + Q[n-1][0][n-1] )
                         + lambda_zz() * ( Q[n-1][n-1][n-2] + Q[n-1][n-1][0] );
      // CASE C4: (i=0, j=n-1, k=n-1)
        P[0][n-1][n-1] = eta()*Q[0][n-1][n-1]
                       + lambda_xx() * ( Q[n-1][n-1][n-1] + Q[1][n-1][n-1] )
                       + lambda_yy() * ( Q[0][n-2][n-1]   + Q[0][0][n-1] )
                       + lambda_zz() * ( Q[0][n-1][n-2]   + Q[0][n-1][0] );
      // CASE C5: (0<i<n-1, j=0, k=n-1)
      for(int i=1; i<n-1; i++) {
        P[i][0][n-1] = eta()*Q[i][0][n-1]
                     + lambda_xx() * ( Q[i-1][0][n-1] + Q[i+1][0][n-1] )
                     + lambda_yy() * ( Q[i][n-1][n-1] + Q[i][1][n-1] )
                     + lambda_zz() * ( Q[i][0][n-2]   + Q[i][0][0] );
      } // i
      // CASE C6: (i=n-1, 0<j<n-1, k=n-1)
      for(int j=1; j<n-1; j++) {
        P[n-1][j][n-1] = eta()*Q[n-1][j][n-1]
                       + lambda_xx() * ( Q[n-2][j][n-1]   + Q[0][j][n-1] )
                       + lambda_yy() * ( Q[n-1][j-1][n-1] + Q[n-1][j+1][n-1] )
                       + lambda_zz() * ( Q[n-1][j][n-2]   + Q[n-1][j][0] );
      } // j
      // CASE C7: (0<i<n-1, j=n-1, k=n-1)
      for(int i=1; i<n-1; i++) {
        P[i][n-1][n-1] = eta()*Q[i][n-1][n-1]
                       + lambda_xx() * ( Q[i-1][n-1][n-1] + Q[i+1][n-1][n-1] )
                       + lambda_yy() * ( Q[i][n-2][n-1]   + Q[i][0][n-1] )
                       + lambda_zz() * ( Q[i][n-1][n-2]   + Q[i][n-1][0] );
      } // i
      // CASE C8: (i=0, 0<j<n-1, k=n-1)
      for(int j=1; j<n-1; j++) {
        P[0][j][n-1] = eta()*Q[0][j][n-1]
                     + lambda_xx() * ( Q[n-1][j][n-1] + Q[1][j][n-1] )
                     + lambda_yy() * ( Q[0][j-1][n-1] + Q[0][j+1][n-1] )
                     + lambda_zz() * ( Q[0][j][n-2]   + Q[0][j][0] );
      } // j
      // CASE C9: (0<i<n-1, 0<j<n-1, k=n-1)
      for(int i=1; i<n-1; i++) {
      for(int j=1; j<n-1; j++) {
        P[i][j][n-1] = eta()*Q[i][j][n-1]
                     + lambda_xx() * ( Q[i-1][j][n-1] + Q[i+1][j][n-1] )
                     + lambda_yy() * ( Q[i][j-1][n-1] + Q[i][j+1][n-1] )
                     + lambda_zz() * ( Q[i][j][n-2]   + Q[i][j][0] );
      } // j
      } // i

    return;
}


void TSDCube::Q2P()
{
    for(int i=0; i<NOfDivisions(); i++) {
    for(int j=0; j<NOfDivisions(); j++) {
    for(int k=0; k<NOfDivisions(); k++) {
      P[i][j][k]=Q[i][j][k];
    }
    }
    }
}

void TSDCube::P2Q()
{
    for(int i=0; i<NOfDivisions(); i++) {
    for(int j=0; j<NOfDivisions(); j++) {
    for(int k=0; k<NOfDivisions(); k++) {
      Q[i][j][k]=P[i][j][k];
    }
    }
    }
}

void TSDCube::run()
{

  forever
  {
    if(stopped) return;

    QFile file(outputFileName());
    QTextStream out(&file);
    if(!file.open(QIODevice::WriteOnly | QIODevice::Text))
    {
        emit sendMessage("Failed to open the output file.");
        stopped=true;
        return;
    }
    file.close();

    double p=0;
    double t=0;
    int counter_obs=nObs();
    int counter_DNP=nDNP();
    prepareCalc();
    currentState=K1;
    nextState=L3;

    double rd;


    while(currentState!=E)
    {
      if(stopped) {return;}
      switch(currentState)
      {
      case ERR:
          emit sendMessage(errorMessage());
          currentState=E;
          break;
      case K1:
//          qDebug()<<"cObs: " << counter_obs << "/" <<nObs() << " cDNP: "<< counter_DNP << "/" <<nDNP();
          if(counter_obs==nObs())
          {
              counter_obs=0;
              currentState=L1;
          }
          else if(counter_DNP==nDNP())
          {
              counter_DNP=0;
              currentState=M1;
          }
          else
          {
              currentState=M2;
          }
          break;

      case L1: // record value;
          p=averageP();

          file.open(QIODevice::WriteOnly | QIODevice::Append);
            out << QString::number(t) << " " << QString::number(p) << "\r\n";
          file.close();

          emit sendMessage(QString::number(t)+": "+QString::number(p));
          emit sendData(QString::number(t) +" " +QString::number(p));

          currentState=K1;
          break;

      case M1:
          rd=(double) rand()/RAND_MAX;
//          qDebug()<<"M1: " << rd << " " << exchangeProbability();
          if(rd < exchangeProbability())
          {
            P[0][0][0] += (sourcePolarization()-P[0][0][0])/nVoxel;
//            qDebug() << "!";
          }
          currentState=M2;
          break;

      case M2:
          if(polarizationScaling())
          {
            evolve_P(); // P-scaled spin diffusion
          }
          else
          {
            evolve_NP(); // spin diffusion without polarization scaling
          }
          spinLatticeRelaxation();
          t+=dt();
          if(t>=finalTime()) {currentState=E; break;}
          counter_obs++;
          counter_DNP++;
          currentState=K1;
          break;

      case M3:

          // blank for future extension

          break;
      case L2: // set up polarization
          if(nextState==L3)
          {
              P[0][0][0] += (sourcePolarization()-P[0][0][0])/nVoxel;
          }
          else if(nextState==L4)
          {
              Q[0][0][0] += (sourcePolarization()-Q[0][0][0])/nVoxel;
          }
          currentState=nextState;
          break;
      case L3:
          // spin diffusion (P->Q)
          evolve();
          // spin-lattice relaxation (Q)
          for(int i=0; i<NOfDivisions(); i++) {
          for(int j=0; j<NOfDivisions(); j++) {
          for(int k=0; k<NOfDivisions(); k++) {
            Q[i][j][k] -= R1()*dt()*Q[i][j][k];
            // Here we assume the thermal polarization is negligibly small.
          }
          }
          }

          t+=dt();
          if(t>=finalTime()) {currentState=E; break;}
          counter_obs++;
          nextState=L4;
          if(counter_obs==nObs()) {counter_obs=0; currentState=L1;}
          else {currentState=L2;}
          break;
      case L4:
          // spin diffusion (Q->P)
          evolve2();
          // spin-lattice relaxation (P)
          for(int i=0; i<NOfDivisions(); i++) {
          for(int j=0; j<NOfDivisions(); j++) {
          for(int k=0; k<NOfDivisions(); k++) {
            P[i][j][k] -= R1()*dt()*P[i][j][k];
            // Here we assume the thermal polarization is negligibly small.
          }
          }
          }
          t+=dt();
          if(t>=finalTime()) {currentState=E; break;}
          counter_obs++;
          nextState=L3;
          if(counter_obs==nObs()) {counter_obs=0; currentState=L1;}
          else {currentState=L2;}
          break;
      case E:

          break;

      } // switch
    } // while



    emit sendMessage("Complete!");
    emit calcComplete();
    mutex.lock();
    condition.wait(&mutex); // We let the thread sleep.
    mutex.unlock();
  } // forever
}
