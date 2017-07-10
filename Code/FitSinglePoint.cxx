/*
  Author: Valerio Bertone
 */

// LHAPDF libs
#include "LHAPDF/LHAPDF.h"

// LHAPDF libs
#include "APFEL/APFEL.h"

// APFEL++ libs
#include "apfel/grid.h"
#include "apfel/alphaqcd.h"
#include "apfel/tabulateobject.h"
#include <apfel/dglapbuilder.h>
#include "apfel/structurefunctionbuilder.h"
#include "apfel/timer.h"
#include "apfel/rotations.h"

using namespace std;
using namespace apfel;

int main() {
  // Timer
  Timer t;
  t.start();

  // Open LHAPDF set.
  const string PDFset = "NNPDF30_nlo_as_0118";
  LHAPDF::PDF* PDFs = LHAPDF::mkPDF(PDFset);

  // Retrieve evolution parameters.
  const int    pto   = PDFs->orderQCD();
  const double Qref  = 91.1876;
  //const double asref = PDFs->alphasQ(Qref);
  const double mc    = PDFs->quarkThreshold(4);
  const double mb    = PDFs->quarkThreshold(5);
  const double mt    = PDFs->quarkThreshold(6);
  const double mu0   = PDFs->qMin();
  const vector<double> Masses = {0, 0, 0, mc, mb, mt};
  const vector<double> Thresholds = Masses;

  // x-space grid
  const Grid g{{SubGrid{100,1e-5,3}, SubGrid{60,1e-1,3}, SubGrid{50,6e-1,3}, SubGrid{50,8e-1,3}}};

  // Initialize
  const auto DglapObj = InitializeDglapObjectsQCD(g);
  const auto F2Obj    = InitializeF2ObjectsZM(g);
  const auto FLObj    = InitializeFLObjectsZM(g);
  const auto F3Obj    = InitializeF3ObjectsZM(g);

  // Initial scale PDFs
  const auto fPDFsLH = [&] (double const& x, double const& Q) -> map<int,double>{ return PDFs->xfxQ(x,Q); };
  const auto fPDFs   = [&] (double const& x, double const& Q) -> unordered_map<int,double>{ return PhysToQCDEv(x,Q,fPDFsLH); };

  // Relevant constants for the computation of the EW charges.
  const double MZ         = Qref;
  const double MZ2        = MZ * MZ;
  const double Sin2ThetaW = 0.23126;
  const double VD         = - 0.5 + 2 * Sin2ThetaW / 3;
  const double VU         = + 0.5 - 4 * Sin2ThetaW / 3;
  const vector<double> Vq = {VD, VU, VD, VU, VD, VU};
  const double AD         = - 0.5;
  const double AU         = + 0.5;
  const vector<double> Aq = {AD, AU, AD, AU, AD, AU};

  // Unpolarized positron
  const int ie     = 1; // Positron
  const double Ve  = - 0.5 + 2 * Sin2ThetaW;
  const double Ae  = - 0.5;
  const double pol = 0;   // No polarization

  // Effective charges.
  function<vector<double>(double const&)> fBq = [=] (double const& Q) -> vector<double>
    {
      const double Q2  = Q * Q;
      const double PZ  = Q2 / ( Q2 + MZ2 ) / ( 4 * Sin2ThetaW * ( 1 - Sin2ThetaW ) );
      const double PZ2 = PZ * PZ;
      vector<double> Bq;
      for (auto i = 0; i < (int) Thresholds.size(); i++)
	{
	  const double b = QCh2[i]
	    - 2 * QCh[i] * Vq[i] * ( Ve + ie * pol * Ae ) * PZ
	    + ( Ve * Ve + Ae * Ae + ie * pol * 2 * Ve * Ae )
	    * ( Vq[i] * Vq[i] + Aq[i] * Aq[i] ) * PZ2;
	  Bq.push_back((Q > Thresholds[i] ? b : 0));
	}
      return Bq;
    };
  function<vector<double>(double const&)> fDq = [=] (double const& Q) -> vector<double>
    {
      const double Q2  = Q * Q;
      const double PZ  = Q2 / ( Q2 + MZ2 ) / ( 4 * Sin2ThetaW * ( 1 - Sin2ThetaW ) );
      const double PZ2 = PZ * PZ;
      vector<double> Dq;
      for (auto i = 0; i < (int) Thresholds.size(); i++)
	{
	  const double d = - 2 * QCh[i] * Aq[i] * ( Ae + ie * pol * Ve ) * PZ
	    + 2 * Vq[i] * Aq[i] * ( 2 * Ve * Ae
			    + ie * pol * ( Ve * Ve + Ae * Ae ) ) * PZ2;
	  Dq.push_back((Q > Thresholds[i] ? d : 0));
	}
      return Dq;
    };


  // Cross section
  const auto sigmaNC = [&] (double const& x, double const& Q, double const& y, double const& asref, double const& data) -> double
    {
      // Running coupling
      AlphaQCD Alphas{asref, Qref, Masses, pto};
      const auto as = [&] (double const& Q) -> double{ return Alphas.Evaluate(Q); };

      // Running coupling
      //AlphaQCD AlphasPDF{0.118, Qref, Masses, pto};
      AlphaQCD AlphasPDF{asref, Qref, Masses, pto};
      const auto asPDF = [&] (double const& Q) -> double{ return AlphasPDF.Evaluate(Q); };

      // PDFs
      auto EvolvedPDFs = DglapBuild(DglapObj, fPDFs, mu0, Masses, Thresholds, pto, asPDF);
      const TabulateObject<Set<Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 1, 1000, 3};
      const auto TabPDFs = [&] (double const& x, double const& Q) -> unordered_map<int,double>{ return TabulatedPDFs.EvaluateMapxQ(x,Q); };

      // Structure functions
      const auto F2 = StructureFunctionBuildNC(F2Obj, TabPDFs, Masses, pto, as, fBq);
      const auto FL = StructureFunctionBuildNC(FLObj, TabPDFs, Masses, pto, as, fBq);
      const auto F3 = StructureFunctionBuildNC(F3Obj, TabPDFs, Masses, pto, as, fDq);

      // Definitions
      const double y2 = y * y;
      const double yp = 1 + pow(1 - y, 2);
      const double ym = 1 - pow(1 - y, 2);

      // Compute obsevales
      const double f2 = F2.at(0).Evaluate(Q).Evaluate(x);
      const double fL = FL.at(0).Evaluate(Q).Evaluate(x);
      const double f3 = F3.at(0).Evaluate(Q).Evaluate(x);

      return f2 + ie * ym / yp * f3 - y2 / yp * fL - data;
    };

  // open a file in read mode.
  ifstream infile; 
  infile.open("HERA1+2_NC_e+p_920.dat");
  cout << scientific;
  cout << endl;
  cout << "      Q2      "
       << "      x       "
       << "      y       "
       << endl;
  for (auto i = 0; i < 15; i++)
    {
      double Q2, Q, x, y, Sigma, stat, uncor, unc;
      infile >> Q2 >> x >> y >> Sigma >> stat >> uncor;
      Q = sqrt(Q2);
      //unc = sqrt(stat*stat+uncor*uncor);
      unc = uncor * Sigma / 100;
      
      const auto RootFind = [Q,x,y,sigmaNC] (double const& sigma) -> double
	{
	  double asmin = 0.1;
	  double asmax = 0.13;

	  const int nstepmax = 100;
	  const double acc = 1e-5;

	  double asd = asmin;
	  double asu = asmax;

	enlarge:
	  double fd = sigmaNC(x,Q,y,asd,sigma);
	  double fu = sigmaNC(x,Q,y,asu,sigma);
	  // "fd" and and "fd" must have opposite sign, otherwise you
	  // need to adjust "asmin" and "asmax"
	  if (fd * fu > 0)
	    {
	      double hrange = ( asu - asd ) / 2;
	      asu += hrange;
	      asd -= hrange;
	      //cout << "New range [" << asd << ":" << asu << "]" << endl;
	      //exit(-10);
	      goto enlarge;
	    }

	  for (auto i = 1; i < nstepmax; i++)
	    {
	      double asc = ( asu + asd ) / 2;
	      double fc = sigmaNC(x,Q,y,asc,sigma);
	      if (fd * fc < 0)
		{
		  asu = asc;
		  fu  = fc;
		}
	      else
		{
		  asd = asc;
		  fd  = fc;
		}
	      double reldif = abs(fd - fu);
	      if (reldif < acc)
		return ( asu + asd ) / 2;
	    }
	  return ( asu + asd ) / 2;
	};
      double sigmac  = RootFind(Sigma);
      double dsigmal = sigmac - RootFind(Sigma-unc);
      double dsigmau = RootFind(Sigma+unc) - sigmac;
 
      cout << Q2 << "  " << x << "  " << y << ", alpha_s(MZ) = " << sigmac << " + " << dsigmal << " - " << dsigmau << endl; 
    }
  cout << endl;
  infile.close();

  t.stop();

  return 0;
}
