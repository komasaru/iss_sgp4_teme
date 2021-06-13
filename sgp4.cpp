#include "sgp4.hpp"

namespace iss_sgp4_teme {

// 定数
static constexpr double kE9          = 1.0e9;
static constexpr char   kWgs72Old[]  = "wgs72old";
static constexpr char   kWgs72[]     = "wgs72";
static constexpr char   kWgs84[]     = "wgs84";
static constexpr double kWgs72OldR   =   6378.135;
static constexpr double kWgs72OldMu  = 398600.79964;
static constexpr double kWgs72OldXke =      0.0743669161;
static constexpr double kWgs72OldJ2  =      0.001082616;
static constexpr double kWgs72OldJ3  =     -0.00000253881;
static constexpr double kWgs72OldJ4  =     -0.00000165597;
static constexpr double kWgs72R      =   6378.135;
static constexpr double kWgs72Mu     = 398600.8;
static constexpr double kWgs72J2     =      0.001082616;
static constexpr double kWgs72J3     =     -0.00000253881;
static constexpr double kWgs72J4     =     -0.00000165597;
static constexpr double kWgs84R      =   6378.137;
static constexpr double kWgs84Mu     = 398600.5;
static constexpr double kWgs84J2     =      0.00108262998905;
static constexpr double kWgs84J3     =     -0.00000253215306;
static constexpr double kWgs84J4     =     -0.00000161098761;
static constexpr double kE5          = 1.0e5;
static constexpr double kE7          = 1.0e7;
static constexpr double kPi          = atan(1.0) * 4.0;      // 円周率
static constexpr double kPi2         = kPi * 2.0;            // 円周率 * 2
static constexpr double kMinD        = 1440.0;               // Minutes per day
static constexpr double kXpdotp      = kMinD / (2.0 * kPi);  // 229.1831180523293
static constexpr double kDeg2Rad     = kPi / 180.0;          // 0.0174532925199433

/*
 * @brief      コンストラクタ
 *
 * @param[in]  UT1 (timespec)
 * @param[in]  TLE (vector<string>)
 * @param[in]  測地系 (string; optional)
 */
Sgp4::Sgp4(struct timespec ut1, std::vector<std::string> tle, std::string wgs) {
  this->ut1 = ut1;
  this->tle = tle;
  cst = get_gravconst(wgs);
}

/*
 * @brief   ISS 初期位置・速度の取得
 *          * this procedure is the sgp4 prediction model from space command. 
 *            this is an updated and combined version of sgp4 and sdp4, which 
 *            were originally published separately in spacetrack report #3. 
 *            this version follows the methodology from the aiaa paper (2006) 
 *            describing the history and development of the code.
 *
 * @param   <none>
 * @return  衛星情報 (Satellite)
 */
Satellite Sgp4::twoline2rv(bool afspc_mode) {
  int          nexp;
  int          ibexp;
  unsigned int two_digit_year;
  unsigned int year;
  DateTime     dt;
  Satellite    sat;

  try {
    sat.opsmode = 'i';
    if (afspc_mode) { sat.opsmode = 'a'; }

    // TLE reading
    // 1 NNNNNC NNNNNAAA NNNNN.NNNNNNNN +.NNNNNNNN +NNNNN-N +NNNNN-N N NNNNN
    // 2 NNNNN NNN.NNNN NNN.NNNN NNNNNNN NNN.NNNN NNN.NNNN NN.NNNNNNNNNNNNNN
    // (1st line)
    sat.satnum     = stoi(tle[0].substr( 2,  5));
    two_digit_year = stoi(tle[0].substr(18,  2));
    sat.epochdays  = stod(tle[0].substr(20, 12));
    sat.ndot       = stod(tle[0].substr(33, 10));
    sat.nddot      = stod(tle[0].substr(44,  6)) / kE5;
    nexp           = stoi(tle[0].substr(50,  2));
    sat.bstar      = stod(tle[0].substr(53,  6)) / kE5;
    ibexp          = stoi(tle[0].substr(59,  2));
    // (2nd line)
    sat.inclo = stod(tle[1].substr( 8,  8));
    sat.nodeo = stod(tle[1].substr(17,  8));
    sat.ecco  = stoi(tle[1].substr(26,  7)) / kE7;
    sat.argpo = stod(tle[1].substr(34,  8));
    sat.mo    = stod(tle[1].substr(43,  8));
    sat.no    = stod(tle[1].substr(52, 11)) / kXpdotp;

    // ---- find no, ndot, nddot ----
    sat.nddot *= std::pow(10.0, nexp);
    sat.bstar *= std::pow(10.0, ibexp);

    // ---- convert to sgp4 units ----
    sat.ndot  = sat.ndot  / (kXpdotp * kMinD);
    sat.nddot = sat.nddot / (kXpdotp * kMinD * kMinD);

    // ---- find standard orbital elements ----
    sat.inclo *= kDeg2Rad;
    sat.nodeo *= kDeg2Rad;
    sat.argpo *= kDeg2Rad;
    sat.mo    *= kDeg2Rad;

    // ----------------------------------------------------------------
    // find sgp4epoch time of element set
    // remember that sgp4 uses units of days from 0 jan 1950 (sgp4epoch)
    // and minutes from the epoch (time)
    // ----------------------------------------------------------------

    // ---- temp fix for years from 1957-2056 ----
    // ---- correct fix will occur when year is 4-digit in tle ----
    if (two_digit_year < 57) {
      year = two_digit_year + 2000;
    } else {
      year = two_digit_year + 1900;
    }

    // ---- year + days -> year + month + day + hour + minute + second ----
    dt = days2ymdhms(year, sat.epochdays);
    sat.epochyr    = year;
    sat.jdsatepoch = jday(dt);

    // ---- initialize the orbit at sgp4epoch ----
    sgp4init(sat);
  } catch (...) {
    throw;
  }

  return sat;
}  // twoline2rv

/*
 * @brief       指定 UT1 の ISS 位置・速度の取得
 *              * Return a position and velocity vector for a given date and 
 *                time.
 *
 * @param[ref]  sat (Satellite)
 * @return      位置・速度 (PvTeme)
 */
PvTeme Sgp4::propagate(Satellite& sat) {
  std::string ut1_str;
  DateTime    dt;
  double      j;
  double      m;
  PvTeme      teme = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

  try {
    ut1_str = gen_time_str(ut1);
    dt.year   = stoi(ut1_str.substr( 0, 4));
    dt.month  = stoi(ut1_str.substr( 5, 2));
    dt.day    = stoi(ut1_str.substr( 8, 2));
    dt.hour   = stoi(ut1_str.substr(11, 2));
    dt.minute = stoi(ut1_str.substr(14, 2));
    dt.second = stod(ut1_str.substr(17, 2)) + ut1.tv_nsec / kE9;
    j = jday(dt);
    m = (j - sat.jdsatepoch) * kMinD;
    teme = sgp4(m, sat);
  } catch (...) {
    throw;
  }

  return teme;
}

/********************************************
 **** 以下、 private function/procedures ****
 ********************************************/

/*
 * @brief      定数取得
 *
 * @param[in]  測地系 ("wgs72old"|"wgs72"|"wgs84") (string)
 *             (default: "wgs84")
 * @return     定数セット (Const)
 */
Const Sgp4::get_gravconst(std::string wgs) {
  Const cst;

  try {
    if (wgs == kWgs72Old) {
      cst.r   = kWgs72OldR;
      cst.mu  = kWgs72OldMu;
      cst.xke = kWgs72OldXke;
      cst.j2  = kWgs72OldJ2;
      cst.j3  = kWgs72OldJ3;
      cst.j4  = kWgs72OldJ4;
    } else if (wgs == kWgs72) {
      cst.r   = kWgs72R;
      cst.mu  = kWgs72Mu;
      cst.xke = 60.0 / sqrt(cst.r * cst.r * cst.r / cst.mu);
      cst.j2  = kWgs72J2;
      cst.j3  = kWgs72J3;
      cst.j4  = kWgs72J4;
    } else if (wgs == kWgs84) {
      cst.r   = kWgs84R;
      cst.mu  = kWgs84Mu;
      cst.xke = 60.0 / sqrt(cst.r * cst.r * cst.r / cst.mu);
      cst.j2  = kWgs84J2;
      cst.j3  = kWgs84J3;
      cst.j4  = kWgs84J4;
    }
    cst.tumin = 1.0 / cst.xke;
    cst.j3oj2 = cst.j3 / cst.j2;
  } catch (...) {
    throw;
  }

  return cst;
}

/*
 * @brief       SGP4 initialization
 *              * this procedure initializes variables for sgp4.
 *
 * @param[ref]  sat (Satellite)
 * @return      <none>
 */
void Sgp4::sgp4init(Satellite& sat) {
  double epoch;
  double ss;
  double qzms2ttemp;
  double qzms2t;
  double x2o3;
  double ainv;
  double ao;
  double con42;
  double cosio;
  double cosio2;
  double cc1sq;
  double eccsq;
  double omeosq;
  double posq;
  double rp;
  double rteosq;
  double sinio;
  double perige;
  double qzms24;
  double sfour;
  double qzms24temp;
  double pinvsq;
  double coef;
  double coef1;
  double eeta;
  double etasq;
  double psisq;
  double tsi;
  double cc2;
  double cc3;
  double cosio4;
  double temp;
  double temp1;
  double temp2;
  double temp3;
  double temp4;
  double xhdot1;
  double xpidot;
  double delmotemp;
  double tc;
  double inclm;
  double argpm;
  double nodem;
  double mm;
  double dndt;
  double snodm;
  double cnodm;
  double sinim;
  double cosim;
  double sinomm;
  double cosomm;
  double day;
  double em;
  double emsq;
  double gam;
  double rtemsq;
  double s1;
  double s2;
  double s3;
  double s4;
  double s5;
  double s6;
  double s7;
  double ss1;
  double ss2;
  double ss3;
  double ss4;
  double ss5;
  double ss6;
  double ss7;
  double sz1;
  double sz2;
  double sz3;
  double sz11;
  double sz12;
  double sz13;
  double sz21;
  double sz22;
  double sz23;
  double sz31;
  double sz32;
  double sz33;
  double nm;
  double z1;
  double z2;
  double z3;
  double z11;
  double z12;
  double z13;
  double z21;
  double z22;
  double z23;
  double z31;
  double z32;
  double z33;
  PvTeme teme;

  try {
    // ---- initialization ----
    // sgp4fix divisor for divide by zero check on inclination
    // the old check used 1.0 + cos(pi-1.0e-9), but { compared it to
    // 1.5 e-12, so the threshold was changed to 1.5e-12 for consistency
    epoch = sat.jdsatepoch - 2433281.5;
    temp4 = 1.5e-12;

    // ---- earth constants ----
    // sgp4fix identify constants and allow alternate values
    ss = 78.0 / cst.r + 1.0;
    // sgp4fix use multiply for speed instead of pow
    qzms2ttemp = (120.0 - 78.0) / cst.r;
    qzms2t     = qzms2ttemp * qzms2ttemp * qzms2ttemp * qzms2ttemp;
    x2o3       = 2.0 / 3.0;

    sat.init = 'y';
    sat.t    = 0.0;

    initl(epoch,
        ainv, ao, sat.con41, con42, cosio, cosio2, eccsq,
        omeosq, posq, rp, rteosq, sinio, sat);

    // sgp4fix remove this check as it is unnecessary
    // the mrt check in sgp4 handles decaying satellite cases even if the starting
    // condition is below the surface of te earth
    if (rp < 1.0) {
      sat.error = 5;
      return;
    }

    if (omeosq >= 0.0 || sat.no >= 0.0) {
      sat.isimp = 0;
      if (rp < 220.0 / cst.r + 1.0) { sat.isimp = 1; }
      sfour  = ss;
      qzms24 = qzms2t;
      perige = (rp - 1.0) * cst.r;

      // for perigees below 156 km, s and qoms2t are altered -
      if (perige < 156.0) {
        sfour = perige - 78.0;
        if (perige < 98.0) { sfour = 20.0;}
        // sgp4fix use multiply for speed instead of pow
        qzms24temp = (120.0 - sfour) / cst.r;
        qzms24    = qzms24temp * qzms24temp * qzms24temp * qzms24temp;
        sfour     = sfour / cst.r + 1.0;
      }

      pinvsq  = 1.0 / posq;
      tsi     = 1.0 / (ao - sfour);
      sat.eta = ao * sat.ecco * tsi;
      etasq   = sat.eta * sat.eta;
      eeta    = sat.ecco * sat.eta;
      psisq   = fabs(1.0 - etasq);
      coef    = qzms24 * pow(tsi, 4.0);
      coef1   = coef / pow(psisq, 3.5);
      cc2     = coef1 * sat.no
              * (ao * (1.0 + 1.5 * etasq + eeta * (4.0 + etasq))
              + 0.375 * cst.j2 * tsi / psisq * sat.con41
              * (8.0 + 3.0 * etasq * (8.0 + etasq)));
      sat.cc1 = sat.bstar * cc2;
      cc3     = 0.0;
      if (sat.ecco > 1.0e-4) {
        cc3 = -2.0 * coef * tsi * cst.j3oj2 * sat.no * sinio / sat.ecco;
      }
      sat.x1mth2 = 1.0 - cosio2;
      sat.cc4 = 2.0 * sat.no * coef1 * ao * omeosq
              * (sat.eta * (2.0 + 0.5 * etasq)
              + sat.ecco * (0.5 + 2.0 * etasq)
              - cst.j2 * tsi / (ao * psisq)
              * (-3.0 * sat.con41 * (1.0 - 2.0 * eeta
              + etasq * (1.5 - 0.5 * eeta))
              + 0.75 * sat.x1mth2 * (2.0 * etasq
              - eeta * (1.0 + etasq)) * cos(2.0 * sat.argpo)));
      sat.cc5 = 2.0 * coef1 * ao * omeosq
              * (1.0 + 2.75 * (etasq + eeta) + eeta * etasq);
      cosio4  = cosio2 * cosio2;
      temp1   = 1.5 * cst.j2 * pinvsq * sat.no;
      temp2   = 0.5 * temp1 * cst.j2 * pinvsq;
      temp3   = -0.46875 * cst.j4 * pinvsq * pinvsq * sat.no;
      sat.mdot    = sat.no + 0.5 * temp1 * rteosq * sat.con41
                  + 0.0625 * temp2 * rteosq
                  * (13.0 - 78.0 * cosio2 + 137.0 * cosio4);
      sat.argpdot = (-0.5 * temp1 * con42 + 0.0625 * temp2
                  * (7.0 - 114.0 * cosio2 + 395.0 * cosio4)
                  + temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4));
      xhdot1      = -temp1 * cosio;
      sat.nodedot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosio2)
                  + 2.0 * temp3 * (3.0 - 7.0 * cosio2)) * cosio;
      xpidot      = sat.argpdot + sat.nodedot;
      sat.omgcof  = sat.bstar * cc3 * cos(sat.argpo);
      sat.xmcof   = 0.0;
      if (sat.ecco > 1.0e-4) { sat.xmcof = -x2o3 * coef * sat.bstar / eeta; }
      sat.nodecf  = 3.5 * omeosq * xhdot1 * sat.cc1;
      sat.t2cof   = 1.5 * sat.cc1;
      // sgp4fix for divide by zero with xinco = 180 deg
      if (fabs(cosio + 1.0) > 1.5e-12) {
        sat.xlcof = -0.25 * cst.j3oj2 * sinio
                  * (3.0 + 5.0 * cosio) / (1.0 + cosio);
      } else {
        sat.xlcof = -0.25 * cst.j3oj2 * sinio * (3.0 + 5.0 * cosio) / temp4;
      }
      sat.aycof = -0.5 * cst.j3oj2 * sinio;
      // sgp4fix use multiply for speed instead of pow
      delmotemp  = 1.0 + sat.eta * cos(sat.mo);
      sat.delmo  = delmotemp * delmotemp * delmotemp;
      sat.sinmao = sin(sat.mo);
      sat.x7thm1 = 7.0 * cosio2 - 1.0;

      // ---- deep space initialization ----
      if (kPi2 / sat.no >= 225.0) {
        sat.method = 'd';
        sat.isimp  = 1;
        tc         =  0.0;
        inclm      = sat.inclo;

        dscom(
          epoch, tc,
          snodm, cnodm, sinim, cosim, sinomm, cosomm, day,
          em, emsq, gam, rtemsq, s1, s2, s3, s4, s5,
          s6, s7, ss1, ss2, ss3, ss4, ss5, ss6, ss7, sz1,
          sz2, sz3, sz11, sz12, sz13, sz21, sz22, sz23,
          sz31, sz32, sz33, nm, z1, z2, z3,
          z11, z12, z13, z21, z22, z23, z31, z32, z33,
          sat);

        argpm = 0.0;
        nodem = 0.0;
        mm    = 0.0;

        dsinit(
          cosim, emsq, s1, s2, s3, s4, s5, sinim,
          ss1, ss2, ss3, ss4, ss5,
          sz1, sz3, sz11, sz13, sz21,
          sz23, sz31, sz33, tc, xpidot,
          z1, z3, z11, z13, z21, z23, z31, z33,
          eccsq, em, argpm, inclm, mm, nm, nodem, dndt,
          sat);
      }

      // ---- set variables if not deep space ----
      if (sat.isimp != 1) {
        cc1sq     = sat.cc1 * sat.cc1;
        sat.d2    = 4.0 * ao * tsi * cc1sq;
        temp      = sat.d2 * tsi * sat.cc1 / 3.0;
        sat.d3    = (17.0 * ao + sfour) * temp;
        sat.d4    = 0.5 * temp * ao * tsi * (221.0 * ao + 31.0 * sfour)
                  * sat.cc1;
        sat.t3cof = sat.d2 + 2.0 * cc1sq;
        sat.t4cof = 0.25 * (3.0 * sat.d3 + sat.cc1
                  * (12.0 * sat.d2 + 10.0 * cc1sq));
        sat.t5cof =  0.2 * (3.0 * sat.d4 + 12.0 * sat.cc1 * sat.d3
                  +  6.0 * sat.d2 * sat.d2 + 15.0 * cc1sq
                  * (2.0 * sat.d2 + cc1sq));
      }
    }

    teme = sgp4(0.0, sat);

    sat.init = 'n';
  } catch (...) {
    throw;
  }
}  // sgp4init

/*
 * @brief       SGP4 propagator initialization
 *              * this procedure initializes the spg4 propagator. all the 
 *                initialization is consolidated here instead of having 
 *                multiple loops inside other routines.
 *
 * @param[in]       epoch (double)
 * @param[ref]       ainv (double)
 * @param[ref]         ao (double)
 * @param[ref]      con41 (double)
 * @param[ref]      con42 (double)
 * @param[ref]      cosio (double)
 * @param[ref]     cosio2 (double)
 * @param[ref]      eccsq (double)
 * @param[ref]     omeosq (double)
 * @param[ref]       posq (double)
 * @param[ref]         rp (double)
 * @param[ref]     rteosq (double)
 * @param[ref]      sinio (double)
 * @param[ref]        sat (Satellite)
 * @return      <none>
 */
void Sgp4::initl(
    double epoch,
    double& ainv, double& ao, double& con41, double& con42, double& cosio,
    double& cosio2, double& eccsq, double& omeosq, double& posq, double& rp,
    double& rteosq, double& sinio, Satellite& sat) {
  int    ds70;
  double x2o3;
  double ak;
  double d1;
  double del_;
  double adel;
  double po;
  double ts70;
  double tfrac;
  double c1;
  double thgr70;
  double fk5r;
  double c1p2p;

  try {
    // sgp4fix use old way of finding gst

    // ---- earth constants ----
    x2o3 = 2.0 / 3.0;

    // ---- calculate auxillary epoch quantities ----
    eccsq  = sat.ecco * sat.ecco;
    omeosq = 1.0 - eccsq;
    rteosq = sqrt(omeosq);
    cosio  = cos(sat.inclo);
    cosio2 = cosio * cosio;

    // ---- un-kozai the mean motion ----
    ak     = pow(cst.xke / sat.no, x2o3);
    d1     = 0.75 * cst.j2 * (3.0 * cosio2 - 1.0) / (rteosq * omeosq);
    del_   = d1 / (ak * ak);
    adel   = ak * (1.0 - del_ * del_ - del_
           * (1.0 / 3.0 + 134.0 * del_ * del_ / 81.0));
    del_   = d1 / (adel * adel);
    sat.no /= 1.0 + del_;
    ao     = pow(cst.xke / sat.no, x2o3);
    sinio  = sin(sat.inclo);
    po     = ao * omeosq;
    con42  = 1.0 - 5.0 * cosio2;
    con41  = - con42 - cosio2 - cosio2;
    ainv   = 1.0 / ao;
    posq   = po * po;
    rp     = ao * (1.0 - sat.ecco);
    sat.method = 'n';

    // ---- sgp4fix modern approach to finding sidereal time ----
    if (sat.opsmode == 'a') {
      // sgp4fix use old way of finding gst
      // count integer number of days from 0 jan 1970
      ts70  = epoch - 7305.0;
      ds70  = int(ts70 + 1.0e-8);
      tfrac = ts70 - double(ds70);
      // find greenwich location at epoch
      c1     = 1.72027916940703639e-2;
      thgr70 = 1.7321343856509374;
      fk5r   = 5.07551419432269442e-15;
      c1p2p  = c1 + kPi2;
      sat.gsto = fmod(thgr70 + c1 * ds70 + c1p2p * tfrac + ts70 * ts70 * fk5r,
                      kPi2);
      if (sat.gsto < 0.0) { sat.gsto += kPi2; }
    } else {
      sat.gsto = gstime(epoch + 2433281.5);
    }
  } catch (...) {
    throw;
  }
}  // initl

/*
 * @brief       Deep space common items
 *              * this procedure provides deep space common items used by both 
 *                the secular and periodics subroutines.  input is provided as 
 *                shown. this routine used to be called dpper, but the 
 *                functions inside weren't well organized.
 *
 * @param[in]   epoch (double)
 * @param[in]   tc (double)
 * @param[ref]  snodm, cnodm, sinim, cosim, sinomm, cosomm,
 *              day, em, emsq, gam, rtemsq,
 *              s1, s2, s3, s4, s5, s6, s7,
 *              ss1, ss2, ss3, ss4, ss5, ss6, ss7,
 *              sz1, sz2, sz3,
 *              sz11, sz12, sz13,
 *              sz21, sz22, sz23,
 *              sz31, sz32, sz33, nm,
 *              z1, z2, z3,
 *              z11, z12, z13,
 *              z21, z22, z23,
 *              z31, z32, z33 (double)
 * @param[ref]  sat (Satellite)
 */
void Sgp4::dscom(
    double  epoch,  double tc,
    double& snodm,  double& cnodm,  double& sinim, double& cosim,
    double& sinomm, double& cosomm, double& day,
    double& em,     double& emsq,   double& gam,   double& rtemsq,
    double& s1,     double& s2,     double& s3,    double& s4,     double& s5,
    double& s6,     double& s7,
    double& ss1,    double& ss2,    double& ss3,   double& ss4,    double& ss5,
    double& ss6,    double& ss7,
    double& sz1,    double& sz2,    double& sz3,
    double& sz11,   double& sz12,   double& sz13,
    double& sz21,   double& sz22,   double& sz23,
    double& sz31,   double& sz32,   double& sz33,  double& nm,
    double& z1,     double& z2,     double& z3,
    double& z11,    double& z12,    double& z13,
    double& z21,    double& z22,    double& z23,
    double& z31,    double& z32,    double& z33,
    Satellite& sat) {
  // constants
  static constexpr double kZes    =  0.01675;
  static constexpr double kZel    =  0.05490;
  static constexpr double kC1ss   =  2.9864797e-6;
  static constexpr double kC1l    =  4.7968065e-7;
  static constexpr double kZsinis =  0.39785416;
  static constexpr double kZcosis =  0.91744867;
  static constexpr double kZcosgs =  0.1945905;
  static constexpr double kZsings = -0.98088458;
  // variables
  unsigned int lsflg;
  double betasq;
  double ctem;
  double stem;
  double xnodce;
  double zx;
  double zy;
  double zcosgl;
  double zcoshl;
  double zcosil;
  double zsingl;
  double zsinhl;
  double zsinil;
  double cc;
  double xnoi;
  double zcosg;
  double zcosh;
  double zcosi;
  double zsing;
  double zsinh;
  double zsini;
  double a[10];
  double x[8];

  try {
    // ---- local variables ----
    nm     = sat.no;
    em     = sat.ecco;
    snodm  = sin(sat.nodeo);
    cnodm  = cos(sat.nodeo);
    sinomm = sin(sat.argpo);
    cosomm = cos(sat.argpo);
    sinim  = sin(sat.inclo);
    cosim  = cos(sat.inclo);
    emsq   = em * em;
    betasq = 1.0 - emsq;
    rtemsq = sqrt(betasq);

    // ---- initialize lunar solar terms ----
    sat.peo   = 0.0;
    sat.pinco = 0.0;
    sat.plo   = 0.0;
    sat.pgho  = 0.0;
    sat.pho   = 0.0;
    day    = epoch + 18261.5 + tc / 1440.0;
    xnodce = fmod(4.5236020 - 9.2422029e-4 * day, kPi2);
    stem   = sin(xnodce);
    ctem   = cos(xnodce);
    zcosil = 0.91375164 - 0.03568096 * ctem;
    zsinil = sqrt(1.0 - zcosil * zcosil);
    zsinhl = 0.089683511 * stem / zsinil;
    zcoshl = sqrt(1.0 - zsinhl * zsinhl);
    gam    = 5.8351514 + 0.0019443680 * day;
    zx     = 0.39785416 * stem / zsinil;
    zy     = zcoshl * ctem + 0.91744867 * zsinhl * stem;
    zx     = atan2(zx, zy);
    zx     = gam + zx - xnodce;
    zcosgl = cos(zx);
    zsingl = sin(zx);

    // ---- do solar terms ----
    zcosg = kZcosgs;
    zsing = kZsings;
    zcosi = kZcosis;
    zsini = kZsinis;
    zcosh = cnodm;
    zsinh = snodm;
    cc    = kC1ss;
    xnoi  = 1.0 / nm;

    for (lsflg = 0; lsflg < 2; ++lsflg) {
      a[0] =  zcosg * zcosh + zsing * zcosi * zsinh;
      a[2] = -zsing * zcosh + zcosg * zcosi * zsinh;
      a[6] = -zcosg * zsinh + zsing * zcosi * zcosh;
      a[7] =  zsing * zsini;
      a[8] =  zsing * zsinh + zcosg * zcosi * zcosh;
      a[9] =  zcosg * zsini;
      a[1] =  cosim * a[6] + sinim * a[7];
      a[3] =  cosim * a[8] + sinim * a[9];
      a[4] = -sinim * a[6] + cosim * a[7];
      a[5] = -sinim * a[8] + cosim * a[9];
      x[0] =  a[0] * cosomm + a[1] * sinomm;
      x[1] =  a[2] * cosomm + a[3] * sinomm;
      x[2] = -a[0] * sinomm + a[1] * cosomm;
      x[3] = -a[2] * sinomm + a[3] * cosomm;
      x[4] =  a[4] * sinomm;
      x[5] =  a[5] * sinomm;
      x[6] =  a[4] * cosomm;
      x[7] =  a[5] * cosomm;
      z31 = 12.0 * x[0] * x[0] - 3.0 * x[2] * x[2];
      z32 = 24.0 * x[0] * x[1] - 6.0 * x[2] * x[3];
      z33 = 12.0 * x[1] * x[1] - 3.0 * x[3] * x[3];
      z1  =  3.0 * (a[0] * a[0] + a[1] * a[1]) + z31 * emsq;
      z2  =  6.0 * (a[0] * a[2] + a[1] * a[3]) + z32 * emsq;
      z3  =  3.0 * (a[2] * a[2] + a[3] * a[3]) + z33 * emsq;
      z11 = -6.0 * a[0] * a[4]
          + emsq * (-24.0 * x[0] * x[6] -6.0 * x[2] * x[4]);
      z12 = -6.0 * (a[0] * a[5] + a[2] * a[4])
          + emsq * (-24.0 * (x[1] * x[6] + x[0] * x[7])
          -  6.0 * (x[2] * x[5] + x[3] * x[4]));
      z13 = -6.0 * a[2] * a[5]
          + emsq * (-24.0 * x[1] * x[7] - 6.0 * x[3] * x[5]);
      z21 =  6.0 * a[1] * a[4]
          + emsq * ( 24.0 * x[0] * x[4] - 6.0 * x[2] * x[6]);
      z22 =  6.0 * (a[3] * a[4] + a[1] * a[5])
          + emsq * ( 24.0 * (x[1] * x[4] + x[0] * x[5])
          -  6.0 * (x[3] * x[6] + x[2] * x[7]));
      z23 =  6.0 * a[3] * a[5]
          + emsq * ( 24.0 * x[1] * x[5] - 6.0 * x[3] * x[7]);
      z1  += z1 + betasq * z31;
      z2  += z2 + betasq * z32;
      z3  += z3 + betasq * z33;
      s3  = cc * xnoi;
      s2  = -0.5 * s3 / rtemsq;
      s4  = s3 * rtemsq;
      s1  = -15.0 * em * s4;
      s5  = x[0] * x[2] + x[1] * x[3];
      s6  = x[1] * x[2] + x[0] * x[3];
      s7  = x[1] * x[3] - x[0] * x[2];

      // ---- do lunar terms ----
      if (lsflg == 0) {
        ss1   = s1;
        ss2   = s2;
        ss3   = s3;
        ss4   = s4;
        ss5   = s5;
        ss6   = s6;
        ss7   = s7;
        sz1   = z1;
        sz2   = z2;
        sz3   = z3;
        sz11  = z11;
        sz12  = z12;
        sz13  = z13;
        sz21  = z21;
        sz22  = z22;
        sz23  = z23;
        sz31  = z31;
        sz32  = z32;
        sz33  = z33;
        zcosg = zcosgl;
        zsing = zsingl;
        zcosi = zcosil;
        zsini = zsinil;
        zcosh = zcoshl * cnodm + zsinhl * snodm;
        zsinh = snodm * zcoshl - cnodm * zsinhl;
        cc    = kC1l;
      }
    }

    sat.zmol = fmod(4.7199672 + 0.22997150  * day - gam, kPi2);
    sat.zmos = fmod(6.2565837 + 0.017201977 * day, kPi2);

    // ---- do solar terms ----
    sat.se2  =   2.0 * ss1 * ss6;
    sat.se3  =   2.0 * ss1 * ss7;
    sat.si2  =   2.0 * ss2 * sz12;
    sat.si3  =   2.0 * ss2 * (sz13 - sz11);
    sat.sl2  =  -2.0 * ss3 * sz2;
    sat.sl3  =  -2.0 * ss3 * (sz3 - sz1);
    sat.sl4  =  -2.0 * ss3 * (-21.0 - 9.0 * emsq) * kZes;
    sat.sgh2 =   2.0 * ss4 * sz32;
    sat.sgh3 =   2.0 * ss4 * (sz33 - sz31);
    sat.sgh4 = -18.0 * ss4 * kZes;
    sat.sh2  =  -2.0 * ss2 * sz22;
    sat.sh3  =  -2.0 * ss2 * (sz23 - sz21);

    // ---- do lunar terms ----
    sat.ee2  =   2.0 * s1 * s6;
    sat.e3   =   2.0 * s1 * s7;
    sat.xi2  =   2.0 * s2 * z12;
    sat.xi3  =   2.0 * s2 * (z13 - z11);
    sat.xl2  =  -2.0 * s3 * z2;
    sat.xl3  =  -2.0 * s3 * (z3 - z1);
    sat.xl4  =  -2.0 * s3 * (-21.0 - 9.0 * emsq) * kZel;
    sat.xgh2 =   2.0 * s4 * z32;
    sat.xgh3 =   2.0 * s4 * (z33 - z31);
    sat.xgh4 = -18.0 * s4 * kZel;
    sat.xh2  =  -2.0 * s2 * z22;
    sat.xh3  =  -2.0 * s2 * (z23 - z21);
  } catch (...) {
    throw;
  }
}  // dscom

/*
 * @brief       Deep space contributions
 *              * this procedure provides deep space contributions to mean 
 *                motion dot due to geopotential resonance with half day and one 
 *                day orbits.
 *
 * @param[in]   cosim, emsq (double)
 * @param[in]   s1,   s2,   s3,   s4,   s5, sinim,
 *              ss1,  ss2,  ss3,  ss4,  ss5,
 *              sz1,  sz3,  sz11, sz13, sz21,
 *              sz23, sz31, sz33, tc,   xpidot,
 *              z1,   z3,   z11,  z13,
 *              z21,  z23,  z31,  z33,  eccsq (double)
 * @param[ref]  em,  argpm, inclm, mm,  nm, nodem, dndt (double)
 * @param[ref]  sat (Satellite )
 */
void Sgp4::dsinit(
    double cosim,  double emsq,
    double s1,   double s2,     double s3,     double s4,   double s5,  double sinim,
    double ss1,  double ss2,    double ss3,    double ss4,  double ss5,
    double sz1,  double sz3,    double sz11,   double sz13, double sz21,
    double sz23, double sz31,   double sz33,   double tc,   double xpidot,
    double z1,   double z3,     double z11,    double z13,
    double z21,  double z23,    double z31,    double z33,  double eccsq,
    double& em,  double& argpm, double& inclm, double& mm,  double& nm,
    double& nodem, double& dndt,
    Satellite& sat) {
  double q22;
  double q31;
  double q33;
  double root22;
  double root32;
  double root44;
  double root52;
  double root54;
  double rptim;
  double x2o3;
  double znl;
  double zns;
  double ses;
  double sghs;
  double sgs;
  double shs;
  double sis;
  double sls;
  double sghl;
  double shll;
  double theta;
  double aonv;
  double cosisq;
  double emo;
  double emsqo;
  double eoc;
  double g201;
  double g200;
  double g211;
  double g300;
  double g310;
  double g322;
  double g410;
  double g422;
  double g520;
  double g521;
  double g532;
  double g533;
  double f220;
  double f221;
  double f311;
  double f321;
  double f322;
  double f330;
  double f441;
  double f442;
  double f522;
  double f523;
  double f542;
  double f543;
  double sini2;
  double ainv2;
  double xno2;
  double temp;
  double temp1;

  try {
    q22    = 1.7891679e-6;
    q31    = 2.1460748e-6;
    q33    = 2.2123015e-7;
    root22 = 1.7891679e-6;
    root44 = 7.3636953e-9;
    root54 = 2.1765803e-9;
    rptim  = 4.37526908801129966e-3;  // equates to 7.29211514668855e-5 rad/sec
    root32 = 3.7393792e-7;
    root52 = 1.1428639e-7;
    x2o3   = 2.0 / 3.0;
    znl    = 1.5835218e-4;
    zns    = 1.19459e-5;

    // ---- deep space initialization ----
    sat.irez = 0;
    if (0.0034906585 < nm && nm < 0.0052359877) {
      sat.irez = 1;
    } else if (8.26e-3 <= nm && nm <= 9.24e-3 && em >= 0.5) {
      sat.irez = 2;
    }

    // ---- do solar terms ----
    ses  =  ss1 * zns * ss5;
    sis  =  ss2 * zns * (sz11 + sz13);
    sls  = -zns * ss3 * (sz1 + sz3 - 14.0 - 6.0 * emsq);
    sghs =  ss4 * zns * (sz31 + sz33 - 6.0);
    shs  = -zns * ss2 * (sz21 + sz23);
    // sgp4fix for 180 deg incl
    if (inclm < 5.2359877e-2 || inclm > kPi - 5.2359877e-2) {
      shs = 0.0;
    } else if (sinim != 0.0) {
      shs /= sinim;
    }
    sgs  = sghs - cosim * shs;

    // ---- do lunar terms ----
    sat.dedt = ses + s1 * znl * s5;
    sat.didt = sis + s2 * znl * (z11 + z13);
    sat.dmdt = sls - znl * s3 * (z1 + z3 - 14.0 - 6.0 * emsq);
    sghl     = s4 * znl * (z31 + z33 - 6.0);
    shll     = -znl * s2 * (z21 + z23);
    // sgp4fix for 180 deg incl
    if (inclm < 5.2359877e-2 || inclm > kPi - 5.2359877e-2) {
      shll = 0.0;
    }
    sat.domdt = sgs + sghl;
    sat.dnodt = shs;
    if (sinim != 0.0) {
      sat.domdt -= cosim / sinim * shll;
      sat.dnodt += shll  / sinim;
    }

    // ---- calculate deep space resonance effects ----
    dndt   = 0.0;
    theta  = fmod(sat.gsto + tc * rptim, kPi2);
    em    += sat.dedt  * sat.t;
    inclm += sat.didt  * sat.t;
    argpm += sat.domdt * sat.t;
    nodem += sat.dnodt * sat.t;
    mm    += sat.dmdt  * sat.t;

    // ---- initialize the resonance terms ----
    if (sat.irez != 0) { aonv = pow((nm / cst.xke), x2o3); }

    // ---- geopotential resonance for 12 hour orbits ----
    if (sat.irez == 2) {
      cosisq = cosim * cosim;
      emo    = em;
      em     = sat.ecco;
      emsqo  = emsq;
      emsq   = eccsq;
      eoc    = em * emsq;
      g201   = -0.306 - (em - 0.64) * 0.440;

      if (em <= 0.65) {
        g211 =    3.616  -   13.2470 * em +   16.2900 * emsq;
        g310 =  -19.302  +  117.3900 * em -  228.4190 * emsq +  156.5910 * eoc;
        g322 =  -18.9068 +  109.7927 * em -  214.6334 * emsq +  146.5816 * eoc;
        g410 =  -41.122  +  242.6940 * em -  471.0940 * emsq +  313.9530 * eoc;
        g422 = -146.407  +  841.8800 * em - 1629.014  * emsq + 1083.4350 * eoc;
        g520 = -532.114  +  3017.977 * em - 5740.032  * emsq + 3708.2760 * eoc;
      } else {
        g211 =   -72.099 +   331.819 * em -   508.738 * emsq +   266.724 * eoc;
        g310 =  -346.844 +  1582.851 * em -  2415.925 * emsq +  1246.113 * eoc;
        g322 =  -342.585 +  1554.908 * em -  2366.899 * emsq +  1215.972 * eoc;
        g410 = -1052.797 +  4758.686 * em -  7193.992 * emsq +  3651.957 * eoc;
        g422 = -3581.690 + 16178.110 * em - 24462.770 * emsq + 12422.520 * eoc;
        if (em > 0.715) {
          g520 = -5149.66 + 29936.92 * em - 54087.36  * emsq + 31324.56  * eoc;
        } else {
          g520 =  1464.74 -  4664.75 * em +  3763.64  * emsq;
        }
      }

      if (em < 0.7) {
        g533 = -919.22770 + 4988.6100 * em - 9064.7700 * emsq + 5542.21   * eoc;
        g521 = -822.71072 + 4568.6173 * em - 8491.4146 * emsq + 5337.524  * eoc;
        g532 = -853.66600 + 4690.2500 * em - 8624.7700 * emsq + 5341.4    * eoc;
      } else {
        g533 = -37995.780 + 161616.52 * em - 229838.20 * emsq + 109377.94 * eoc;
        g521 = -51752.104 + 218913.95 * em - 309468.16 * emsq + 146349.42 * eoc;
        g532 = -40023.880 + 170470.89 * em - 242699.48 * emsq + 115605.82 * eoc;
      }

      sini2     = sinim * sinim;
      f220      = 0.75 * (1.0 + 2.0 * cosim + cosisq);
      f221      = 1.5 * sini2;
      f321      =  1.875 * sinim * (1.0 - 2.0 * cosim - 3.0 * cosisq);
      f322      = -1.875 * sinim * (1.0 + 2.0 * cosim - 3.0 * cosisq);
      f441      = 35.0 * sini2 * f220;
      f442      = 39.3750 * sini2 * sini2;
      f522      = 9.84375 * sinim * (sini2 * (1.0 - 2.0 * cosim - 5.0 * cosisq)
                + 0.33333333 * (-2.0 + 4.0 * cosim + 6.0 * cosisq));
      f523      = sinim * (4.92187512 * sini2
                * (-2.0 - 4.0 * cosim + 10.0 * cosisq)
                + 6.56250012 * (1.0 + 2.0 * cosim - 3.0 * cosisq));
      f542      = 29.53125 * sinim * (2.0 - 8.0 * cosim
                + cosisq * (-12.0 + 8.0 * cosim + 10.0 * cosisq));
      f543      = 29.53125 * sinim * (-2.0 - 8.0 * cosim
                + cosisq * (12.0 + 8.0 * cosim - 10.0 * cosisq));
      xno2      = nm * nm;
      ainv2     = aonv * aonv;
      temp1     = 3.0 * xno2 * ainv2;
      temp      = temp1 * root22;
      sat.d2201 = temp * f220 * g201;
      sat.d2211 = temp * f221 * g211;
      temp1     = temp1 * aonv;
      temp      = temp1 * root32;
      sat.d3210 = temp * f321 * g310;
      sat.d3222 = temp * f322 * g322;
      temp1     = temp1 * aonv;
      temp      = 2.0 * temp1 * root44;
      sat.d4410 = temp * f441 * g410;
      sat.d4422 = temp * f442 * g422;
      temp1     = temp1 * aonv;
      temp      = temp1 * root52;
      sat.d5220 = temp * f522 * g520;
      sat.d5232 = temp * f523 * g532;
      temp      = 2.0 * temp1 * root54;
      sat.d5421 = temp * f542 * g521;
      sat.d5433 = temp * f543 * g533;
      sat.xlamo = fmod(sat.mo + sat.nodeo + sat.nodeo - theta - theta, kPi2);
      sat.xfact = sat.mdot + sat.dmdt + 2.0 * (sat.nodedot + sat.dnodt - rptim)
                - sat.no;
      em        = emo;
      emsq      = emsqo;

    // ---- synchronous resonance terms ----
    } else if (sat.irez == 1) {
      g200      = 1.0 + emsq * (-2.5 + 0.8125 * emsq);
      g310      = 1.0 + 2.0 * emsq;
      g300      = 1.0 + emsq * (-6.0 + 6.60937 * emsq);
      f220      = 0.75 * (1.0 + cosim) * (1.0 + cosim);
      f311      = 0.9375 * sinim * sinim * (1.0 + 3.0 * cosim)
                - 0.75 * (1.0 + cosim);
      f330      = 1.0 + cosim;
      f330      = 1.875 * f330 * f330 * f330;
      sat.del1  = 3.0 * nm * nm * aonv * aonv;
      sat.del2  = 2.0 * sat.del1 * f220 * g200 * q22;
      sat.del3  = 3.0 * sat.del1 * f330 * g300 * q33 * aonv;
      sat.del1 *= f311 * g310 * q31 * aonv;
      sat.xlamo = fmod(sat.mo + sat.nodeo + sat.argpo - theta, kPi2);
      sat.xfact = sat.mdot + xpidot - rptim + sat.dmdt + sat.domdt + sat.dnodt
                - sat.no;
    }

    // ---- for sgp4, initialize the integrator ----
    sat.xli   = sat.xlamo;
    sat.xni   = sat.no;
    sat.atime = 0.0;
    nm        = sat.no + dndt;
  } catch (...) {
    throw;
  }
}  // dsinit

/*
 * @brief       SGP4 prediction model
 *              * this procedure is the sgp4 prediction model from space command.
 *                this is an updated and combined version of sgp4 and sdp4, 
 *                which were originally published separately in spacetrack 
 *                report #3. this version follows the methodology from the aiaa 
 *                paper (2006) describing the history and development of the code.
 *
 * @param[in]   tsince (double)
 * @param[ref]  sat (Satellite )
 * @return      位置・速度 (PvTeme)
 */
PvTeme Sgp4::sgp4(double tsince, Satellite& sat) {
  int    ktr;
  double mrt;
  double temp;
  double temp4;
  double tempa;
  double tempe;
  double templ;
  double x2o3;
  double vkmpersec;
  double xmdf;
  double argpdf;
  double argpm;
  double am;
  double em;
  double mm;
  double nm;
  double nodedf;
  double nodem;
  double delomg;
  double delm;
  double delmtemp;
  double inclm;
  double t2;
  double t3;
  double t4;
  double tc;
  double dndt;
  double emsq;
  double xlm;
  double argpp;
  double ep;
  double mp;
  double nodep;
  double xincp;
  double sinim;
  double sinip;
  double cosim;
  double cosip;
  double axnl;
  double aynl;
  double xl;
  double eo1;
  double tem5;
  double u;
  double coseo1;
  double sineo1;
  double ecose;
  double esine;
  double el2;
  double pl;
  double betal;
  double rdotl;
  double rl;
  double rvdotl;
  double cosu;
  double cos2u;
  double sinu;
  double sin2u;
  double su;
  double temp1;
  double temp2;
  double cosisq;
  double mvt;
  double rvdot;
  double xinc;
  double xnode;
  double ux;
  double uy;
  double uz;
  double vx;
  double vy;
  double vz;
  double cnod;
  double cosi;
  double cossu;
  double sini;
  double sinsu;
  double snod;
  double xmx;
  double xmy;
  double mr;
  PvTeme teme = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

  try {
    sat.error = 0;
    mrt       = 0.0;

    // ---- set mathematical constants ----
    // sgp4fix divisor for divide by zero check on inclination
    // the old check used 1.0 + cos(pi-1.0e-9), but then compared it to
    // 1.5 e-12, so the threshold was changed to 1.5e-12 for consistency
    temp4 = 1.5e-12;
    x2o3  = 2.0 / 3.0;
    // sgp4fix identify constants and allow alternate values
    vkmpersec = cst.r * cst.xke / 60.0;

    // ---- clear sgp4 error flag ----
    sat.t = tsince;

    // ---- update for secular gravity and atmospheric drag ----
    xmdf   = sat.mo + sat.mdot * sat.t;
    argpdf = sat.argpo + sat.argpdot * sat.t;
    nodedf = sat.nodeo + sat.nodedot * sat.t;
    argpm  = argpdf;
    mm     = xmdf;
    t2     = sat.t * sat.t;
    nodem  = nodedf + sat.nodecf * t2;
    tempa  = 1.0 - sat.cc1 * sat.t;
    tempe  = sat.bstar * sat.cc4 * sat.t;
    templ  = sat.t2cof * t2;

    if (sat.isimp != 1) {
      delomg = sat.omgcof * sat.t;
      // sgp4fix use mutliply for speed instead of pow
      delmtemp = 1.0 + sat.eta * cos(xmdf);
      delm     = sat.xmcof * (delmtemp * delmtemp * delmtemp - sat.delmo);
      temp     = delomg + delm;
      mm       = xmdf + temp;
      argpm    = argpdf - temp;
      t3       = t2 * sat.t;
      t4       = t3 * sat.t;
      tempa   -= sat.d2 * t2 + sat.d3 * t3 + sat.d4 * t4;
      tempe   += sat.bstar * sat.cc5 * (sin(mm) - sat.sinmao);
      templ   += sat.t3cof * t3 + t4 * (sat.t4cof + sat.t * sat.t5cof);
    }

    nm    = sat.no;
    em    = sat.ecco;
    inclm = sat.inclo;
    if (sat.method == 'd') {
      tc = sat.t;
      dspace(tc, em, argpm, inclm, mm, nodem, nm, dndt, sat);
    }

    if (nm <= 0.0) {
      sat.error = 2;
      return teme;
    }

    // mean motion less than 0.0
    am  = pow(cst.xke / nm, x2o3) * tempa * tempa;
    nm  = cst.xke / pow(am, 1.5);
    em -= tempe;

    // fix tolerance for error recognition
    if (em >= 1.0 || em < -0.001 || am < 0.95) {
      sat.error = 1;
      return teme;
    }

    // sgp4fix fix tolerance to avoid a divide by zero
    if (em < 1.0e-6) { em = 1.0e-6; }
    mm  += sat.no * templ;
    xlm  = mm + argpm + nodem;
    emsq = em * em;
    temp = 1.0 - emsq;

    nodem = fmod(nodem, kPi2);
    argpm = fmod(argpm, kPi2);
    xlm   = fmod(xlm, kPi2);
    mm    = fmod(xlm - argpm - nodem, kPi2);
    while (mm < 0.0) { mm += kPi2; }

    // ---- compute extra mean quantities ----
    sinim = sin(inclm);
    cosim = cos(inclm);

    // ---- add lunar-solar periodics ----
    ep    = em;
    xincp = inclm;
    argpp = argpm;
    nodep = nodem;
    mp    = mm;
    sinip = sinim;
    cosip = cosim;
    if (sat.method == 'd') {
      dpper('n', ep, xincp, nodep, argpp, mp, sat);

      if (xincp < 0.0) {
        xincp *= -1.0;
        nodep += kPi;
        argpp -= kPi;
      }

      if (ep < 0.0 || ep > 1.0) {
        sat.error = 3;
        return teme;
      }
    }

    // ---- long period periodics ----
    if (sat.method == 'd') {
      sinip = sin(xincp);
      cosip = cos(xincp);
      sat.aycof = -0.5 * cst.j3oj2 * sinip;
      // sgp4fix for divide by zero for xincp = 180 deg
      if (fabs(cosip + 1.0) > 1.5e-12) {
        sat.xlcof = -0.25 * cst.j3oj2 * sinip * (3.0 + 5.0 * cosip)
                  / (1.0 + cosip);
      } else {
        sat.xlcof = -0.25 * cst.j3oj2 * sinip * (3.0 + 5.0 * cosip)
                  / temp4;
      }
    }

    axnl = ep * cos(argpp);
    temp = 1.0 / (am * (1.0 - ep * ep));
    aynl = ep * sin(argpp) + temp * sat.aycof;
    xl   = mp + argpp + nodep + temp * sat.xlcof * axnl;

    // ---- solve kepler's equation ----
    u    = fmod(xl - nodep, kPi2);
    eo1  = u;
    tem5 = 9999.9;
    ktr  = 1;
    // sgp4fix for kepler iteration
    // the following iteration needs better limits on corrections
    while (fabs(tem5) >= 1.0e-12 && ktr <= 10) {
      sineo1 = sin(eo1);
      coseo1 = cos(eo1);
      tem5   = 1.0 - coseo1 * axnl - sineo1 * aynl;
      tem5   = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5;
      if (fabs(tem5) >= 0.95) {
        if (tem5 > 0.0) {
          tem5 = 0.95;
        } else {
          tem5 = -0.95;
        }
      }
      eo1 += tem5;
      ++ktr;
    }

    // ---- short period preliminary quantities ----
    ecose = axnl * coseo1 + aynl * sineo1;
    esine = axnl * sineo1 - aynl * coseo1;
    el2   = axnl * axnl + aynl * aynl;
    pl    = am * (1.0 - el2);
    if (pl < 0.0) {
      sat.error = 4;
      return teme;
    } else {
      rl     = am * (1.0 - ecose);
      rdotl  = sqrt(am) * esine/rl;
      rvdotl = sqrt(pl) / rl;
      betal  = sqrt(1.0 - el2);
      temp   = esine / (1.0 + betal);
      sinu   = am / rl * (sineo1 - aynl - axnl * temp);
      cosu   = am / rl * (coseo1 - axnl + aynl * temp);
      su     = atan2(sinu, cosu);
      sin2u  = (cosu + cosu) * sinu;
      cos2u  = 1.0 - 2.0 * sinu * sinu;
      temp   = 1.0 / pl;
      temp1  = 0.5 * cst.j2 * temp;
      temp2  = temp1 * temp;

      // ---- update for short period periodics ----
      if (sat.method == 'd') {
        cosisq = cosip * cosip;
        sat.con41  = 3.0 * cosisq - 1.0;
        sat.x1mth2 = 1.0 - cosisq;
        sat.x7thm1 = 7.0 * cosisq - 1.0;
      }

      mrt   = rl * (1.0 - 1.5 * temp2 * betal * sat.con41)
            + 0.5 * temp1 * sat.x1mth2 * cos2u;
      su   -= 0.25 * temp2 * sat.x7thm1 * sin2u;
      xnode = nodep + 1.5 * temp2 * cosip * sin2u;
      xinc  = xincp + 1.5 * temp2 * cosip * sinip * cos2u;
      mvt   = rdotl - nm * temp1 * sat.x1mth2 * sin2u / cst.xke;
      rvdot = rvdotl + nm * temp1 * (sat.x1mth2 * cos2u
            + 1.5 * sat.con41) / cst.xke;

      // ---- orientation vectors ----
      sinsu =  sin(su);
      cossu =  cos(su);
      snod  =  sin(xnode);
      cnod  =  cos(xnode);
      sini  =  sin(xinc);
      cosi  =  cos(xinc);
      xmx   = -snod * cosi;
      xmy   =  cnod * cosi;
      ux    =  xmx * sinsu + cnod * cossu;
      uy    =  xmy * sinsu + snod * cossu;
      uz    =  sini * sinsu;
      vx    =  xmx * cossu - cnod * sinsu;
      vy    =  xmy * cossu - snod * sinsu;
      vz    =  sini * cossu;

      // ---- position and velocity (in km and km/sec) ----
      mr = mrt * cst.r;
      teme.r.x = mr * ux;
      teme.r.y = mr * uy;
      teme.r.z = mr * uz;
      teme.v.x = (mvt * ux + rvdot * vx) * vkmpersec;
      teme.v.y = (mvt * uy + rvdot * vy) * vkmpersec;
      teme.v.z = (mvt * uz + rvdot * vz) * vkmpersec;
    }

    // sgp4fix for decaying satellites
    if (mrt < 1.0) { sat.error = 6; }
  } catch (...) {
    throw;
  }

  return teme;
}  // sgp4

/*
 * @brief       Deep space contributions
 *              * this procedure provides deep space contributions to mean 
 *                elements for perturbing third body.  these effects have been 
 *                averaged over one revolution of the sun and moon.  for earth 
 *                resonance effects, the effects have been averaged over no 
 *                revolutions of the satellite.  (mean motion)
 *
 * @param[in]   tc (double)
 * @param[ref]  em, argpm, inclm, mm, nodem, nm, dndt (double)
 * @param[ref]  sat (Satellite)
 */
void Sgp4::dspace(
    double tc,
    double& em, double& argpm, double& inclm, double& mm,
    double& nodem, double& nm, double& dndt,
    Satellite& sat) {
  //int    iret;
  int    iretn;
  double fasx2;
  double fasx4;
  double fasx6;
  double g22;
  double g32;
  double g44;
  double g52;
  double g54;
  double rptim;
  double stepp;
  double stepn;
  double step2;
  double theta;
  double ft;
  double delt;
  double x2li;
  double xomi;
  double x2omi;
  double xl;
  double xldot;
  double xnddt;
  double xndt;

  try {
    fasx2 = 0.13130908;
    fasx4 = 2.8843198;
    fasx6 = 0.37448087;
    g22   = 5.7686396;
    g32   = 0.95240898;
    g44   = 1.8014998;
    g52   = 1.0508330;
    g54   = 4.4108898;
    rptim = 4.37526908801129966e-3;  // equates to 7.29211514668855e-5 rad/sec
    stepp =    720.0;
    stepn =   -720.0;
    step2 = 259200.0;

    // ---- calculate deep space resonance effects ----
    dndt   = 0.0;
    theta  = fmod(sat.gsto + tc * rptim, kPi2);
    em    += sat.dedt * sat.t;
    inclm += sat.didt * sat.t;
    argpm += sat.domdt * sat.t;
    nodem += sat.dnodt * sat.t;
    mm    += sat.dmdt * sat.t;

    // ---- epoch restart ----
    // sgp4fix for propagator problems
    // the following integration works for negative time steps and periods
    // the specific changes are unknown because the original code was so convoluted
    //
    // sgp4fix take out atime = 0.0 and fix for faster operation
    ft = 0.0;
    if (sat.irez != 0) {
      // sgp4fix streamline check
      if (sat.atime == 0.0 || sat.t * sat.atime <= 0.0 ||
          abs(sat.t) < abs(sat.atime)) {
        sat.atime = 0.0;
        sat.xni   = sat.no;
        sat.xli   = sat.xlamo;
      }

      // sgp4fix move check outside loop
      delt = stepn;
      if (sat.t > 0.0) { delt = stepp; }

      iretn = 381;  // added for do loop
      while (iretn == 381) {
        // ---- dot terms calculated ------------------
        // ---- near - synchronous resonance terms ----
        if (sat.irez != 2) {
          xndt  = sat.del1 * sin(sat.xli - fasx2)
                + sat.del2 * sin(2.0 * (sat.xli - fasx4))
                + sat.del3 * sin(3.0 * (sat.xli - fasx6));
          xldot = sat.xni + sat.xfact;
          xnddt = sat.del1 * cos(sat.xli - fasx2)
                + 2.0 * sat.del2 * cos(2.0 * (sat.xli - fasx4))
                + 3.0 * sat.del3 * cos(3.0 * (sat.xli - fasx6));
          xnddt *= xldot;
        } else {
          // ---- near - half-day resonance terms ----
          xomi  = sat.argpo + sat.argpdot * sat.atime;
          x2omi = xomi + xomi;
          x2li  = sat.xli + sat.xli;
          xndt  = (sat.d2201 * sin(x2omi + sat.xli - g22)
                +  sat.d2211 * sin(        sat.xli - g22)
                +  sat.d3210 * sin( xomi + sat.xli - g32)
                +  sat.d3222 * sin(-xomi + sat.xli - g32)
                +  sat.d4410 * sin(x2omi + x2li    - g44)
                +  sat.d4422 * sin(        x2li    - g44)
                +  sat.d5220 * sin( xomi + sat.xli - g52)
                +  sat.d5232 * sin(-xomi + sat.xli - g52)
                +  sat.d5421 * sin( xomi + x2li    - g54)
                +  sat.d5433 * sin(-xomi + x2li    - g54));
          xldot = sat.xni + sat.xfact;
          xnddt = (sat.d2201 * cos(x2omi + sat.xli - g22)
                +  sat.d2211 * cos(        sat.xli - g22)
                +  sat.d3210 * cos( xomi + sat.xli - g32)
                +  sat.d3222 * cos(-xomi + sat.xli - g32)
                +  sat.d5220 * cos( xomi + sat.xli - g52)
                +  sat.d5232 * cos(-xomi + sat.xli - g52)
                +  2.0
                * (sat.d4410 * cos(x2omi + x2li - g44)
                +  sat.d4422 * cos(        x2li - g44)
                +  sat.d5421 * cos( xomi + x2li - g54)
                +  sat.d5433 * cos(-xomi + x2li - g54)));
          xnddt *= xldot;
        }

        // ---- integrator ----
        // sgp4fix move end checks to end of routine
        if (fabs(sat.t - sat.atime) >= stepp) {
          iretn = 381;
        } else {
          ft    = sat.t - sat.atime;
          iretn = 0;
        }
        if (iretn == 381) {
          sat.xli   += xldot * delt + xndt * step2;
          sat.xni   += xndt * delt + xnddt * step2;
          sat.atime += delt;
        }
      }

      nm = sat.xni + xndt  * ft + xnddt * ft * ft * 0.5;
      xl = sat.xli + xldot * ft + xndt  * ft * ft * 0.5;
      if (sat.irez != 1) {
        mm   = xl - 2.0 * nodem + 2.0 * theta;
        dndt = nm - sat.no;
      } else {
        mm   = xl - nodem - argpm + theta;
        dndt = nm - sat.no;
      }

      nm = sat.no + dndt;
    }
  } catch (...) {
    throw;
  }
} // dspace

/*
 * @brief       Deep space long period periodic contributions
 *              * this procedure provides deep space long period periodic 
 *                contributions to the mean elements.  by design, these periodics 
 *                are zero at epoch.  This used to be dscom which included 
 *                initialization, but it's really a recurring function.
 *
 * @param[in]   opsmode (char)
 * @param[ref]  ep, xincp, nodep, argpp, mp, sat (double)
 * @param[ref]  sat (Satellite)
 */
void Sgp4::dpper(
    char opsmode,
    double& ep, double& xincp, double& nodep, double& argpp, double& mp,
    Satellite& sat) {
  double e3;
  double ee2;
  double peo;
  double pgho;
  double pho;
  double pinco;
  double plo;
  double se2;
  double se3;
  double sgh2;
  double sgh3;
  double sgh4;
  double sh2;
  double sh3;
  double si2;
  double si3;
  double sl2;
  double sl3;
  double sl4;
  double t;
  double xgh2;
  double xgh3;
  double xgh4;
  double xh2;
  double xh3;
  double xi2;
  double xi3;
  double xl2;
  double xl3;
  double xl4;
  double xls;
  double xnoh;
  double zmol;
  double zmos;
  double zns;
  double zes;
  double znl;
  double zel;
  double zm;
  double f2;
  double f3;
  double sel;
  double ses;
  double sghl;
  double sghs;
  double shll;
  double shs;
  double sil;
  double sinzf;
  double sis;
  double sll;
  double sls;
  double pe;
  double pgh;
  double ph;
  double pinc;
  double pl;
  double zf;
  double cosip;
  double cosop;
  double sinip;
  double sinop;
  double inclp;
  double alfdp;
  double betdp;
  double dalf;
  double dbet;

  try {
    // Copy satellite attributes into local variables for convenience
    // and symmetry in writing formulae.
    e3    = sat.e3;
    ee2   = sat.ee2;
    peo   = sat.peo;
    pgho  = sat.pgho;
    pho   = sat.pho;
    pinco = sat.pinco;
    plo   = sat.plo;
    se2   = sat.se2;
    se3   = sat.se3;
    sgh2  = sat.sgh2;
    sgh3  = sat.sgh3;
    sgh4  = sat.sgh4;
    sh2   = sat.sh2;
    sh3   = sat.sh3;
    si2   = sat.si2;
    si3   = sat.si3;
    sl2   = sat.sl2;
    sl3   = sat.sl3;
    sl4   = sat.sl4;
    t     = sat.t;
    xgh2  = sat.xgh2;
    xgh3  = sat.xgh3;
    xgh4  = sat.xgh4;
    xh2   = sat.xh2;
    xh3   = sat.xh3;
    xi2   = sat.xi2;
    xi3   = sat.xi3;
    xl2   = sat.xl2;
    xl3   = sat.xl3;
    xl4   = sat.xl4;
    zmol  = sat.zmol;
    zmos  = sat.zmos;

    // ---------------------- constants -----------------------------
    zns = 1.19459e-5;
    zes = 0.01675;
    znl = 1.5835218e-4;
    zel = 0.05490;

    // --------------- calculate time varying periodics -----------
    zm    = zmos + zns * t;
    // be sure that the initial call has time set to zero
    if (sat.init == 'y') { zm = zmos; }
    zf    = zm + 2.0 * zes * sin(zm);
    sinzf = sin(zf);
    f2    =  0.5 * sinzf * sinzf - 0.25;
    f3    = -0.5 * sinzf * cos(zf);
    ses   = se2 * f2 + se3 * f3;
    sis   = si2 * f2 + si3 * f3;
    sls   = sl2 * f2 + sl3 * f3 + sl4 * sinzf;
    sghs  = sgh2 * f2 + sgh3 * f3 + sgh4 * sinzf;
    shs   = sh2 * f2 + sh3 * f3;
    zm    = zmol + znl * t;
    if (sat.init == 'y') { zm = zmol; }
    zf    = zm + 2.0 * zel * sin(zm);
    sinzf = sin(zf);
    f2    =  0.5 * sinzf * sinzf - 0.25;
    f3    = -0.5 * sinzf * cos(zf);
    sel   = ee2 * f2 + e3  * f3;
    sil   = xi2 * f2 + xi3 * f3;
    sll   = xl2 * f2 + xl3 * f3 + xl4 * sinzf;
    sghl  = xgh2 * f2 + xgh3 * f3 + xgh4 * sinzf;
    shll  = xh2 * f2 + xh3 * f3;
    pe    = ses + sel;
    pinc  = sis + sil;
    pl    = sls + sll;
    pgh   = sghs + sghl;
    ph    = shs + shll;

    if (sat.init == 'n') {
      pe    -= peo;
      pinc  -= pinco;
      pl    -= plo;
      pgh   -= pgho;
      ph    -= pho;
      inclp  = 0.0;  // inclp += pinc;
      ep    += pe;
      sinip  = sin(inclp);
      cosip  = cos(inclp);

      // ---- apply periodics directly ----
      // sgp4fix for lyddane choice
      // strn3 used original inclination - this is technically feasible
      // gsfc used perturbed inclination - also technically feasible
      // probably best to readjust the 0.2 limit value and limit discontinuity
      // 0.2 rad = 11.45916 deg
      // use next line for original strn3 approach and original inclination
      // if (inclo >= 0.2)
      // use next line for gsfc version and perturbed inclination
      if (inclp >= 0.2) {
        ph    /= sinip;
        pgh   -= cosip * ph;
        argpp += pgh;
        nodep += ph;
        mp    += pl;
      } else {
        // ---- apply periodics with lyddane modification ----
        sinop  = sin(nodep);
        cosop  = cos(nodep);
        alfdp  = sinip * sinop;
        betdp  = sinip * cosop;
        dalf   =  ph * cosop + pinc * cosip * sinop;
        dbet   = -ph * sinop + pinc * cosip * cosop;
        alfdp += dalf;
        betdp += dbet;
        nodep  = fmod(nodep, kPi2);
        // sgp4fix for afspc written intrinsic functions
        // nodep used without a trigonometric function ahead
        if (nodep < 0.0 && opsmode == 'a') { nodep += kPi2; }
        xls = mp + argpp + pl + pgh + (cosip - pinc * sinip) * nodep;
        xnoh  = nodep;
        nodep = atan2(alfdp, betdp);
        // sgp4fix for afspc written intrinsic functions
        // nodep used without a trigonometric function ahead
        if (nodep < 0.0 && opsmode == 'a') { nodep += kPi2; }
        if (fabs(xnoh - nodep) > kPi) {
          if (nodep < xnoh) {
            nodep += kPi2;
          } else {
            nodep -= kPi2;
          }
        }
        mp += pl;
        argpp = xls - mp - cosip * nodep;
      }
    }
  } catch (...) {
    throw;
  }
}  // dpper

}  // namespace iss_sgp4_teme

