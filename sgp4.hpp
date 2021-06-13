#ifndef ISS_SGP4_TEME_SGP4_HPP_
#define ISS_SGP4_TEME_SGP4_HPP_

#include "time.hpp"
#include "tle.hpp"

#include <iomanip>
#include <string>
#include <vector>

namespace iss_sgp4_teme {

// 定数セット構造体
struct Const {
  double r;
  double mu;
  double xke;
  double tumin;
  double j2;
  double j3;
  double j4;
  double j3oj2;
};
// 衛星情報構造体
struct Satellite {
  // Unique satellite number given in the TLE file.
  int    satnum     = 0;
  // Full four-digit year of this element set's epoch moment.
  int    epochyr    = 0;
  // Fractional days into the year of the epoch moment.
  double epochdays  = 0.0;
  // Julian date of the epoch (computed from epochyr and epochdays).
  double jdsatepoch = 0.0;
  // First time derivative of the mean motion (ignored by SGP4).
  double ndot       = 0.0;
  // Second time derivative of the mean motion (ignored by SGP4).
  double nddot      = 0.0;
  // Ballistic drag coefficient B* in inverse earth radii.
  double bstar      = 0.0;
  // Inclination in radians.
  double inclo      = 0.0;
  // Right ascension of ascending node in radians.
  double nodeo      = 0.0;
  // Eccentricity.
  double ecco       = 0.0;
  // Argument of perigee in radians.
  double argpo      = 0.0;
  // Mean anomaly in radians.
  double mo         = 0.0;
  // Mean motion in radians per minute.
  double no         = 0.0;
  //
  // Near Earth
  char   method  = 'n';
  char   opsmode = 'i';
  char   init    = 'y';
  int    isimp   = 0;
  double aycof   = 0.0;
  double con41   = 0.0;
  double cc1     = 0.0;
  double cc4     = 0.0;
  double cc5     = 0.0;
  double d2      = 0.0;
  double d3      = 0.0;
  double d4      = 0.0;
  double delmo   = 0.0;
  double eta     = 0.0;
  double argpdot = 0.0;
  double omgcof  = 0.0;
  double sinmao  = 0.0;
  double t       = 0.0;
  double t2cof   = 0.0;
  double t3cof   = 0.0;
  double t4cof   = 0.0;
  double t5cof   = 0.0;
  double x1mth2  = 0.0;
  double x7thm1  = 0.0;
  double mdot    = 0.0;
  double nodedot = 0.0;
  double xlcof   = 0.0;
  double xmcof   = 0.0;
  double nodecf  = 0.0;
  //
  // Deep space
  int    irez  = 0;
  double d2201 = 0.0;
  double d2211 = 0.0;
  double d3210 = 0.0;
  double d3222 = 0.0;
  double d4410 = 0.0;
  double d4422 = 0.0;
  double d5220 = 0.0;
  double d5232 = 0.0;
  double d5421 = 0.0;
  double d5433 = 0.0;
  double dedt  = 0.0;
  double del1  = 0.0;
  double del2  = 0.0;
  double del3  = 0.0;
  double didt  = 0.0;
  double dmdt  = 0.0;
  double dnodt = 0.0;
  double domdt = 0.0;
  double e3    = 0.0;
  double ee2   = 0.0;
  double peo   = 0.0;
  double pgho  = 0.0;
  double pho   = 0.0;
  double pinco = 0.0;
  double plo   = 0.0;
  double se2   = 0.0;
  double se3   = 0.0;
  double sgh2  = 0.0;
  double sgh3  = 0.0;
  double sgh4  = 0.0;
  double sh2   = 0.0;
  double sh3   = 0.0;
  double si2   = 0.0;
  double si3   = 0.0;
  double sl2   = 0.0;
  double sl3   = 0.0;
  double sl4   = 0.0;
  double gsto  = 0.0;
  double xfact = 0.0;
  double xgh2  = 0.0;
  double xgh3  = 0.0;
  double xgh4  = 0.0;
  double xh2   = 0.0;
  double xh3   = 0.0;
  double xi2   = 0.0;
  double xi3   = 0.0;
  double xl2   = 0.0;
  double xl3   = 0.0;
  double xl4   = 0.0;
  double xlamo = 0.0;
  double zmol  = 0.0;
  double zmos  = 0.0;
  double atime = 0.0;
  double xli   = 0.0;
  double xni   = 0.0;
  // error
  int    error = 0;
};
// 座標構造体
struct Coord{
  double x;
  double y;
  double z;
};
// 位置・速度構造体(TEME)
struct PvTeme {
  Coord r;  // 位置
  Coord v;  // 速度
};

class Sgp4 {
  Const cst;  // 定数(gravconst)

public:
  Sgp4(struct timespec, std::vector<std::string>, std::string = "wgs84");
                                                  // コンストラクタ
  Satellite twoline2rv(bool afspc_mode = false);  // ISS 初期位置・速度の取得
  PvTeme propagate(Satellite&);                   // 指定 UT1 の ISS 位置・速度の取得

private:
  struct timespec ut1;               // UT1
  std::vector<std::string> tle;      // TLE
  Const get_gravconst(std::string);  // 定数取得
  void sgp4init(Satellite&);         // SGP4 初期化
  void initl(
      double,
      double&, double&, double&, double&, double&, double&,
      double&, double&, double&, double&, double&, double&, Satellite&);
                                     // SGP4 propagator 初期化
  void dscom(
      double,  double,
      double&, double&, double&, double&, double&, double&, double&, double&,
      double&, double&, double&, double&, double&, double&, double&, double&,
      double&, double&, double&, double&, double&, double&, double&, double&,
      double&, double&, double&, double&, double&, double&, double&, double&,
      double&, double&, double&, double&, double&, double&, double&, double&,
      double&, double&, double&, double&, double&, double&, double&, double&,
      double&, double&,
      Satellite&);                   // Deep space common items
  void dsinit(
      double,  double,  double,  double,  double,  double,  double,  double,
      double,  double,  double,  double,  double,  double,  double,  double,
      double,  double,  double,  double,  double,  double,  double,  double,
      double,  double,  double,  double,  double,  double,  double,  double,
      double&, double&, double&, double&, double&, double&, double&,
      Satellite&);                   // Deep space contributions 初期化
  PvTeme sgp4(double, Satellite&);   // SGP4 prediction model
  void dspace(
      double,
      double&, double&, double&, double&, double&, double&, double&,
      Satellite&);                   // Deep space contributions
  void dpper(
      char,
      double&, double&, double&, double&, double&,
      Satellite&);                   // Deep space long period periodic contributions

};

}  // namespace iss_sgp4_teme

#endif

