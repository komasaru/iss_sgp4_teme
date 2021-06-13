/***********************************************************
  TLE から ISS の位置・速度を計算
  : 但し、座標系は TEME(True Equator, Mean Equinox; 真赤道面平均春分点)
  : 指定 JST における BLH 座標の計算は iss_sgp4_blh

    DATE        AUTHOR       VERSION
    2021.05.13  mk-mode.com  1.00 新規作成

  Copyright(C) 2021 mk-mode.com All Rights Reserved.

  引数 : UT1（世界時1）
           書式：最大23桁の数字
                 （先頭から、西暦年(4), 月(2), 日(2), 時(2), 分(2), 秒(2),
                             1秒未満(9)（小数点以下9桁（ナノ秒）まで））
                 無指定なら現在(システム日時)と判断。
***********************************************************/
#include "sgp4.hpp"
#include "time.hpp"
#include "tle.hpp"

#include <cstdlib>   // for EXIT_XXXX
#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char* argv[]) {
  namespace ns = iss_sgp4_teme;
  std::string tm_str;            // time string
  unsigned int s_tm;             // size of time string
  int s_nsec;                    // size of nsec string
  int ret;                       // return of functions
  struct tm t = {};              // for work
  struct timespec ut1;           // UT1
  std::vector<std::string> tle;  // TLE
  ns::Satellite sat;             // 衛星情報
  ns::PvTeme    teme;            // 位置・速度

  try {
    // 現在日時(UT1) 取得
    if (argc > 1) {
      // コマンドライン引数より取得
      tm_str = argv[1];
      s_tm = tm_str.size();
      if (s_tm > 23) { 
        std::cout << "[ERROR] Over 23-digits!" << std::endl;
        return EXIT_FAILURE;
      }
      s_nsec = s_tm - 14;
      std::istringstream is(tm_str);
      is >> std::get_time(&t, "%Y%m%d%H%M%S");
      ut1.tv_sec  = mktime(&t);
      ut1.tv_nsec = 0;
      if (s_tm > 14) {
        ut1.tv_nsec = std::stod(
            tm_str.substr(14, s_nsec) + std::string(9 - s_nsec, '0'));
      }
    } else {
      // 現在日時の取得
      ret = std::timespec_get(&ut1, TIME_UTC);
      if (ret != 1) {
        std::cout << "[ERROR] Could not get now time!" << std::endl;
        return EXIT_FAILURE;
      }
    }

    // TLE 読み込み, gravconst 取得
    ns::Tle o_t(ut1);
    tle = o_t.get_tle();

    // ISS 初期位置・速度の取得
    ns::Sgp4 o_s(ut1, tle);
    sat = o_s.twoline2rv();

    // 指定 UT1 の ISS 位置・速度の取得
    teme = o_s.propagate(sat);

    // Calculation & display
    std::cout << "[" << ns::gen_time_str(ut1) << " UT1]" << std::endl;
    std::cout << "TLE:" << tle[0] << std::endl;
    std::cout << "    " << tle[1] << std::endl;
    std::cout << std::fixed << std::setprecision(8);
    std::cout << "TEME POS:["
              << std::setw(16) << teme.r.x << ", "
              << std::setw(16) << teme.r.y << ", "
              << std::setw(16) << teme.r.z
              << "]" << std::endl;
    std::cout << "     VEL:["
              << std::setw(16) << teme.v.x << ", "
              << std::setw(16) << teme.v.y << ", "
              << std::setw(16) << teme.v.z
              << "]" << std::endl;
  } catch (...) {
      std::cerr << "EXCEPTION!" << std::endl;
      return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

