#ifndef ISS_SGP4_TEME_TIME_HPP_
#define ISS_SGP4_TEME_TIME_HPP_

#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace iss_sgp4_teme {

struct DateTime {
  unsigned int year;
  unsigned int month;
  unsigned int day;
  unsigned int hour;
  unsigned int minute;
  double       second ;
};

std::string gen_time_str(struct timespec ts);     // 日時文字列生成
struct timespec ts_add(struct timespec, double);  // 時刻 + 秒数
DateTime days2ymdhms(unsigned int, double);       // 年+経過日数 => 年月日時分秒
double jday(DateTime);                            // 年月日時分秒 => ユリウス日
double gstime(double);                            // Greenwich sidereal time calculation

}  // namespace iss_sgp4_teme

#endif

