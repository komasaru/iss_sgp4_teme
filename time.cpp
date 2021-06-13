#include "time.hpp"

namespace iss_sgp4_teme {

static constexpr double kE9      = 1.0e9;
static constexpr double kEm6     = 1.0e-6;
static constexpr double kPi      = atan(1.0) * 4.0;      // 円周率
static constexpr double kPi2     = kPi * 2.0;            // 円周率 * 2
static constexpr double kDeg2Rad = kPi / 180.0;          // 0.0174532925199433

/*
 * @brief      日時文字列生成
 *
 * @param[in]  日時 (timespec)
 * @return     日時文字列 (string)
 */
std::string gen_time_str(struct timespec ts) {
  struct tm t;
  std::stringstream ss;
  std::string str_tm;

  try {
    localtime_r(&ts.tv_sec, &t);
    ss << std::setfill('0')
       << std::setw(4) << t.tm_year + 1900 << "-"
       << std::setw(2) << t.tm_mon + 1     << "-"
       << std::setw(2) << t.tm_mday        << " "
       << std::setw(2) << t.tm_hour        << ":"
       << std::setw(2) << t.tm_min         << ":"
       << std::setw(2) << t.tm_sec         << "."
       << std::setw(3) << round(ts.tv_nsec * kEm6);
    return ss.str();
  } catch (...) {
    throw;
  }
}

/*
 * @brief      時刻 + 秒数
 *
 * @param[in]  時刻 (timespec)
 * @param[in]  秒数 (double)
 * @return     時刻 (timespec)
 */
struct timespec ts_add(struct timespec ts_src, double s) {
  struct timespec ts;

  try {
    ts.tv_sec  = ts_src.tv_sec + int(s);
    ts.tv_nsec = ts_src.tv_nsec + (s - int(s)) * kE9;
    while (ts.tv_nsec > kE9) {
      ++ts.tv_sec;
      ts.tv_nsec -= kE9;
    }
    while (ts.tv_nsec < 0) {
      --ts.tv_sec;
      ts.tv_nsec += kE9;
    }
  } catch (...) {
    throw;
  }

  return ts;
}

/*
 * @brief       年+経過日数 => 年月日時分秒
 *              * this procedure converts the day of the year, days, to the 
 *                equivalent month, day, hour, minute and second.
 *              * algorithm: set up array for the number of days per month 
 *                           find leap year - use 1900 because 2000 is a leap 
 *                           year loop through a temp value while the value 
 *                           is < the days perform int conversions to the 
 *                           correct day and month convert remainder into 
 *                           h m s using type conversions
 *
 * @param[in]  年 (unsigned int)
 * @param[in]  日数 (double)
 * @return     日時 (DateTime)
 */
DateTime days2ymdhms(unsigned int year, double days) {
  unsigned int lmonth[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  unsigned int dayofyr;
  unsigned int i;
  unsigned int inttemp;
  double       temp;
  DateTime     dt;

  try {
    dayofyr = int(days);

    // ----------------- find month and day of month ----------------
    if (year % 4 == 0) { lmonth[1] = 29; }
    i = 0;
    inttemp = 0;
    while (dayofyr > inttemp + lmonth[i] && i < 12) {
      inttemp += lmonth[i];
      ++i;
    }
    dt.year  = year;
    dt.month = i + 1;
    dt.day   = dayofyr - inttemp;

    // ----------------- find hours minutes and seconds -------------
    temp   = (days - dayofyr) * 24.0;
    dt.hour   = int(temp);
    temp   = (temp - dt.hour) * 60.0;
    dt.minute = int(temp);
    dt.second = (temp - dt.minute) * 60.0;
  } catch (...) {
    throw;
  }

  return dt;
}

/*
 * @brief      年月日時分秒 => ユリウス日
 *             * this procedure finds the julian date given the year, month, day, and time.
 *               the julian date is defined by each elapsed day since noon, jan 1, 4713 bc.
 *             * algorithm: calculate the answer in one step for efficiency
 *
 * @param[in]  日時 (DateTime)
 * @return     ユリウス日 (double)
 */
double jday(DateTime dt) {
  double jd;

  try {
    jd = (367.0 * dt.year
       - int(7.0 * (dt.year + int((dt.month + 9.0) / 12.0)) * 0.25)
       + int(275.0 * dt.month / 9.0) + dt.day + 1721013.5
       + ((dt.second / 60.0 + dt.minute) / 60.0 + dt.hour) / 24.0);
  } catch (...) {
    throw;
  }

  return jd;
}

/*
 * @brief      Greenwich sidereal time calculation
 *
 * @param[in]  ユリウス日(UT1) (double)
 * @return     GST (double
 */
double gstime(double jdut1) {
  double tut1;
  double gst;

  try {
    tut1 = (jdut1 - 2451545.0) / 36525.0;
    gst = 67310.54841
        + ((876600.0 * 3600.0 + 8640184.812866)
        + (0.093104
        - 6.2e-6
        * tut1) * tut1) * tut1;
    gst = fmod(gst * kDeg2Rad / 240.0, kPi2);
    if (gst < 0.0) { gst += kPi2; }
  } catch (...) {
    throw;
  }

  return gst;
}

}  // namespace iss_sgp4_teme

