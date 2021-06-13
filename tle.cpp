#include "tle.hpp"

namespace iss_sgp4_teme {

// 定数
static constexpr char         kFTle[]  = "tle.txt";
static constexpr unsigned int kSecDay  =  86400;  // Seconds in a day

/*
 * @brief      コンストラクタ
 *
 * @param[in]  UT1 (timespec)
 */
Tle::Tle(struct timespec ut1) {
  this->ut1 = ut1;
}

/*
 * @brief   TLE 読み込み
 *
 * @param   <none>
 * @return  TLE(2行) (vector<string>)
 */
std::vector<std::string> Tle::get_tle() {
  std::string              f(kFTle);  // ファイル名
  std::vector<std::string> data;      // TLE 一覧（全データ）
  struct tm                t;         // 時刻構造体
  std::string              buf;       // 1行分バッファ
  std::vector<std::string> tle_p(2);  // TLE（退避用）
  unsigned int             y;         // year
  double                   d;         // day
  struct timespec          utc;       // UTC
  unsigned int             i;         // loop index
  std::vector<std::string> tle(2, "");  // TLE

  try {
    // ファイル OPEN
    std::ifstream ifs(f);
    if (!ifs) throw;  // 読み込み失敗

    // ファイル READ（一旦、全て読み込む）
    while (getline(ifs, buf)) {
      std::istringstream iss(buf);   // 文字列ストリーム
      if (buf[0] != '1' && buf[0] != '2') { continue; }
      data.push_back(buf);
    }

    // 最新 TLE 検索
    for (auto l: data) {
      if (l[0] == '1') {
        y = 2000 + stoi(l.substr(18, 2));
        d = stod(l.substr(20, 12));
        // y 年の 01-01 00:00:00
        std::stringstream ss;
        ss << std::setfill('0')
           << std::setw(4) << y
           << std::setw(2) << 1
           << std::setw(2) << 1
           << std::setw(2) << 0
           << std::setw(2) << 0
           << std::setw(2) << 0;
        std::istringstream is(ss.str());
        is >> std::get_time(&t, "%Y%m%d%H%M%S");
        utc.tv_sec = mktime(&t);
        utc.tv_nsec = 0;
        // utc = y 年の 01-01 00:00:00 に経過日数 d を加算した日時
        utc = ts_add(utc, d * kSecDay);
        if (utc.tv_sec > ut1.tv_sec) {
          tle = tle_p;
          break;
        }
      }
      tle_p[stoi(l.substr(0, 1)) - 1] = l;
    }

    // 上記の処理で該当レコードが得られなかった場合は、最初の2行
    if (tle[0] == "") {
      for (i = 0; i < 2; ++i) { tle[i] = data[i]; }
    }
  } catch (...) {
    throw;
  }

  return tle;
}

}  // namespace iss_sgp4_teme

