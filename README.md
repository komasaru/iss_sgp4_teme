# iss_sgp4_teme

概要
====

* TLE から ISS の位置・速度を計算  
  （但し、座標系は TEME(True Equator, Mean Equinox; 真赤道面平均春分点)）  
  + 指定 JST における BLH 座標の計算は `iss_sgp4_blh`

ビルド方法
==========

`make`

（やり直す場合は、 `make clean` をしてから）

準備
====

* TLE（Two-line elements; 2行軌道要素形式）データを [こちら](https://celestrak.com/NORAD/elements/supplemental/iss.txt") からダウンロード。実行プログラムと同じディレクトリ内に配置する。（ファイル名: `tle.txt`）
* 当 GitHub リポジトリにアップされている `tle.txt` は古いものなので、注意。

実行方法
========

`./iss_sgp4_teme [YYYYMMDDHHMMSSMMMMMMMMM]`

* コマンドライン引数には UT1（世界時1） を指定する。
* UT1（世界時1）は「年・月・日・時・分・秒・ナノ秒」を最大23桁で指定する。
* UT1（世界時1）を指定しない場合は、システム日時を UT1 とみなす。
* UT1（世界時1）を先頭から部分的に指定した場合は、指定していない部分を 0 とみなす。

