2022.3.8 水本
(mizumoto.misaki.n68@kyoto-u.jp) <- 半永久的に使えるアドレスです

fits
figures
README.txt
という構成になっています。
fits以下で作業をする想定で、figures以下に画像が出力されます。
今figures以下には実際に以下のスクリプトを走らせた場合の画像が入っています。

以下、current directory = fits と置き換えてください。


# 準備
current directory以下に
spec/spec_all.txt
という名前の2列49800行のascii fileがある。
1列目はenergy, 2列目はcountを示している。100eVから25000eVまで0.5eV区切りとなっている。
今は全ピクセルを積分したスペクトルが入っている。

current directory以下にCalLinesというディレクトリが存在する。
元素ごと、Ka or Kbごとにenergy, FWHM, relative amplitudeが書かれているファイル（＝ラインリスト）と、recommended*.txtというファイルがある。後者の方はどのラインリストを使うかを指定するためのものである。

# ラインフィット
今、spec/spec_all.txtをフィットすることを考える。
python3 linefitting_all.pyを実行すると、
./CalLines/recommended_TC4use.txt
に記載されているラインを順番にフィットする。
Lorentzianで自然幅の分だけ広がっているラインを、検出器がGaussianでさらに広げている。
フィット結果はlinefitting_result_all.txtに書かれる。
python3 linefitting_all.py -> Ka1, Ka2の同時フィット、あるいはKb1, Kb2の同時フィット
python3 linefitting_all2.py -> Ka1, Ka2を個別でフィット
python3 linefitting_all3.py -> V KaにTi Kbが隣接しているのでそれらをまとめてフィット
python3 linefitting_all4.py -> Cr KaにV Kbが混ざっているのでそれらをまとめてフィット

なので、使い方としては
rm -f linefitting_result_all.txt
python3 linefitting_all.py
python3 linefitting_all2.py
python3 linefitting_all3.py
python3 linefitting_all4.py
とすればよい。

# 結果のまとめ
さらに
python3 FWHMtrend_all.py
python3 gaintrend_all.py
とすると、linefitting_result_all.txtのデータを読み込んで
FWHMとgain errorのenergy依存性を出力してくれる。
