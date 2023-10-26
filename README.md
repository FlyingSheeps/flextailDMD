# 第28回スカイスポーツシンポジウムの関連資料(まだ未完成です)

## 目的
　このリポジトリは元論文におけるDMD (Dynamic Mode Decomposition)によるモード抽出の具体的な方法，およびそれに利用する胴体弾性の計算方法を補完するために作られました．

## 使い方
　まずこのリポジトリをクローンして，その中に入ります．
```
git clone https://github.com/FlyingSheeps/flextailDMD.git
cd flextailDMD
```
　計算コードを実行します．pythonのコードを実行してください，ここではpythonで実行しますが，環境によってはpyやpython3だったりするかもしれません．
```
python tailflexDMD.py
```
　実行した結果はvideoとPhiに出力されます．videoがコマ撮りのたわみ運動，Phiがモードの実部です．

　計算にはmunpy,scipy,control，描画にmatplotlibを使います．適宜pipなどで導入してください．